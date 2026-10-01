"""Encoding-site selection and the MS readback simulation."""
import numpy as np
import resconst as rc
import kernel as K

# ---------------------------------------------------------------- SASA
def _sphere_points(n=200):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(theta) * np.sin(phi),
                 np.sin(theta) * np.sin(phi), np.cos(phi)]


def sasa(prot, probe=1.4, n_points=200):
    """Shrake-Rupley solvent accessible surface area, per residue (A^2)."""
    idx = np.argwhere(prot.atom_mask > 0.5)
    xyz = prot.atom_positions[idx[:, 0], idx[:, 1]].astype(np.float64)
    rad = np.array([rc.vdw.get(rc.atom_types[j][0], 1.7) for j in idx[:, 1]]) + probe
    sph = _sphere_points(n_points)
    n_res = prot.aatype.shape[0]
    out = np.zeros(n_res)
    # neighbour lists
    for a in range(len(xyz)):
        d = np.linalg.norm(xyz - xyz[a], axis=1)
        nb = np.where((d < rad + rad[a]) & (d > 0))[0]
        pts = xyz[a] + rad[a] * sph
        acc = np.ones(len(pts), dtype=bool)
        for b in nb:
            acc &= np.linalg.norm(pts - xyz[b], axis=1) >= rad[b]
        out[idx[a, 0]] += 4 * np.pi * rad[a] ** 2 * acc.mean()
    return out


def relative_sasa(prot):
    a = sasa(prot)
    seq = K.seq_of(prot)
    return np.array([a[i] / rc.max_asa[seq[i]] for i in range(len(seq))])


# ---------------------------------------------------------------- sites
def candidate_sites(prot, rsa_min=0.35):
    """Positions eligible to carry information before the separation filter."""
    rsa = relative_sasa(prot)
    seq = K.seq_of(prot)
    ss = K.disulfides(prot)
    cys_bonded = {i for pair in ss for i in pair}
    cands = []
    for i, aa in enumerate(seq):
        if rsa[i] < rsa_min:
            continue
        if aa in ("G", "P"):
            continue           # backbone-defining, substitution alters the fold
        if i in cys_bonded:
            continue           # structural disulfide
        if prot.atom_mask[i, rc.atom_order["CB"]] < 0.5:
            continue
        cands.append(i)
    return cands, rsa


def select_sites(prot, cands, alphabet, delta=3.0, eps=0.05, greedy_order="rsa"):
    """Maximal set of sites whose envelopes are pairwise separated, i.e. that
    satisfy the hypothesis of the Lean composition theorem.

    Envelope radius r_i = max side-chain reach over the alphabet.
    Required separation: d(CB_i, CB_j) > r_i + r_j + delta + 2*eps.
    """
    rmax = max(rc.sidechain_reach[a] for a in alphabet)
    cb = rc.atom_order["CB"]
    pos = {i: prot.atom_positions[i, cb].astype(np.float64) for i in cands}
    rsa = relative_sasa(prot)
    order = sorted(cands, key=lambda i: -rsa[i]) if greedy_order == "rsa" else list(cands)
    chosen = []
    need = 2 * rmax + delta + 2 * eps
    for i in order:
        if all(np.linalg.norm(pos[i] - pos[j]) > need for j in chosen):
            chosen.append(i)
    return sorted(chosen), need


def site_alphabet(prot, i, base_alphabet, rsa):
    """Per-site alphabet after the structural filter. Proxy-based: no folding
    model is available offline, so this uses exposure-compatible hydropathy and
    envelope fit rather than a learned stability predictor."""
    out = []
    for a in base_alphabet:
        if rc.sidechain_reach[a] > 6.5:
            continue
        # a highly exposed site should not be given a strongly buried-type
        # residue; a partially exposed one tolerates more
        if rsa[i] > 0.5 and rc.kd[a] > 3.0:
            continue
        out.append(a)
    return out


# ---------------------------------------------------------------- MS model
def mass_int(seq):
    return sum(rc.res_mass_i[a] for a in seq)


def simulate_readout(seq, p_drop, rng):
    """Return the observed set of prefix masses (b-ion ladder) after dropout.
    Termini are always observed."""
    n = len(seq)
    prefix = np.cumsum([rc.res_mass_i[a] for a in seq])
    observed = [0]
    for k in range(1, n):
        if rng.random() >= p_drop:
            observed.append(int(prefix[k - 1]))
    observed.append(int(prefix[-1]))
    return sorted(set(observed))


def decode(observed, alphabet, tol, site_positions=None):
    """Read residues from consecutive observed prefix-mass gaps.

    Returns (calls, n_ambiguous) where calls maps a gap to either a single
    residue letter, None if the gap matches no single letter (a DETECTED
    dropout), or a set if several letters match (detected ambiguity).
    """
    calls = []
    for a, b in zip(observed, observed[1:]):
        gap = b - a
        matches = [r for r in alphabet if abs(rc.res_mass_i[r] - gap) <= tol]
        if len(matches) == 1:
            calls.append((a, b, matches[0]))
        elif len(matches) == 0:
            calls.append((a, b, None))          # detected: no letter fits
        else:
            calls.append((a, b, set(matches)))  # detected: ambiguous
    return calls


def error_rates(seq, alphabet, p_drop, tol, rng, trials=2000):
    """Silent (undetected, wrong) vs detected error rate per residue call."""
    silent = 0
    detected = 0
    total = 0
    prefix = np.cumsum([rc.res_mass_i[a] for a in seq])
    truth = {0: None}
    for k in range(len(seq)):
        truth[int(prefix[k])] = seq[k]
    for _ in range(trials):
        obs = simulate_readout(seq, p_drop, rng)
        for a, b, call in decode(obs, alphabet, tol):
            total += 1
            spans_one = (b in truth and a in truth and
                         _residues_between(prefix, a, b) == 1)
            if spans_one:
                true_res = seq[_index_of_prefix(prefix, b)]
                if isinstance(call, str):
                    if call != true_res:
                        silent += 1
                else:
                    detected += 1
            else:
                # gap spans >1 residue: a correct decoder must NOT emit a
                # confident single-residue call here
                if isinstance(call, str):
                    silent += 1     # aliased: reads one wrong residue silently
                else:
                    detected += 1
    return silent / total, detected / total, total


def _residues_between(prefix, a, b):
    ia = 0 if a == 0 else int(np.searchsorted(prefix, a)) + 1
    ib = int(np.searchsorted(prefix, b)) + 1
    return ib - ia


def _index_of_prefix(prefix, b):
    return int(np.searchsorted(prefix, b))
