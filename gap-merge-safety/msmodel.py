"""Realistic tandem MS readout model.

Replaces the idealized model (b ions only, independent uniform dropout, exact
masses, no noise) with:
  - both b and y ion series
  - proline-directed cleavage bias (enhanced N-terminal, suppressed C-terminal)
  - basic-residue charge retention affecting y-ion survival
  - Gaussian mass measurement error at a stated ppm
  - spurious noise peaks that can split a genuine gap into two false gaps
"""
import numpy as np
import resconst as rc

PROTON = 1.00728
H2O = 18.010565

# Integer units of 1e-5 Da throughout, matching the Lean kernel.
PROTON_I = int(round(PROTON * 1e5))
H2O_I = int(round(H2O * 1e5))


def prefix_masses(seq):
    return np.cumsum([rc.res_mass_i[a] for a in seq])


def cleavage_prob(seq, k, p_base):
    """Survival probability of the fragment cleaved between residue k-1 and k.

    Proline effects: cleavage N-terminal to Pro is strongly enhanced (so that
    fragment is well observed); cleavage C-terminal to Pro is suppressed.
    """
    p_keep = 1.0 - p_base
    left = seq[k - 1]
    right = seq[k] if k < len(seq) else None
    if right == "P":          # N-terminal to proline: enhanced
        p_keep = min(1.0, p_keep * 1.6)
    if left == "P":           # C-terminal to proline: suppressed
        p_keep = p_keep * 0.35
    return float(np.clip(p_keep, 0.0, 1.0))


def simulate_spectrum(seq, p_base, rng, ppm=5.0, n_noise=0, both_series=True):
    """Return observed prefix masses (integer 1e-5 Da), reconstructed from the
    b and y ladders, with measurement error and optional noise peaks."""
    n = len(seq)
    pre = prefix_masses(seq)
    total = int(pre[-1])
    obs = set([0, total])

    for k in range(1, n):
        pk = cleavage_prob(seq, k, p_base)
        b_seen = rng.random() < pk
        # y ion for the same cleavage; charge retention favours the side with a
        # basic residue, so y survival is boosted when the C-term half has one
        cterm_basic = any(a in "KRH" for a in seq[k:])
        py = pk * (1.15 if cterm_basic else 0.85)
        y_seen = both_series and (rng.random() < min(py, 1.0))
        if b_seen or y_seen:
            obs.add(int(pre[k - 1]))

    out = []
    for m in obs:
        if m in (0, total):
            out.append(m)
        else:
            sigma = ppm * 1e-6 * m
            out.append(int(round(m + rng.normal(0, sigma))))

    for _ in range(n_noise):
        out.append(int(rng.integers(1, max(total, 2))))
    return sorted(set(out))


def decode(observed, alphabet, tol):
    calls = []
    for a, b in zip(observed, observed[1:]):
        gap = b - a
        m = [r for r in alphabet if abs(rc.res_mass_i[r] - gap) <= tol]
        calls.append((a, b, m[0] if len(m) == 1 else (None if not m else set(m))))
    return calls


def error_rates(seq, alphabet, p_base, tol, rng, trials=2000,
                ppm=5.0, n_noise=0, both_series=True):
    """Silent (confident and wrong) vs detected error rate per residue call."""
    pre = prefix_masses(seq)
    true_at = {int(pre[k]): seq[k] for k in range(len(seq))}
    true_at[0] = None
    silent = detected = total = 0
    for _ in range(trials):
        obs = simulate_spectrum(seq, p_base, rng, ppm, n_noise, both_series)
        for a, b, call in decode(obs, alphabet, tol):
            total += 1
            # how many true residues does this gap actually span?
            ia = int(np.searchsorted(pre, a + tol, side="right"))
            ib = int(np.searchsorted(pre, b + tol, side="right"))
            spans = ib - ia
            if spans == 1:
                # find the true residue for this single-residue gap
                idx = ib - 1
                true_res = seq[idx] if 0 <= idx < len(seq) else None
                if isinstance(call, str):
                    if call != true_res:
                        silent += 1
                else:
                    detected += 1
            else:
                # gap spans 0 or >1 residues: any confident call is wrong
                if isinstance(call, str):
                    silent += 1
                else:
                    detected += 1
    return silent / total, detected / total, total
