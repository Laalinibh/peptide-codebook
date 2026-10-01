"""Re-run the readout experiments under a fragmentation process calibrated on
identified public spectra, rather than the stated parameters of the paper.

The paper's Tables 1-2 use the idealized process of sites.py: b ions only,
independent uniform dropout, exact masses. This script keeps the exact-mass
decoding and the silent/detected accounting of sites.py unchanged, and replaces
only the *dropout process* with one whose parameters were measured by
calibrate.py over identified spectra:

  - both b and y series contribute a prefix observation
  - cleavage N-terminal to Pro enhanced, C-terminal to Pro suppressed
  - y-ion survival depends on whether the C-terminal half carries K/R/H

The question is whether gap-merge safety still eliminates silent errors when the
dropout process is sequence-dependent and empirically parameterized, rather than
independent and uniform. Exact masses are retained so that this isolates the
fragmentation process; measurement error is the separate axis of Table 3.
"""
import argparse, random
import numpy as np

import resconst as rc
import sites as S

TAU = 1000                              # 0.01 Da in units of 1e-5 Da
VERIFIED = list("ASCTVPLND")            # gap-merge safe
NONISO = list("ASCTVPLNDG")             # non-isobaric only; violates Definition 2

# Measured by calibrate.py over 5,112 cleavage sites from 423 identified
# unmodified peptides across two independent public sets.
MEASURED = dict(pro_enh=1.10, pro_sup=0.79, y_basic=1.00, y_nonbasic=0.55,
                both_series=True)
# As stated in msmodel.py.
STATED = dict(pro_enh=1.60, pro_sup=0.35, y_basic=1.15, y_nonbasic=0.85,
              both_series=True)
# The idealized process of sites.py, for reference.
IDEAL = dict(pro_enh=1.0, pro_sup=1.0, y_basic=0.0, y_nonbasic=0.0,
             both_series=False)


def simulate_readout(seq, p_drop, rng, par):
    """Observed prefix masses under a parameterized fragmentation process.

    Exact masses; termini always observed. With both_series=False and unit
    proline factors this reduces exactly to sites.simulate_readout.
    """
    n = len(seq)
    prefix = np.cumsum([rc.res_mass_i[a] for a in seq])
    obs = [0]
    for k in range(1, n):
        p_keep = 1.0 - p_drop
        if seq[k] == "P":
            p_keep = min(1.0, p_keep * par["pro_enh"])
        if seq[k - 1] == "P":
            p_keep = p_keep * par["pro_sup"]
        b_seen = rng.random() < p_keep
        y_seen = False
        if par["both_series"]:
            basic = any(a in "KRH" for a in seq[k:])
            py = p_keep * (par["y_basic"] if basic else par["y_nonbasic"])
            y_seen = rng.random() < min(py, 1.0)
        if b_seen or y_seen:
            obs.append(int(prefix[k - 1]))
    obs.append(int(prefix[-1]))
    return sorted(set(obs))


def error_rates(seq, alphabet, p_drop, tol, rng, par, trials=2000):
    """sites.error_rates with the fragmentation process swapped."""
    silent = detected = total = 0
    prefix = np.cumsum([rc.res_mass_i[a] for a in seq])
    truth = {0: None}
    for k in range(len(seq)):
        truth[int(prefix[k])] = seq[k]
    for _ in range(trials):
        obs = simulate_readout(seq, p_drop, rng, par)
        for a, b, call in S.decode(obs, alphabet, tol):
            total += 1
            spans_one = (b in truth and a in truth and
                         S._residues_between(prefix, a, b) == 1)
            if spans_one:
                true_res = seq[S._index_of_prefix(prefix, b)]
                if isinstance(call, str):
                    if call != true_res:
                        silent += 1
                else:
                    detected += 1
            else:
                if isinstance(call, str):
                    silent += 1
                else:
                    detected += 1
    return silent / total, detected / total, total


def rand_pep(pool, n, r):
    return "".join(r.choice(pool) for _ in range(n))


def table1(par, trials, seed=0, length=18):
    """Silent vs detected for the verified and non-isobaric-only alphabets.

    Peptides are drawn over the decoder's own alphabet, matching the paper.
    """
    rows = []
    for p in (0.0, 0.1, 0.2, 0.3, 0.4):
        row = [p]
        for alpha in (VERIFIED, NONISO):
            rng = np.random.default_rng(seed)
            r = random.Random(seed)
            s_t = d_t = n_t = 0
            for _ in range(trials // 50):
                seq = rand_pep(alpha, length, r)
                s, d, n = error_rates(seq, alpha, p, TAU, rng, par, trials=50)
                s_t += s * n; d_t += d * n; n_t += n
            row += [s_t / n_t, d_t / n_t]
        rows.append(row)
    return rows


def table2(par, trials, seed=0):
    """The two summed-gap collisions, planted, as in the paper's motif table."""
    cases = [
        ("GG", list("ASCTVPLNDG"), list("ASCTVPLNDG"), "contains N"),
        ("GA", list("ASCTVPLNDG"), list("ASCTVPLNDG"), "no Q"),
        ("GA", list("ASCTVPLNDGQ"), list("ASCTVPLNDGQ"), "contains Q"),
        # G in the scaffold, removed from the decoder alphabet only
        ("GG", list("ASCTVPLND"), list("ASCTVPLNDG"), "G removed from alphabet only"),
        # G absent from scaffold and alphabet
        ("AA", list("ASCTVPLND"), list("ASCTVPLND"), "G absent from both"),
    ]
    out = []
    for motif, decoder_alpha, scaffold_pool, note in cases:
        row = [motif, "".join(decoder_alpha), note]
        for p in (0.1, 0.3):
            rng = np.random.default_rng(seed)
            r = random.Random(seed)
            s_t = n_t = 0
            for _ in range(trials // 50):
                seq = list(rand_pep(scaffold_pool, 18, r))
                pos = r.randrange(0, 16)
                seq[pos:pos + 2] = list(motif)
                s, d, n = error_rates("".join(seq), decoder_alpha, p, TAU,
                                      rng, par, trials=50)
                s_t += s * n; n_t += n
            row.append(s_t / n_t)
        out.append(row)
    return out


def show(name, par, trials):
    print(f"\n################ {name}")
    print(f"   {par}")
    print("\n-- Table 1: silent vs detected, 18-mers, tau = 0.01 Da, exact masses")
    print(f"{'p_drop':>7} {'V silent':>10} {'V detect':>10} "
          f"{'NI silent':>10} {'NI detect':>10}")
    for p, vs, vd, ns, nd in table1(par, trials):
        print(f"{p:7.2f} {vs:10.4f} {vd:10.4f} {ns:10.4f} {nd:10.4f}")

    print("\n-- Table 2: planted motifs, silent-error rate")
    print(f"{'motif':>6} {'decoder':>12} {'p=0.1':>8} {'p=0.3':>8}   note")
    for motif, alpha, note, s1, s3 in table2(par, trials):
        print(f"{motif:>6} {alpha:>12} {s1:8.4f} {s3:8.4f}   {note}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=3000)
    a = ap.parse_args()
    show("IDEAL (reproduces the paper's stated process)", IDEAL, a.trials)
    show("STATED (msmodel.py parameters)", STATED, a.trials)
    show("MEASURED (calibrated on identified public spectra)", MEASURED, a.trials)
