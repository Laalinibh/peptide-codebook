"""Calibrate the fragmentation model of msmodel.py against identified public spectra.

Measures, per backbone cleavage site of an identified peptide, whether the
corresponding b and/or y fragment is observed. From those indicators we estimate
the three structural parameters msmodel.py states rather than measures:

  p_base      baseline per-cleavage fragment dropout
  Pro factors enhancement N-terminal to proline, suppression C-terminal to it
  charge      y-ion survival boost when the C-terminal half carries K/R/H

Inputs are two independent public sets of identified spectra (see SOURCES).
No spectrum is used unless its peptide identification is given by the source.
"""
import argparse, math, re, sys
from collections import defaultdict

import numpy as np
import resconst as rc

PROTON = 1.00727646
H2O = 18.0105646

SOURCES = {
    "massivekb": "MassIVE-KB HCD subset distributed with ms2pip (500 identified spectra)",
    "casanovo": "Annotated sample spectra distributed with Casanovo (128 identified spectra)",
}


# ---------------------------------------------------------------- MGF parsing

def parse_mgf(path):
    """Yield (title, charge, peaks) with peaks an ascending m/z array."""
    title, charge, mz = None, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line == "BEGIN IONS":
                title, charge, mz = None, None, []
            elif line == "END IONS":
                if title is not None:
                    yield title, charge, np.array(sorted(mz))
            elif "=" in line:
                k, _, v = line.partition("=")
                if k == "TITLE":
                    title = v
                elif k == "CHARGE":
                    charge = int(re.sub(r"[^0-9]", "", v) or 0)
                elif k == "SEQ":
                    title = title if title is not None else v
            elif line:
                parts = line.split()
                if len(parts) >= 1:
                    try:
                        mz.append(float(parts[0]))
                    except ValueError:
                        pass


def parse_mgf_with_seq(path):
    """Casanovo-style MGF where SEQ= carries the identification."""
    seq, charge, mz = None, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line == "BEGIN IONS":
                seq, charge, mz = None, None, []
            elif line == "END IONS":
                if seq:
                    yield seq, charge, np.array(sorted(mz))
            elif "=" in line:
                k, _, v = line.partition("=")
                if k == "SEQ":
                    seq = v
                elif k == "CHARGE":
                    charge = int(re.sub(r"[^0-9]", "", v) or 0)
            elif line:
                parts = line.split()
                try:
                    mz.append(float(parts[0]))
                except ValueError:
                    pass


def load_peprec(path):
    """spec_id -> (peptide, charge, modifications) from an ms2pip PEPREC file."""
    out = {}
    with open(path) as fh:
        header = fh.readline().split()
        idx = {k: i for i, k in enumerate(header)}
        for line in fh:
            f = line.split()
            if len(f) < len(header):
                continue
            out[f[idx["spec_id"]]] = (
                f[idx["peptide"]], int(f[idx["charge"]]), f[idx["modifications"]])
    return out


# ------------------------------------------------------------- ion bookkeeping

def observed(peaks, target, tol_da):
    """Is there a peak within tol_da of target m/z?"""
    if peaks.size == 0:
        return False
    i = np.searchsorted(peaks, target)
    for j in (i - 1, i):
        if 0 <= j < peaks.size and abs(peaks[j] - target) <= tol_da:
            return True
    return False


def site_indicators(seq, peaks, precursor_charge, tol_da, max_frag_charge=2):
    """For each cleavage site k (between residue k-1 and k), return a record of
    whether b_k and y_{n-k} were observed at any charge up to max_frag_charge."""
    n = len(seq)
    masses = [rc.res_mass[a] for a in seq]
    pre = np.cumsum(masses)                      # pre[k-1] = sum of first k residues
    total = pre[-1] + H2O
    zmax = max(1, min(max_frag_charge, (precursor_charge or 2) - 1))

    recs = []
    for k in range(1, n):
        b_neutral = pre[k - 1]                   # b ion neutral mass
        y_neutral = total - pre[k - 1]           # complementary y ion neutral mass
        b_seen = any(observed(peaks, (b_neutral + z * PROTON) / z, tol_da)
                     for z in range(1, zmax + 1))
        y_seen = any(observed(peaks, (y_neutral + z * PROTON) / z, tol_da)
                     for z in range(1, zmax + 1))
        recs.append(dict(
            k=k,
            left=seq[k - 1],
            right=seq[k],
            b=b_seen,
            y=y_seen,
            either=b_seen or y_seen,
            cterm_basic=any(a in "KRH" for a in seq[k:]),
        ))
    return recs


# --------------------------------------------------------------- aggregation

def wilson(k, n, z=1.96):
    """Wilson score interval, so small strata report honest uncertainty."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def summarize(recs, label):
    n = len(recs)
    if n == 0:
        print(f"[{label}] no usable cleavage sites")
        return None

    either = sum(r["either"] for r in recs)
    b_only = sum(r["b"] for r in recs)
    y_only = sum(r["y"] for r in recs)

    # Baseline strata exclude every proline-adjacent site, so the proline
    # factors are measured against a clean reference.
    base = [r for r in recs if r["left"] != "P" and r["right"] != "P"]
    n_base = len(base)
    keep_base = sum(r["either"] for r in base) / n_base if n_base else float("nan")

    npro = [r for r in recs if r["right"] == "P"]          # N-terminal to Pro
    cpro = [r for r in recs if r["left"] == "P"]           # C-terminal to Pro
    keep_npro = (sum(r["either"] for r in npro) / len(npro)) if npro else float("nan")
    keep_cpro = (sum(r["either"] for r in cpro) / len(cpro)) if cpro else float("nan")

    yb = [r for r in base if r["cterm_basic"]]
    yn = [r for r in base if not r["cterm_basic"]]
    y_basic = (sum(r["y"] for r in yb) / len(yb)) if yb else float("nan")
    y_nonbasic = (sum(r["y"] for r in yn) / len(yn)) if yn else float("nan")

    out = dict(
        label=label, sites=n, peptides=None,
        keep_any=either / n, p_drop=1 - either / n,
        b_rate=b_only / n, y_rate=y_only / n,
        keep_base=keep_base, n_base=n_base,
        keep_npro=keep_npro, n_npro=len(npro),
        keep_cpro=keep_cpro, n_cpro=len(cpro),
        pro_enh=keep_npro / keep_base if n_base and npro else float("nan"),
        pro_sup=keep_cpro / keep_base if n_base and cpro else float("nan"),
        y_basic=y_basic, n_ybasic=len(yb),
        y_nonbasic=y_nonbasic, n_ynonbasic=len(yn),
        charge_ratio=y_basic / y_nonbasic if yn and yb else float("nan"),
    )
    lo, hi = wilson(either, n)
    out["keep_any_ci"] = (lo, hi)
    return out


def report(s):
    print(f"\n=== {s['label']} ===")
    print(f"cleavage sites analysed        {s['sites']}")
    print(f"fragment observed (b or y)     {s['keep_any']:.3f}  "
          f"[95% CI {s['keep_any_ci'][0]:.3f}, {s['keep_any_ci'][1]:.3f}]")
    print(f"  -> empirical p_drop          {s['p_drop']:.3f}")
    print(f"b-ion observation rate         {s['b_rate']:.3f}")
    print(f"y-ion observation rate         {s['y_rate']:.3f}")
    print(f"non-proline baseline keep      {s['keep_base']:.3f}   (n={s['n_base']})")
    print(f"keep N-terminal to Pro         {s['keep_npro']:.3f}   (n={s['n_npro']})"
          f"   factor {s['pro_enh']:.2f}  [model 1.60]")
    print(f"keep C-terminal to Pro         {s['keep_cpro']:.3f}   (n={s['n_cpro']})"
          f"   factor {s['pro_sup']:.2f}  [model 0.35]")
    print(f"y rate, C-term basic           {s['y_basic']:.3f}   (n={s['n_ybasic']})")
    print(f"y rate, C-term not basic       {s['y_nonbasic']:.3f}   (n={s['n_ynonbasic']})")
    print(f"  -> charge-retention ratio    {s['charge_ratio']:.2f}"
          f"   [model 1.15/0.85 = 1.35]")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--massivekb-mgf")
    ap.add_argument("--massivekb-peprec")
    ap.add_argument("--casanovo-mgf")
    ap.add_argument("--tol-da", type=float, default=0.02,
                    help="fragment match tolerance (ms2pip config states 0.02 Da)")
    ap.add_argument("--min-len", type=int, default=6)
    args = ap.parse_args()

    all_recs = {}

    if args.massivekb_mgf and args.massivekb_peprec:
        ids = load_peprec(args.massivekb_peprec)
        recs, used, skipped_mod = [], 0, 0
        for title, charge, peaks in parse_mgf(args.massivekb_mgf):
            hit = ids.get(title)
            if hit is None:
                continue
            pep, z, mods = hit
            # Restrict to unmodified peptides over the standard alphabet: a
            # modification shifts fragment masses and would be scored as a miss.
            if mods not in ("-", "", "None"):
                skipped_mod += 1
                continue
            if len(pep) < args.min_len or any(a not in rc.res_mass for a in pep):
                continue
            recs += site_indicators(pep, peaks, z or charge, args.tol_da)
            used += 1
        print(f"massivekb: {used} identified unmodified peptides used, "
              f"{skipped_mod} modified peptides skipped")
        all_recs["massivekb"] = recs

    if args.casanovo_mgf:
        recs, used = [], 0
        for seq, charge, peaks in parse_mgf_with_seq(args.casanovo_mgf):
            pep = re.sub(r"[^A-Z]", "", seq)
            if len(pep) < args.min_len or any(a not in rc.res_mass for a in pep):
                continue
            if pep != seq:           # carried a modification annotation
                continue
            recs += site_indicators(pep, peaks, charge, args.tol_da)
            used += 1
        print(f"casanovo: {used} identified unmodified peptides used")
        all_recs["casanovo"] = recs

    summaries = []
    for k, recs in all_recs.items():
        s = summarize(recs, f"{k} -- {SOURCES[k]}")
        if s:
            report(s)
            summaries.append(s)

    if len(all_recs) > 1:
        pooled = [r for recs in all_recs.values() for r in recs]
        s = summarize(pooled, "pooled")
        report(s)
        summaries.append(s)

    return summaries


if __name__ == "__main__":
    main()


def dump_indicators(out_csv, **kw):
    """Write the per-site indicator table used by the paper's calibration."""
    import csv
    rows = []
    ids = load_peprec(kw["massivekb_peprec"])
    for title, charge, peaks in parse_mgf(kw["massivekb_mgf"]):
        hit = ids.get(title)
        if not hit:
            continue
        pep, z, mods = hit
        if mods not in ("-", "", "None") or len(pep) < 6:
            continue
        if any(a not in rc.res_mass for a in pep):
            continue
        for r in site_indicators(pep, peaks, z or charge, kw.get("tol_da", 0.02)):
            rows.append(dict(source="massivekb", peptide=pep, **r))
    for seq, charge, peaks in parse_mgf_with_seq(kw["casanovo_mgf"]):
        pep = re.sub(r"[^A-Z]", "", seq)
        if pep != seq or len(pep) < 6 or any(a not in rc.res_mass for a in pep):
            continue
        for r in site_indicators(pep, peaks, charge, kw.get("tol_da", 0.02)):
            rows.append(dict(source="casanovo", peptide=pep, **r))
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return len(rows)
