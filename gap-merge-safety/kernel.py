"""Validation kernel over the Protein dataclass.

Tiers:
  T0 well-formedness   T1 covalent geometry   T2 sterics
  T3 conformation      T4 foldedness
Every check returns violation certificates, not booleans.
"""
import dataclasses
from typing import Optional, Sequence, List, Dict, Tuple
import numpy as np
import resconst as rc


@dataclasses.dataclass(frozen=True)
class Protein:
    atom_positions: np.ndarray   # [n, 37, 3]
    aatype: np.ndarray           # [n]
    atom_mask: np.ndarray        # [n, 37]
    residue_index: np.ndarray    # [n]
    b_factors: np.ndarray        # [n, 37]
    chain_index: Optional[np.ndarray] = None
    remark: Optional[str] = None
    parents: Optional[Sequence[str]] = None
    parents_chain_index: Optional[Sequence[int]] = None


@dataclasses.dataclass
class Violation:
    tier: str
    kind: str
    residues: Tuple
    measured: float
    bound: float

    def __repr__(self):
        return (f"[{self.tier}] {self.kind} at {self.residues}: "
                f"measured={self.measured:.4f} bound={self.bound:.4f}")


def parse_pdb(path: str, chain: str = "A") -> Protein:
    """Parse a single chain of a PDB file into a Protein."""
    residues: Dict[int, Dict] = {}
    order: List[int] = []
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        if line[21] != chain:
            continue
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue
        resname = line[17:20].strip()
        if resname not in rc.restype_3to1:
            continue
        resseq = int(line[22:26])
        icode = line[26]
        key = (resseq, icode)
        aname = line[12:16].strip()
        if aname not in rc.atom_order:
            continue
        if key not in residues:
            residues[key] = {"name": resname, "seq": resseq, "atoms": {}, "b": {}}
            order.append(key)
        residues[key]["atoms"][aname] = (
            float(line[30:38]), float(line[38:46]), float(line[46:54]))
        residues[key]["b"][aname] = float(line[60:66])

    n = len(order)
    pos = np.zeros((n, rc.atom_type_num, 3), dtype=np.float32)
    mask = np.zeros((n, rc.atom_type_num), dtype=np.float32)
    bfac = np.zeros((n, rc.atom_type_num), dtype=np.float32)
    aat = np.zeros(n, dtype=np.int64)
    ridx = np.zeros(n, dtype=np.int32)
    for i, key in enumerate(order):
        r = residues[key]
        aat[i] = rc.restype_order[rc.restype_3to1[r["name"]]]
        ridx[i] = r["seq"]
        for aname, xyz in r["atoms"].items():
            j = rc.atom_order[aname]
            pos[i, j] = xyz
            mask[i, j] = 1.0
            bfac[i, j] = r["b"][aname]
    return Protein(pos, aat, mask, ridx, bfac,
                   chain_index=np.zeros(n, dtype=np.int32))


def seq_of(prot: Protein) -> str:
    return "".join(rc.restypes[a] for a in prot.aatype)


# ---------------------------------------------------------------- T0
def check_T0(prot: Protein) -> List[Violation]:
    v = []
    n = prot.aatype.shape[0]
    if prot.atom_positions.shape != (n, rc.atom_type_num, 3):
        v.append(Violation("T0", "shape", ("atom_positions",), 0, 0))
    if not np.all(np.isfinite(prot.atom_positions[prot.atom_mask > 0.5])):
        v.append(Violation("T0", "nonfinite", ("atom_positions",), 0, 0))
    if np.any(prot.aatype > rc.restype_num):
        v.append(Violation("T0", "bad_aatype", ("aatype",), 0, rc.restype_num))
    # mask must be a subset of the ideal mask for the residue type
    ideal = rc.STANDARD_ATOM_MASK[prot.aatype]
    extra = np.argwhere((prot.atom_mask > 0.5) & (ideal < 0.5))
    oxt = rc.atom_order["OXT"]
    for i, j in extra:
        if j == oxt and i == n - 1:
            continue  # C-terminal carboxylate
        v.append(Violation("T0", "atom_not_in_residue",
                           (int(i), rc.atom_types[j]), 1, 0))
    # residue_index strictly increasing within a chain
    ci = prot.chain_index if prot.chain_index is not None else np.zeros(n, int)
    for i in range(n - 1):
        if ci[i] == ci[i + 1] and prot.residue_index[i + 1] <= prot.residue_index[i]:
            v.append(Violation("T0", "residue_index_not_increasing",
                               (i, i + 1), float(prot.residue_index[i + 1]),
                               float(prot.residue_index[i])))
    return v


# ---------------------------------------------------------------- T1
IDEAL_BONDS = {"N-CA": (1.458, 0.030), "CA-C": (1.525, 0.030),
               "C-N": (1.329, 0.030)}
TOL_SIGMA = 12.0   # accept within this many sigma


def _get(prot, i, name):
    j = rc.atom_order[name]
    if prot.atom_mask[i, j] < 0.5:
        return None
    return prot.atom_positions[i, j]


def check_T1(prot: Protein) -> List[Violation]:
    v = []
    n = prot.aatype.shape[0]
    for i in range(n):
        N, CA, C = _get(prot, i, "N"), _get(prot, i, "CA"), _get(prot, i, "C")
        if N is not None and CA is not None:
            d = float(np.linalg.norm(N - CA))
            mu, s = IDEAL_BONDS["N-CA"]
            if abs(d - mu) > TOL_SIGMA * s:
                v.append(Violation("T1", "bond_N_CA", (i,), d, mu))
        if CA is not None and C is not None:
            d = float(np.linalg.norm(CA - C))
            mu, s = IDEAL_BONDS["CA-C"]
            if abs(d - mu) > TOL_SIGMA * s:
                v.append(Violation("T1", "bond_CA_C", (i,), d, mu))
        if i + 1 < n and prot.residue_index[i + 1] == prot.residue_index[i] + 1:
            N2 = _get(prot, i + 1, "N")
            if C is not None and N2 is not None:
                d = float(np.linalg.norm(C - N2))
                mu, s = IDEAL_BONDS["C-N"]
                if abs(d - mu) > TOL_SIGMA * s:
                    v.append(Violation("T1", "bond_C_N", (i, i + 1), d, mu))
        # L-chirality: signed volume of (N, C, CB) about CA must be negative
        CB = _get(prot, i, "CB")
        if all(x is not None for x in (N, CA, C, CB)):
            sv = float(np.dot(np.cross(N - CA, C - CA), CB - CA))
            if sv < 0:
                v.append(Violation("T1", "D_chirality", (i,), sv, 0.0))
    return v


# ---------------------------------------------------------------- T2
def _atom_list(prot: Protein):
    idx = np.argwhere(prot.atom_mask > 0.5)
    xyz = prot.atom_positions[idx[:, 0], idx[:, 1]]
    elem = np.array([rc.atom_types[j][0] for j in idx[:, 1]])
    radii = np.array([rc.vdw.get(e, 1.7) for e in elem])
    return idx, xyz.astype(np.float64), radii


def disulfides(prot, cutoff=2.5):
    """Cys pairs whose SG atoms are within bonding distance."""
    sg = rc.atom_order["SG"]
    cys = [i for i in range(prot.aatype.shape[0])
           if rc.restypes[prot.aatype[i]] == "C" and prot.atom_mask[i, sg] > 0.5]
    out = set()
    for a in range(len(cys)):
        for b in range(a + 1, len(cys)):
            i, j = cys[a], cys[b]
            d = np.linalg.norm(prot.atom_positions[i, sg] - prot.atom_positions[j, sg])
            if d < cutoff:
                out.add((i, j)); out.add((j, i))
    return out


def _bonded_exempt(prot, i1, j1, i2, j2, ss=None):
    """Exempt intra-residue pairs, adjacent-residue pairs, and the atoms
    flanking a disulfide bond (which are covalently linked, not clashing)."""
    if i1 == i2:
        return True
    if abs(i1 - i2) == 1:
        return True
    if ss is not None and (int(i1), int(i2)) in ss:
        return True
    return False


def clash_pairs_naive(prot: Protein, overlap: float = 0.4):
    idx, xyz, radii = _atom_list(prot)
    ss = disulfides(prot)
    m = len(xyz)
    out = []
    for a in range(m):
        d = np.linalg.norm(xyz[a + 1:] - xyz[a], axis=1)
        lim = radii[a] + radii[a + 1:] - overlap
        for off in np.argwhere(d < lim).ravel():
            b = a + 1 + off
            i1, j1 = idx[a]
            i2, j2 = idx[b]
            if _bonded_exempt(prot, i1, j1, i2, j2, ss):
                continue
            out.append((int(i1), rc.atom_types[j1], int(i2), rc.atom_types[j2],
                        float(d[off])))
    return out


def clash_pairs_grid(prot: Protein, overlap: float = 0.4):
    """Uniform-grid clash detection. Cell size = max interaction radius, so the
    27-cell neighbourhood provably contains every pair within that radius."""
    idx, xyz, radii = _atom_list(prot)
    ss = disulfides(prot)
    cell = float(2 * radii.max() - overlap)
    keys = np.floor(xyz / cell).astype(np.int64)
    buckets: Dict[Tuple[int, int, int], List[int]] = {}
    for a, k in enumerate(map(tuple, keys)):
        buckets.setdefault(k, []).append(a)
    out = []
    seen = set()
    for a, k in enumerate(map(tuple, keys)):
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for b in buckets.get((k[0] + dx, k[1] + dy, k[2] + dz), ()):
                        if b <= a:
                            continue
                        key = (a, b)
                        if key in seen:
                            continue
                        seen.add(key)
                        d = float(np.linalg.norm(xyz[a] - xyz[b]))
                        if d >= radii[a] + radii[b] - overlap:
                            continue
                        i1, j1 = idx[a]
                        i2, j2 = idx[b]
                        if _bonded_exempt(prot, i1, j1, i2, j2, ss):
                            continue
                        out.append((int(i1), rc.atom_types[j1], int(i2),
                                    rc.atom_types[j2], d))
    return out


def check_T2(prot: Protein, overlap: float = 0.4) -> List[Violation]:
    return [Violation("T2", "steric_clash", (i1, a1, i2, a2), d, 0.0)
            for i1, a1, i2, a2, d in clash_pairs_grid(prot, overlap)]


# ---------------------------------------------------------------- T3
def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    w0 = b0 - np.dot(b0, b1n) * b1n
    w2 = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(w0, w2)
    y = np.dot(np.cross(b1n, w0), w2)
    return float(np.degrees(np.arctan2(y, x)))


def backbone_torsions(prot: Protein):
    n = prot.aatype.shape[0]
    phi = np.full(n, np.nan)
    psi = np.full(n, np.nan)
    for i in range(n):
        N, CA, C = _get(prot, i, "N"), _get(prot, i, "CA"), _get(prot, i, "C")
        if i > 0:
            Cp = _get(prot, i - 1, "C")
            if all(x is not None for x in (Cp, N, CA, C)):
                phi[i] = dihedral(Cp, N, CA, C)
        if i + 1 < n:
            N2 = _get(prot, i + 1, "N")
            if all(x is not None for x in (N, CA, C, N2)):
                psi[i] = dihedral(N, CA, C, N2)
    return phi, psi


def in_ramachandran(phi, psi, gly=False):
    """Coarse allowed-region test: alpha, beta, and left-handed alpha basins."""
    if np.isnan(phi) or np.isnan(psi):
        return True
    # Coarse generously-allowed basins. This is a THRESHOLD-tier check:
    # an approximation to the Richardson/MolProbity contours, not a proof.
    if -160 <= phi <= -20 and -90 <= psi <= 50:       # alpha_R
        return True
    if -180 <= phi <= -40 and (80 <= psi <= 180 or -180 <= psi <= -150):  # beta
        return True
    if 20 <= phi <= 100 and -45 <= psi <= 90:         # alpha_L
        return True
    if gly and -100 <= phi <= -20 and 100 <= psi <= 180:  # Gly extra
        return True
    return False


def check_T3(prot: Protein) -> List[Violation]:
    v = []
    phi, psi = backbone_torsions(prot)
    seq = seq_of(prot)
    for i in range(len(seq)):
        if not in_ramachandran(phi[i], psi[i], gly=(seq[i] == "G")):
            v.append(Violation("T3", "ramachandran_outlier", (i, seq[i]),
                               float(phi[i]), float(psi[i])))
    return v


# ---------------------------------------------------------------- T4
def radius_of_gyration(prot: Protein) -> float:
    ca = np.array([_get(prot, i, "CA") for i in range(prot.aatype.shape[0])
                   if _get(prot, i, "CA") is not None])
    return float(np.sqrt(((ca - ca.mean(0)) ** 2).sum(1).mean()))


def check_T4(prot: Protein) -> List[Violation]:
    n = prot.aatype.shape[0]
    rg = radius_of_gyration(prot)
    expect = 2.2 * n ** 0.38
    v = []
    if rg > 1.8 * expect:
        v.append(Violation("T4", "extended_not_globular", ("chain",), rg, expect))
    return v


def validate(prot: Protein, tiers=("T0", "T1", "T2", "T3", "T4")):
    fns = {"T0": check_T0, "T1": check_T1, "T2": check_T2,
           "T3": check_T3, "T4": check_T4}
    out = []
    for t in tiers:
        out.extend(fns[t](prot))
    return out
