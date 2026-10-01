"""Explicit side-chain construction and the physical substitution oracle.

Replaces the hydropathy/reach proxy with fixed-backbone rotamer repacking:
for a candidate substitution we build the side chain from ideal internal
coordinates over a discrete chi-rotamer set and accept the residue only if some
rotamer is clash-free against the (fixed) remainder of the structure.

Limitations, stated up front: rigid backbone, neighbouring side chains held
fixed, coarse rotamer set (chi1/chi2 on 60-degree grids), no electrostatics or
solvation energy, no minimisation. This is a steric feasibility oracle, not a
free-energy calculation.
"""
import numpy as np
import resconst as rc
import kernel as K

# ---------------------------------------------------------------- NeRF
def place_atom(a, b, c, bond, angle_deg, dihedral_deg):
    """Place atom D given A-B-C and internal coordinates (NeRF)."""
    ang = np.radians(angle_deg)
    tor = np.radians(dihedral_deg)
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d2 = np.array([-bond * np.cos(ang),
                   bond * np.sin(ang) * np.cos(tor),
                   bond * np.sin(ang) * np.sin(tor)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


# Side-chain build trees beyond CB.
# (atom, parent_a, parent_b, parent_c, bond, angle, chi_spec)
# chi_spec: ("chi", k) uses chi_k; ("chi_off", k, off) uses chi_k + off;
#           ("fix", v) uses a fixed dihedral.
BUILD = {
    "A": [],
    "S": [("OG", "N", "CA", "CB", 1.417, 110.8, ("chi", 1))],
    "C": [("SG", "N", "CA", "CB", 1.808, 114.0, ("chi", 1))],
    "T": [("OG1", "N", "CA", "CB", 1.433, 109.6, ("chi", 1)),
          ("CG2", "N", "CA", "CB", 1.521, 110.5, ("chi_off", 1, -120.0))],
    "V": [("CG1", "N", "CA", "CB", 1.521, 110.5, ("chi", 1)),
          ("CG2", "N", "CA", "CB", 1.521, 110.5, ("chi_off", 1, 120.0))],
    "I": [("CG1", "N", "CA", "CB", 1.530, 110.4, ("chi", 1)),
          ("CG2", "N", "CA", "CB", 1.521, 110.5, ("chi_off", 1, -122.0)),
          ("CD1", "CA", "CB", "CG1", 1.521, 113.8, ("chi", 2))],
    "L": [("CG", "N", "CA", "CB", 1.530, 116.3, ("chi", 1)),
          ("CD1", "CA", "CB", "CG", 1.521, 110.7, ("chi", 2)),
          ("CD2", "CA", "CB", "CG", 1.521, 110.7, ("chi_off", 2, 122.0))],
    "N": [("CG", "N", "CA", "CB", 1.516, 112.6, ("chi", 1)),
          ("OD1", "CA", "CB", "CG", 1.231, 120.8, ("chi", 2)),
          ("ND2", "CA", "CB", "CG", 1.328, 116.4, ("chi_off", 2, 180.0))],
    "D": [("CG", "N", "CA", "CB", 1.516, 112.6, ("chi", 1)),
          ("OD1", "CA", "CB", "CG", 1.249, 118.4, ("chi", 2)),
          ("OD2", "CA", "CB", "CG", 1.249, 118.4, ("chi_off", 2, 180.0))],
    "M": [("CG", "N", "CA", "CB", 1.520, 114.1, ("chi", 1)),
          ("SD", "CA", "CB", "CG", 1.803, 112.7, ("chi", 2)),
          ("CE", "CB", "CG", "SD", 1.791, 100.9, ("fix", -60.0))],
    "E": [("CG", "N", "CA", "CB", 1.520, 114.1, ("chi", 1)),
          ("CD", "CA", "CB", "CG", 1.516, 112.6, ("chi", 2)),
          ("OE1", "CB", "CG", "CD", 1.249, 118.4, ("fix", 0.0)),
          ("OE2", "CB", "CG", "CD", 1.249, 118.4, ("fix", 180.0))],
    "Q": [("CG", "N", "CA", "CB", 1.520, 114.1, ("chi", 1)),
          ("CD", "CA", "CB", "CG", 1.516, 112.6, ("chi", 2)),
          ("OE1", "CB", "CG", "CD", 1.231, 120.8, ("fix", 0.0)),
          ("NE2", "CB", "CG", "CD", 1.328, 116.4, ("fix", 180.0))],
    "K": [("CG", "N", "CA", "CB", 1.520, 114.1, ("chi", 1)),
          ("CD", "CA", "CB", "CG", 1.520, 111.3, ("chi", 2)),
          ("CE", "CB", "CG", "CD", 1.520, 111.3, ("fix", 180.0)),
          ("NZ", "CG", "CD", "CE", 1.489, 111.9, ("fix", 180.0))],
}
BUILDABLE = sorted(BUILD.keys())

CHI_GRID = [-60.0, 60.0, 180.0]


def _rotamers(res):
    """Discrete rotamer set for a residue: chi1 (and chi2 where used)."""
    spec = BUILD[res]
    uses2 = any(s[6][0] in ("chi", "chi_off") and s[6][1] == 2 for s in spec)
    if not spec:
        return [{}]
    if uses2:
        return [{1: c1, 2: c2} for c1 in CHI_GRID for c2 in CHI_GRID]
    return [{1: c1} for c1 in CHI_GRID]


def build_sidechain(prot, i, res, chi):
    """Build heavy side-chain atoms of `res` at site i on the fixed backbone."""
    coords = {}
    for nm in ("N", "CA", "CB"):
        j = rc.atom_order[nm]
        if prot.atom_mask[i, j] < 0.5:
            return None
        coords[nm] = prot.atom_positions[i, j].astype(np.float64)
    out = []
    for (atom, pa, pb, pc, bond, angle, spec) in BUILD[res]:
        if spec[0] == "chi":
            dih = chi[spec[1]]
        elif spec[0] == "chi_off":
            dih = chi[spec[1]] + spec[2]
        else:
            dih = spec[1]
        if pa not in coords or pb not in coords or pc not in coords:
            return None
        p = place_atom(coords[pa], coords[pb], coords[pc], bond, angle, dih)
        coords[atom] = p
        out.append((atom, p))
    return out


# ---------------------------------------------------------------- oracle
def _environment(prot, exclude):
    """Atoms of the structure excluding residue indices in `exclude`."""
    idx = np.argwhere(prot.atom_mask > 0.5)
    keep = ~np.isin(idx[:, 0], list(exclude))
    idx = idx[keep]
    xyz = prot.atom_positions[idx[:, 0], idx[:, 1]].astype(np.float64)
    rad = np.array([rc.vdw.get(rc.atom_types[j][0], 1.7) for j in idx[:, 1]])
    return xyz, rad, idx


def substitution_ok(prot, i, res, overlap=0.4, env=None):
    """True if some rotamer of `res` at site i is clash-free on the fixed
    backbone. Excludes residue i and its sequence neighbours from the
    environment (their covalent geometry is not a clash)."""
    if res == "G":
        return True
    if res not in BUILD:
        return None            # not buildable; caller must handle
    if env is None:
        env = _environment(prot, {i - 1, i, i + 1})
    xyz, rad, _ = env
    if not BUILD[res]:         # alanine: CB only, always present
        return True
    for chi in _rotamers(res):
        built = build_sidechain(prot, i, res, chi)
        if built is None:
            return False
        ok = True
        for atom, p in built:
            r = rc.vdw.get(atom[0], 1.7)
            d = np.linalg.norm(xyz - p, axis=1)
            if np.any(d < rad + r - overlap):
                ok = False
                break
        if ok:
            return True
    return False


def site_alphabet_physical(prot, i, base_alphabet, env=None):
    """Per-site alphabet from the rotamer oracle."""
    if env is None:
        env = _environment(prot, {i - 1, i, i + 1})
    return [a for a in base_alphabet
            if a in BUILD and substitution_ok(prot, i, a, env=env)]


# ---------------------------------------------------------------- features
def site_features(prot, i, rsa, phi, psi):
    """Cheap descriptors available without running the oracle."""
    ca = prot.atom_positions[:, rc.atom_order["CA"]].astype(np.float64)
    d = np.linalg.norm(ca - ca[i], axis=1)
    seq = K.seq_of(prot)
    nb8 = int(((d < 8.0).sum()) - 1)
    nb10 = int(((d < 10.0).sum()) - 1)
    nb12 = int(((d < 12.0).sum()) - 1)
    near = np.argsort(d)[1:9]
    return dict(
        rsa=float(rsa[i]),
        nb8=nb8, nb10=nb10, nb12=nb12,
        phi=float(phi[i]) if not np.isnan(phi[i]) else 0.0,
        psi=float(psi[i]) if not np.isnan(psi[i]) else 0.0,
        has_phi=float(not np.isnan(phi[i])),
        kd_native=rc.kd[seq[i]],
        kd_env=float(np.mean([rc.kd[seq[j]] for j in near])),
        depth=float(np.linalg.norm(ca[i] - ca.mean(0))),
        n_res=float(len(seq)),
    )
