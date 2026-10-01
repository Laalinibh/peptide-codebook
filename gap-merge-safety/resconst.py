"""Minimal residue_constants, matching the AlphaFold/OpenFold conventions used by
the Protein dataclass. Written from standard chemistry tables."""
import numpy as np

atom_types = [
    "N", "CA", "C", "CB", "O", "CG", "CG1", "CG2", "OG", "OG1", "SG", "CD",
    "CD1", "CD2", "ND1", "ND2", "OD1", "OD2", "SD", "CE", "CE1", "CE2", "CE3",
    "NE", "NE1", "NE2", "OE1", "OE2", "CH2", "NH1", "NH2", "OH", "CZ", "CZ2",
    "CZ3", "NZ", "OXT",
]
atom_order = {a: i for i, a in enumerate(atom_types)}
atom_type_num = len(atom_types)  # 37

restypes = list("ARNDCQEGHILKMFPSTWYV")
restype_order = {r: i for i, r in enumerate(restypes)}
restype_num = len(restypes)  # 20

restype_1to3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}
restype_3to1 = {v: k for k, v in restype_1to3.items()}

# Heavy atoms present in each residue type.
residue_atoms = {
    "ALA": ["N", "CA", "C", "O", "CB"],
    "ARG": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["N", "CA", "C", "O", "CB", "CG", "ND2", "OD1"],
    "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
    "CYS": ["N", "CA", "C", "O", "CB", "SG"],
    "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE2", "OE1"],
    "GLU": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],
    "GLY": ["N", "CA", "C", "O"],
    "HIS": ["N", "CA", "C", "O", "CB", "CG", "CD2", "ND1", "CE1", "NE2"],
    "ILE": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"],
    "LEU": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],
    "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
    "MET": ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"],
    "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD"],
    "SER": ["N", "CA", "C", "O", "CB", "OG"],
    "THR": ["N", "CA", "C", "O", "CB", "CG2", "OG1"],
    "TRP": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE2", "CE3",
            "NE1", "CH2", "CZ2", "CZ3"],
    "TYR": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2",
            "CZ", "OH"],
    "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2"],
}

STANDARD_ATOM_MASK = np.zeros((21, atom_type_num), dtype=np.float32)
for _r in restypes:
    for _a in residue_atoms[restype_1to3[_r]]:
        STANDARD_ATOM_MASK[restype_order[_r], atom_order[_a]] = 1.0

# van der Waals radii (Angstrom), by element.
vdw = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}

# Monoisotopic residue masses (Da).
res_mass = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276,
    "V": 99.06841, "T": 101.04768, "C": 103.00919, "L": 113.08406,
    "I": 113.08406, "N": 114.04293, "D": 115.02694, "Q": 128.05858,
    "K": 128.09496, "E": 129.04259, "M": 131.04049, "H": 137.05891,
    "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}
# Integer masses in units of 1e-5 Da -- exact, matching the Lean kernel.
res_mass_i = {k: int(round(v * 1e5)) for k, v in res_mass.items()}

# Kyte-Doolittle hydropathy.
kd = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
    "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
    "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
    "Y": -1.3, "V": 4.2,
}

# Max side-chain extent from CB (Angstrom), used to size encoding envelopes.
# Computed as the largest CB-to-sidechain-atom distance in ideal geometry.
sidechain_reach = {
    "G": 0.0, "A": 0.0, "S": 1.42, "C": 1.81, "T": 1.52, "V": 1.53,
    "P": 2.30, "L": 2.61, "I": 2.60, "N": 2.48, "D": 2.47, "M": 3.90,
    "E": 3.77, "Q": 3.78, "K": 5.05, "H": 3.55, "R": 6.20, "F": 4.31,
    "Y": 5.70, "W": 5.30,
}

# Max solvent-accessible area of residue X in Gly-X-Gly, for relative SASA.
max_asa = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}
