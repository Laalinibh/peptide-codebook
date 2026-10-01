# Supporting Information: code and derived data

Code and derived data accompanying *Non-Isobaric Is Not Enough: Gap-Merge Safety
for Unambiguous Peptide Readout* (L. Bhogadi).

Everything reported in the paper runs on CPU. No structure predictor is executed
anywhere in this archive.

## Layout

    resconst.py        residue constants: exact integer masses (1e-5 Da), atom
                       types, van der Waals radii, side-chain reach, max ASA
    msmodel.py         tandem MS readout: b/y ladders, proline-directed cleavage
                       bias, charge retention, ppm mass error, noise peaks
    kernel.py          validation tiers T0-T4 over a Protein representation,
                       grid and exhaustive clash detection
    sites.py           SASA, encoding-site eligibility, separation-constrained
                       site selection, idealized readout and decoding
    rotamer.py         fixed-backbone side-chain construction (NeRF) and the
                       steric substitution oracle
    calibrate.py       the calibration subsection: measures fragmentation parameters from
                       identified public spectra
    exp_calibrated.py  the calibration subsection: re-runs the decode experiment under
                       idealized, stated and measured parameters
    lean/              machine-checked kernel (Lean 4.15, zero `sorry`)
    tests/             claim-level test suite
    calibration/       outputs of the two scripts above (see below)

## Which file supports which claim

| Paper | Claim | Artifact |
|---|---|---|
| Theorem 1 (decision procedure) | `checkMSReadable I A = none <-> MSReadable I A` | `lean/MassSep.lean` |
| Thm 2 (hereditary) | readability closed under sub-alphabets | `lean/MassSep.lean` |
| Thm 3 (tolerance monotonicity) | verification at wider tau is stronger | `lean/MassSep.lean` |
| Thm 4 (composition) | separated envelopes compose | `lean/Independence.lean` |
| Thm 5 (encoder injectivity) | mixed-radix encoding is injective | `lean/Codebook.lean` |
| Sec 3 | exhaustive alphabet characterization | `tests/test_claims.py` |
| Calibration | fragmentation calibration | `calibrate.py`, `calibration/` |
| Calibration | guarantee survives calibration | `exp_calibrated.py`, `calibration/` |
| Sec 5 | grid vs exhaustive clash detection agree | `kernel.py`, `tests/` |

## Running

Tests (about one second):

    PYTHONPATH=. python -m pytest tests/ -v

Lean (Lean 4.15, no `sorry`):

    cd lean
    lean -o MassSep.olean MassSep.lean
    LEAN_PATH=. lean Check.lean

Re-running the calibration of the calibration subsection requires the two public spectral sets
named below; the derived outputs are included so the section can be checked
without them:

    PYTHONPATH=. python calibrate.py \
        --massivekb-mgf <path>/massivekb_selected_500.mgf \
        --massivekb-peprec <path>/massivekb_selected_500.peprec \
        --casanovo-mgf <path>/sample_preprocessed_spectra.mgf

    PYTHONPATH=. python exp_calibrated.py --trials 3000

## Calibration data provenance

The spectra themselves are not redistributed here. Both sets are public test
data shipped with open-source tools under permissive licences:

- `massivekb_selected_500.mgf` / `.peprec` — a 500-spectrum MassIVE-KB HCD
  subset distributed in the test data of the MS2PIP project.
- `sample_preprocessed_spectra.mgf` — 128 annotated sample spectra distributed
  with the Casanovo project.

Included derived files:

    calibration/per_site_indicators.csv     5,112 rows, one per backbone cleavage
                                            site of the 423 retained identified
                                            peptides: source, peptide, site index,
                                            flanking residues, b/y observation
                                            indicators, C-terminal basicity
    calibration/calibration_output.txt      parameter estimates, per source and
                                            pooled (the calibration tables)
    calibration/exp_calibrated_output.txt   decode experiment under all three
                                            fragmentation processes 

## Panel

The fifteen-chain structural panel of the folded-carrier section consists of
deposited experimental structures, listed by accession in the Supporting Information PDF. The
coordinate files are not redistributed here; they are retrievable from the PDB
by the accessions given. The oracle dataset and learned proposer of the oracle subsection
were computed on an earlier six-structure panel and are carried over unchanged,
as stated in that section.

## Dependencies

`numpy` only, for everything except the Lean development. No GPU.
