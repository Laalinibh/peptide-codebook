# Gap-Merge Safety: Unambiguous Peptide Readout

Code, machine-checked proofs and derived data for the paper

> **Non-Isobaric Is Not Enough: Gap-Merge Safety for Unambiguous Peptide Readout**
> Laalini Bhogadi · Preprint on ChemRxiv · under review at the *Journal of Proteome Research*

## The problem

When peptide sequences carry information, the standard safeguard is to encode
with **non-isobaric** amino acids, so no two residues share a mass. The paper
shows this is not enough.

Tandem MS reads a sequence as the mass *differences* between consecutive
fragment ions. If a single b/y ion is missing, two adjacent residues merge into
one observed gap. When that summed gap lands on the mass of a third residue,
the decoder does not see a dropout. It emits one confident, **wrong** residue.

Non-isobaricity cannot prevent this, because the colliding quantity is a sum,
not a residue.

## The result

The missing condition is **gap-merge safety**:

    |m(a) + m(b) − m(c)| > τ    for every a, b, c in the alphabet

- It is a decidable predicate, with a **machine-checked soundness and completeness proof in Lean 4** (zero `sorry`).
- An exhaustive search of the proteinogenic alphabet finds that, beyond the Leu/Ile isobar, exactly two summed-gap collisions exist: **Gly + Gly = Asn** and **Gly + Ala = Gln**.
- Gap-merge-safe alphabets produce **zero silent errors** at every simulated fragment-dropout rate. Every failure becomes a detectable erasure.
- The guarantee survives calibration against **5,112 cleavage sites from 423 identified public spectra**.
- The verification tolerance must scale with the largest prefix mass, not nominal instrument resolution.
- Applied to carriers that must stay folded proteins, verified codebooks carry **0.31–0.79 bits per residue**.

## Repository layout

| Path | What it is |
|---|---|
| `resconst.py` | Residue constants: exact integer masses (1e-5 Da), atom types, van der Waals radii, side-chain reach, max ASA |
| `msmodel.py` | Tandem MS readout model: b/y ladders, proline-directed cleavage bias, charge retention, ppm mass error, noise peaks |
| `kernel.py` | Structural validation tiers T0–T4, with grid and exhaustive clash detection |
| `sites.py` | SASA, encoding-site eligibility, separation-constrained site selection, readout and decoding |
| `rotamer.py` | Fixed-backbone side-chain construction (NeRF) and the steric substitution oracle |
| `calibrate.py` | Measures fragmentation parameters from identified public spectra |
| `exp_calibrated.py` | Re-runs the decode experiment under idealized, stated and measured parameters |
| `lean/` | Machine-checked kernel (Lean 4.15, zero `sorry`) |
| `tests/` | Claim-level test suite |
| `calibration/` | Derived outputs of the calibration scripts |

## Which file backs which claim

| Paper | Claim | Artifact |
|---|---|---|
| Theorem 1 | `checkMSReadable I A = none ↔ MSReadable I A` | `lean/MassSep.lean`, `lean/Extra.lean` |
| Theorem 2 | Readability is closed under sub-alphabets | `lean/Extra.lean` |
| Theorem 3 | Verifying at a wider tolerance is stronger | `lean/Extra.lean` |
| Theorem 4 | Separated encoding sites compose | `lean/Independence.lean` |
| Theorem 5 | Mixed-radix encoding is injective | `lean/Codebook.lean` |
| Section 3 | Exhaustive alphabet characterization | `tests/test_claims.py` |
| Calibration | Fragmentation parameters and the calibrated decode experiment | `calibrate.py`, `exp_calibrated.py`, `calibration/` |
| Section 5 | Grid and exhaustive clash detection agree | `kernel.py`, `tests/` |

## Quick start

Everything runs on a CPU. The Python code needs only `numpy`.

```bash
pip install -r requirements.txt
PYTHONPATH=. python -m pytest tests/ -v      # 12 tests, about a second
```

### Lean proofs

```bash
cd lean
lean -o MassSep.olean MassSep.lean
LEAN_PATH=. lean Check.lean                  # checks the working alphabet is MS-readable
```

### Re-running the calibration

The spectra are public test data shipped with two open-source projects and are
not redistributed here. The derived outputs in `calibration/` let you check the
results without them.

```bash
PYTHONPATH=. python calibrate.py \
    --massivekb-mgf <path>/massivekb_selected_500.mgf \
    --massivekb-peprec <path>/massivekb_selected_500.peprec \
    --casanovo-mgf <path>/sample_preprocessed_spectra.mgf

PYTHONPATH=. python exp_calibrated.py --trials 3000
```

- `massivekb_selected_500.mgf` / `.peprec`: a 500-spectrum MassIVE-KB HCD subset from the [MS2PIP](https://github.com/compomics/ms2pip) test data.
- `sample_preprocessed_spectra.mgf`: annotated sample spectra from [Casanovo](https://github.com/Noble-Lab/casanovo).

## A note on tolerance

Both predicates take a tolerance parameter. Set it to the measurement
uncertainty at the **largest** prefix mass:

    τ_eff = k · ppm · 1e-6 · M_max

Verifying at nominal instrument resolution is unsound whenever the mass error
at the largest fragment exceeds it. `tests/test_claims.py` checks both
directions.

## Structural panel

The fifteen-chain panel used for folded carriers consists of deposited
experimental structures, listed by PDB accession in the paper's Supporting
Information. Coordinate files are not included; fetch them from the
[RCSB PDB](https://www.rcsb.org/). Residue conventions follow AlphaFold/OpenFold.

## Citation

If you use this work, please cite the preprint (see `CITATION.cff`).

## License

MIT. See `LICENSE`.
