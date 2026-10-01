"""Tests asserting the claims made in the paper.

Each test corresponds to a stated result. Run with: python -m pytest tests/ -v
"""
import itertools
import numpy as np
import pytest

import resconst as rc
import msmodel as M

TAU = 1000        # 0.01 Da in units of 1e-5 Da
ALL = list("GAVLISTCMNQDEKRHFYWP")
VERIFIED = list("ASCTVPLND")


# ---------------------------------------------------------------- predicates
def separated(A, tol):
    return all(abs(rc.res_mass_i[a] - rc.res_mass_i[b]) > tol
               for a, b in itertools.permutations(A, 2))


def gap_safe(A, tol):
    return all(abs(rc.res_mass_i[a] + rc.res_mass_i[b] - rc.res_mass_i[c]) > tol
               for a in A for b in A for c in A)


def ms_readable(A, tol):
    return separated(A, tol) and gap_safe(A, tol)


# ---------------------------------------------- exhaustive characterization
def test_exactly_one_isobaric_pair():
    """Leu/Ile is the only exact isobar in the proteinogenic alphabet."""
    pairs = [(a, b) for a, b in itertools.combinations(ALL, 2)
             if rc.res_mass_i[a] == rc.res_mass_i[b]]
    assert pairs == [("L", "I")] or pairs == [("I", "L")]


def test_gln_lys_is_next_closest():
    d = sorted((abs(rc.res_mass_i[a] - rc.res_mass_i[b]), a, b)
               for a, b in itertools.combinations(ALL, 2)
               if rc.res_mass_i[a] != rc.res_mass_i[b])
    gap, a, b = d[0]
    assert {a, b} == {"Q", "K"}
    assert gap == 3638          # 0.03638 Da


def test_exactly_two_summed_gap_collisions():
    """Gly+Gly = Asn and Gly+Ala = Gln are the only summed-gap collisions."""
    hits = {(a, b, c) for a in ALL for b in ALL for c in ALL
            if abs(rc.res_mass_i[a] + rc.res_mass_i[b] - rc.res_mass_i[c]) <= TAU}
    unordered = {(frozenset([a, b]), c) for a, b, c in hits}
    assert unordered == {(frozenset(["G"]), "N"), (frozenset(["G", "A"]), "Q")}


def test_collisions_involve_glycine():
    hits = {(a, b, c) for a in ALL for b in ALL for c in ALL
            if abs(rc.res_mass_i[a] + rc.res_mass_i[b] - rc.res_mass_i[c]) <= TAU}
    assert all("G" in (a, b) for a, b, c in hits)


# ---------------------------------------------------------------- hereditary
def test_predicates_are_hereditary():
    """Dropping letters can never create a violation."""
    rng = np.random.default_rng(0)
    for _ in range(200):
        k = int(rng.integers(3, 12))
        A = list(rng.choice(ALL, size=k, replace=False))
        if not ms_readable(A, TAU):
            continue
        B = A[:-1]
        assert ms_readable(B, TAU), f"{A} readable but subset {B} is not"


def test_maximal_alphabet_is_18_and_excludes_glycine():
    for k in range(len(ALL), 1, -1):
        found = [A for A in itertools.combinations(ALL, k)
                 if ms_readable(list(A), TAU)]
        if found:
            assert k == 18
            assert all("G" not in A for A in found)
            break


def test_verified_alphabet_is_readable():
    assert ms_readable(VERIFIED, TAU)


def test_adding_G_or_I_breaks_readability():
    assert not ms_readable(VERIFIED + ["I"], TAU)   # L/I isobar
    assert not ms_readable(VERIFIED + ["G"], TAU)   # G+G = N


# ---------------------------------------------------------- mixed radix
def test_mixed_radix_injective():
    sizes = [3, 5, 2, 7, 4]
    seen = set()
    for w in itertools.product(*[range(s) for s in sizes]):
        n = 0
        for a, s in zip(reversed(w), reversed(sizes)):
            n = n * s + a
        assert n not in seen
        seen.add(n)
    assert len(seen) == int(np.prod(sizes))


# ------------------------------------------------- tolerance must scale
def test_nominal_tolerance_fails_under_mass_error():
    """The central correction: verifying at nominal resolution is unsound
    when measurement error at the largest prefix mass exceeds it."""
    rng = np.random.default_rng(11)
    pep = "".join(rng.choice(VERIFIED) for _ in range(17)) + "R"
    silent, _, _ = M.error_rates(pep, VERIFIED + ["R"], 0.15, TAU,
                                 np.random.default_rng(3), 1000,
                                 ppm=5.0, n_noise=0, both_series=True)
    assert silent > 0.005, "expected silent errors at nominal tolerance"


def test_scaled_tolerance_restores_guarantee():
    rng = np.random.default_rng(11)
    pep = "".join(rng.choice(VERIFIED) for _ in range(17)) + "R"
    Mmax = sum(rc.res_mass_i[a] for a in pep)
    tau_eff = int(3 * 5.0e-6 * Mmax)
    assert ms_readable(VERIFIED, tau_eff)
    silent, _, _ = M.error_rates(pep, VERIFIED + ["R"], 0.15, tau_eff,
                                 np.random.default_rng(3), 1000,
                                 ppm=5.0, n_noise=0, both_series=True)
    assert silent == 0.0


def test_scaled_tolerance_costs_at_most_two_letters():
    def max_size(tol):
        for k in range(len(ALL), 1, -1):
            if any(ms_readable(list(A), tol)
                   for A in itertools.combinations(ALL, k)):
                return k
    assert max_size(TAU) == 18
    assert max_size(2800) == 17
    assert max_size(5600) == 16
