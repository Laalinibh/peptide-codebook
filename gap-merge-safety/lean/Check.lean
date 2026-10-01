import MassSep
open PepStore Res
-- The empirically optimal alphabet from the density sweep (reach <= 2.7 A)
def optAlphabet : List Res := [A, S, C, T, V, P, L, N, D]
-- The same alphabet with G and I added (no mass constraint applied)
def unsafeAlphabet : List Res := [A, S, C, T, V, P, L, N, D, G, I]

#eval checkMSReadable orbitrap optAlphabet
#eval checkMSReadable orbitrap unsafeAlphabet
-- machine-checked: the optimal alphabet is MS-readable
theorem opt_readable : MSReadable orbitrap optAlphabet :=
  checkMSReadable_sound orbitrap optAlphabet (by native_decide)
