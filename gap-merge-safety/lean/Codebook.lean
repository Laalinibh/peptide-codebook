/-
  Codebook.lean -- mixed-radix encoding and its injectivity.

  The paper states that injectivity of the encoder "is exactly decodability".
  That makes it a proof obligation rather than a remark, and it is discharged
  here rather than asserted by round-trip test.
-/
namespace PepStore

/-- Per-site alphabet sizes, innermost site first. -/
abbrev Radices := List Nat

/-- A message: one index per site. -/
abbrev Digits := List Nat

/-- Digits are in range when each is below its site's alphabet size. -/
inductive InRange : Radices → Digits → Prop
  | nil : InRange [] []
  | cons {r : Nat} {rs : Radices} {d : Nat} {ds : Digits} :
      d < r → InRange rs ds → InRange (r :: rs) (d :: ds)

/-- Mixed-radix encoding: `n = d0 + r0*(d1 + r1*(...))`. -/
def encode : Radices → Digits → Nat
  | [], _ => 0
  | _, [] => 0
  | r :: rs, d :: ds => d + r * encode rs ds

/-- The encoding of in-range digits is bounded by the product of the radices,
    so distinct sites cannot interfere. This is the key lemma. -/
theorem encode_lt : ∀ {rs : Radices} {ds : Digits},
    InRange rs ds → encode rs ds < rs.foldr (· * ·) 1
  | [], [], _ => by simp [encode]
  | r :: rs, d :: ds, h => by
    cases h with
    | cons hd htail =>
      have ih := encode_lt htail
      simp only [encode, List.foldr]
      calc d + r * encode rs ds
          < r + r * encode rs ds := by omega
        _ = r * (encode rs ds + 1) := by
              rw [Nat.mul_succ]; exact Nat.add_comm r (r * encode rs ds)
        _ ≤ r * rs.foldr (· * ·) 1 := by
              exact Nat.mul_le_mul_left r (by omega)

/-- **Decodability.** The encoder is injective on in-range digit strings.
    Equivalently: no two distinct messages produce the same readout. -/
theorem encode_injective : ∀ {rs : Radices} {ds es : Digits},
    InRange rs ds → InRange rs es → encode rs ds = encode rs es → ds = es
  | [], [], [], _, _, _ => rfl
  | r :: rs, d :: ds, e :: es, hd, he, heq => by
    cases hd with
    | cons hdr hdt =>
      cases he with
      | cons her het =>
        simp only [encode] at heq
        have hlt1 := encode_lt hdt
        have hlt2 := encode_lt het
        -- d and e are both < r, so they are the residues mod r
        have hdeq : d = e := by
          have : d % r = e % r := by
            have h1 : (d + r * encode rs ds) % r = d % r := by
              simp [Nat.add_mul_mod_self_left]
            have h2 : (e + r * encode rs es) % r = e % r := by
              simp [Nat.add_mul_mod_self_left]
            rw [← h1, ← h2, heq]
          rwa [Nat.mod_eq_of_lt hdr, Nat.mod_eq_of_lt her] at this
        subst hdeq
        have hr : 0 < r := Nat.lt_of_le_of_lt (Nat.zero_le d) hdr
        have htail : encode rs ds = encode rs es := by
          have : r * encode rs ds = r * encode rs es := by omega
          exact Nat.eq_of_mul_eq_mul_left hr this
        rw [encode_injective hdt het htail]

/-- Capacity: the number of distinct representable messages. -/
def capacity (rs : Radices) : Nat := rs.foldr (· * ·) 1

example : capacity [8,8,8,8,8,8,8,8] = 16777216 := by native_decide
example : capacity [4,9,8,12] = 3456 := by native_decide

end PepStore
