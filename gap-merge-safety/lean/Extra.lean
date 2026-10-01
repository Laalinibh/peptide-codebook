import MassSep
namespace PepStore
open Res


/-! ## 0. Membership inversion for the enumerations. -/

theorem of_mem_allPairs {A : List Res} {p : Res × Res} (h : p ∈ allPairs A) :
    p.1 ∈ A ∧ p.2 ∈ A := by
  unfold allPairs at h
  rw [List.mem_flatMap] at h
  obtain ⟨x, hx, h2⟩ := h
  rw [List.mem_map] at h2
  obtain ⟨y, hy, heq⟩ := h2
  subst heq
  exact ⟨hx, hy⟩

theorem of_mem_allTriples {A : List Res} {t : Res × Res × Res}
    (h : t ∈ allTriples A) : t.1 ∈ A ∧ t.2.1 ∈ A ∧ t.2.2 ∈ A := by
  unfold allTriples at h
  rw [List.mem_flatMap] at h
  obtain ⟨x, hx, h2⟩ := h
  rw [List.mem_flatMap] at h2
  obtain ⟨y, hy, h3⟩ := h2
  rw [List.mem_map] at h3
  obtain ⟨z, hz, heq⟩ := h3
  subst heq
  exact ⟨hx, hy, hz⟩

/-! ## 1. Completeness: rejection implies the predicate genuinely fails. -/

theorem checkSeparated_complete (I : Instrument) (A : List Res) :
    Separated I A → checkSeparated I A = none := by
  intro h
  unfold checkSeparated
  cases hf : (allPairs A).find? (badPair I) with
  | none => rfl
  | some p =>
    exfalso
    have hmem := List.find?_some hf
    have hin := of_mem_allPairs (List.mem_of_find?_eq_some hf)
    unfold badPair at hmem
    simp only [Bool.and_eq_true, bne_iff_ne, decide_eq_true_eq] at hmem
    exact absurd (h p.1 hin.1 p.2 hin.2 hmem.1) (Nat.not_lt.mpr hmem.2)

theorem checkGapSafe_complete (I : Instrument) (A : List Res) :
    GapSafe I A → checkGapSafe I A = none := by
  intro h
  unfold checkGapSafe
  cases hf : (allTriples A).find? (badTriple I) with
  | none => rfl
  | some t =>
    exfalso
    have hmem := List.find?_some hf
    have hin := of_mem_allTriples (List.mem_of_find?_eq_some hf)
    unfold badTriple at hmem
    simp only [decide_eq_true_eq] at hmem
    exact absurd (h t.1 hin.1 t.2.1 hin.2.1 t.2.2 hin.2.2) (Nat.not_lt.mpr hmem)

/-- The checker decides `MSReadable` exactly: no false accepts, no false rejects. -/
theorem checkMSReadable_iff (I : Instrument) (A : List Res) :
    checkMSReadable I A = none ↔ MSReadable I A := by
  constructor
  · exact checkMSReadable_sound I A
  · intro ⟨hs, hg⟩
    unfold checkMSReadable
    rw [checkSeparated_complete I A hs]
    exact checkGapSafe_complete I A hg

/-! ## 2. Hereditary: readability is preserved under taking sub-alphabets.

The paper uses this to justify enumerating maximal readable alphabets via
minimal forbidden configurations rather than searching 2^20 subsets. It was
previously asserted in a test; here it is proven. -/

theorem Separated_mono {I : Instrument} {A B : List Res}
    (hsub : ∀ x ∈ B, x ∈ A) (h : Separated I A) : Separated I B :=
  fun a ha b hb hne => h a (hsub a ha) b (hsub b hb) hne

theorem GapSafe_mono {I : Instrument} {A B : List Res}
    (hsub : ∀ x ∈ B, x ∈ A) (h : GapSafe I A) : GapSafe I B :=
  fun a ha b hb c hc => h a (hsub a ha) b (hsub b hb) c (hsub c hc)

theorem MSReadable_mono {I : Instrument} {A B : List Res}
    (hsub : ∀ x ∈ B, x ∈ A) (h : MSReadable I A) : MSReadable I B :=
  ⟨Separated_mono hsub h.1, GapSafe_mono hsub h.2⟩

/-! ## 3. Tolerance monotonicity.

Verifying at a wider tolerance is strictly stronger. This is the formal content
of the effective-tolerance rule: an alphabet checked at tau_eff >= tau is also
valid at tau, so a codebook verified for a coarse instrument remains valid on a
finer one -- but not conversely. -/

theorem Separated_tol_mono {I J : Instrument} {A : List Res}
    (hle : J.tol ≤ I.tol) (h : Separated I A) : Separated J A :=
  fun a ha b hb hne => Nat.lt_of_le_of_lt hle (h a ha b hb hne)

theorem GapSafe_tol_mono {I J : Instrument} {A : List Res}
    (hle : J.tol ≤ I.tol) (h : GapSafe I A) : GapSafe J A :=
  fun a ha b hb c hc => Nat.lt_of_le_of_lt hle (h a ha b hb c hc)

theorem MSReadable_tol_mono {I J : Instrument} {A : List Res}
    (hle : J.tol ≤ I.tol) (h : MSReadable I A) : MSReadable J A :=
  ⟨Separated_tol_mono hle h.1, GapSafe_tol_mono hle h.2⟩

/-! ## 4. Scaffold context.

The paper's corrected design consequence: gap-merge safety must hold over the
union of the encoding alphabet and the fixed scaffold residues, because the
merging pair need not consist of two encoding residues. Checking the alphabet
alone is insufficient; this makes the requirement precise. -/

def ReadableInContext (I : Instrument) (alphabet scaffold : List Res) : Prop :=
  MSReadable I (alphabet ++ scaffold)

/-- Context-readability is strictly stronger than alphabet-readability. -/
theorem context_implies_alphabet (I : Instrument) (alphabet scaffold : List Res) :
    ReadableInContext I alphabet scaffold → MSReadable I alphabet :=
  fun h => MSReadable_mono (fun _ hx => List.mem_append.mpr (Or.inl hx)) h

/-- ...and the converse fails: ASCTVPLND is readable alone but not with a
    glycine-bearing scaffold, since A+G aliases Q. -/
example : MSReadable orbitrap [A,S,C,T,V,P,L,N,D] := by
  exact checkMSReadable_sound _ _ (by native_decide)

example : ¬ ReadableInContext orbitrap [A,S,C,T,V,P,L,N,D] [G,Q] := by
  intro h
  have := (checkMSReadable_iff orbitrap ([A,S,C,T,V,P,L,N,D] ++ [G,Q])).mpr h
  revert this
  native_decide

end PepStore
