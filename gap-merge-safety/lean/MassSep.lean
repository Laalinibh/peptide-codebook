namespace PepStore

abbrev Mass := Nat

inductive Res
  | G | A | S | P | V | T | C | L | I | N
  | D | Q | K | E | M | H | F | R | Y | W
  deriving DecidableEq, Repr, Inhabited

open Res

def allRes : List Res := [G,A,S,P,V,T,C,L,I,N,D,Q,K,E,M,H,F,R,Y,W]

/-- Monoisotopic residue masses in units of 1e-5 Da. Exact integers: no Float. -/
def mass : Res → Mass
  | G =>  5702146 | A =>  7103711 | S =>  8703203 | P =>  9705276
  | V =>  9906841 | T => 10104768 | C => 10300919 | L => 11308406
  | I => 11308406 | N => 11404293 | D => 11502694 | Q => 12805858
  | K => 12809496 | E => 12904259 | M => 13104049 | H => 13705891
  | F => 14706841 | R => 15610111 | Y => 16306333 | W => 18607931

def adist (a b : Mass) : Mass := max a b - min a b

structure Instrument where
  tol : Mass
  deriving Repr

def orbitrap : Instrument := ⟨500⟩
def qtof     : Instrument := ⟨2000⟩
def ionTrap  : Instrument := ⟨50000⟩

/-- Pairwise separation: distinct letters must be resolvable. -/
def Separated (I : Instrument) (A : List Res) : Prop :=
  ∀ a ∈ A, ∀ b ∈ A, a ≠ b → adist (mass a) (mass b) > I.tol

/-- Gap-merge safety: a merged two-residue gap must not alias a single letter. -/
def GapSafe (I : Instrument) (A : List Res) : Prop :=
  ∀ a ∈ A, ∀ b ∈ A, ∀ c ∈ A, adist (mass a + mass b) (mass c) > I.tol

def MSReadable (I : Instrument) (A : List Res) : Prop := Separated I A ∧ GapSafe I A

inductive Violation
  | unresolvable (a b : Res) (gap : Mass)
  | gapMerge     (a b c : Res) (gap : Mass)
  deriving Repr, DecidableEq

def allPairs (A : List Res) : List (Res × Res) :=
  A.flatMap (fun a => A.map (fun b => (a, b)))

def allTriples (A : List Res) : List (Res × Res × Res) :=
  A.flatMap (fun a => A.flatMap (fun b => A.map (fun c => (a, b, c))))

def badPair (I : Instrument) (p : Res × Res) : Bool :=
  p.1 != p.2 && decide (adist (mass p.1) (mass p.2) ≤ I.tol)

def badTriple (I : Instrument) (t : Res × Res × Res) : Bool :=
  decide (adist (mass t.1 + mass t.2.1) (mass t.2.2) ≤ I.tol)

def checkSeparated (I : Instrument) (A : List Res) : Option Violation :=
  match (allPairs A).find? (badPair I) with
  | some (a, b) => some (.unresolvable a b (adist (mass a) (mass b)))
  | none => none

def checkGapSafe (I : Instrument) (A : List Res) : Option Violation :=
  match (allTriples A).find? (badTriple I) with
  | some (a, b, c) => some (.gapMerge a b c (adist (mass a + mass b) (mass c)))
  | none => none

/-- Kernel entry point. `none` = readable; `some v` = actionable certificate. -/
def checkMSReadable (I : Instrument) (A : List Res) : Option Violation :=
  match checkSeparated I A with
  | some v => some v
  | none => checkGapSafe I A

/-! ## Membership lemmas -/

theorem mem_allPairs {A : List Res} {a b : Res} (ha : a ∈ A) (hb : b ∈ A) :
    (a, b) ∈ allPairs A := by
  unfold allPairs
  rw [List.mem_flatMap]
  exact ⟨a, ha, List.mem_map.mpr ⟨b, hb, rfl⟩⟩

theorem mem_allTriples {A : List Res} {a b c : Res}
    (ha : a ∈ A) (hb : b ∈ A) (hc : c ∈ A) : (a, b, c) ∈ allTriples A := by
  unfold allTriples
  rw [List.mem_flatMap]
  refine ⟨a, ha, ?_⟩
  rw [List.mem_flatMap]
  exact ⟨b, hb, List.mem_map.mpr ⟨c, hc, rfl⟩⟩

/-! ## Soundness -/

theorem checkSeparated_sound (I : Instrument) (A : List Res) :
    checkSeparated I A = none → Separated I A := by
  intro h a ha b hb hne
  unfold checkSeparated at h
  have hfind : (allPairs A).find? (badPair I) = none := by
    cases hf : (allPairs A).find? (badPair I) with
    | none => rfl
    | some p => rw [hf] at h; cases p; simp at h
  have := List.find?_eq_none.mp hfind (a, b) (mem_allPairs ha hb)
  unfold badPair at this
  simp only [Bool.and_eq_true, bne_iff_ne, decide_eq_true_eq, not_and] at this
  exact Nat.lt_of_not_le (this hne)

theorem checkGapSafe_sound (I : Instrument) (A : List Res) :
    checkGapSafe I A = none → GapSafe I A := by
  intro h a ha b hb c hc
  unfold checkGapSafe at h
  have hfind : (allTriples A).find? (badTriple I) = none := by
    cases hf : (allTriples A).find? (badTriple I) with
    | none => rfl
    | some t => rw [hf] at h; obtain ⟨x, y, z⟩ := t; simp at h
  have := List.find?_eq_none.mp hfind (a, b, c) (mem_allTriples ha hb hc)
  unfold badTriple at this
  simp only [decide_eq_true_eq] at this
  exact Nat.lt_of_not_le this

/-- Main soundness theorem: acceptance implies the semantic predicate. -/
theorem checkMSReadable_sound (I : Instrument) (A : List Res) :
    checkMSReadable I A = none → MSReadable I A := by
  intro h
  unfold checkMSReadable at h
  cases hs : checkSeparated I A with
  | some v => rw [hs] at h; simp at h
  | none =>
    rw [hs] at h
    exact ⟨checkSeparated_sound I A hs, checkGapSafe_sound I A h⟩

end PepStore
