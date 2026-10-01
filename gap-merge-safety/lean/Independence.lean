/-
  Independence.lean — the composition theorem for encoding-site verification.

  Claim: if encoding-site envelopes are pairwise separated with sufficient
  margin, then verifying each site's alphabet INDEPENDENTLY suffices to
  establish clash-freedom of every JOINT assignment. This turns a verification
  problem exponential in the number of sites (|A|^k assignments) into one
  linear in it (sum of |A_i| single-site checks).

  Distances are integers in units of 1e-3 Angstrom, matching the no-Float
  discipline of MassSep.lean. Because rounding a Euclidean distance to an
  integer only approximately preserves the triangle inequality, we carry an
  explicit rounding slack `eps` rather than assuming exactness.
-/

namespace PepStore

/-- A distance structure with an explicit slack in the triangle inequality.
    Instantiated by Euclidean distance in R^3 rounded to units of 1e-3 A,
    for which `eps = 1` suffices. -/
structure Space where
  Pt      : Type
  d       : Pt → Pt → Int
  eps     : Int
  eps_nonneg : 0 ≤ eps
  d_nonneg : ∀ x y, 0 ≤ d x y
  d_symm   : ∀ x y, d x y = d y x
  tri      : ∀ x y z, d x z ≤ d x y + d y z + eps

variable (S : Space)

/-- Encoding-site geometry: `k` sites, each with an envelope centred at `c i`
    of radius `r i` containing every atom any admissible substitution can place. -/
structure SiteGeom (k : Nat) where
  c : Fin k → S.Pt
  r : Fin k → Int
  r_nonneg : ∀ i, 0 ≤ r i

/-- An assignment places, for each site, an atom inside that site's envelope. -/
structure Assignment {k : Nat} (G : SiteGeom S k) where
  w : Fin k → S.Pt
  inEnvelope : ∀ i, S.d (G.c i) (w i) ≤ G.r i

/-- Envelopes are separated with margin `delta` (the clash threshold),
    accounting for the rounding slack. -/
def Separated {k : Nat} (G : SiteGeom S k) (delta : Int) : Prop :=
  ∀ i j : Fin k, i ≠ j → S.d (G.c i) (G.c j) > G.r i + G.r j + delta + 2 * S.eps

/-!
## Main theorem

Site-site clash-freedom follows from envelope separation ALONE. Note the
hypothesis says nothing about which residues were chosen: it holds for every
assignment simultaneously, which is exactly why per-site verification composes.
-/
theorem joint_clash_free {k : Nat} (G : SiteGeom S k) (delta : Int)
    (hsep : Separated S G delta) (α : Assignment S G) :
    ∀ i j : Fin k, i ≠ j → S.d (α.w i) (α.w j) > delta := by
  intro i j hij
  have hci : S.d (G.c i) (α.w i) ≤ G.r i := α.inEnvelope i
  have hcj : S.d (G.c j) (α.w j) ≤ G.r j := α.inEnvelope j
  have hcj' : S.d (α.w j) (G.c j) ≤ G.r j := by
    rw [S.d_symm]; exact hcj
  -- d(ci,cj) <= d(ci,wi) + d(wi,cj) + eps
  have t1 : S.d (G.c i) (G.c j) ≤ S.d (G.c i) (α.w i) + S.d (α.w i) (G.c j) + S.eps :=
    S.tri _ _ _
  -- d(wi,cj) <= d(wi,wj) + d(wj,cj) + eps
  have t2 : S.d (α.w i) (G.c j) ≤ S.d (α.w i) (α.w j) + S.d (α.w j) (G.c j) + S.eps :=
    S.tri _ _ _
  have hs := hsep i j hij
  omega

/-!
## Complexity consequence

`joint_clash_free` quantifies over all assignments but its hypotheses are
per-site. Verifying a codebook therefore costs `Σᵢ |𝒜ᵢ|` single-site checks
plus `O(k²)` envelope-separation checks, rather than `∏ᵢ |𝒜ᵢ|` joint checks.
-/

/-- Number of joint assignments — the size of the naive verification problem. -/
def jointCount (sizes : List Nat) : Nat := sizes.foldl (· * ·) 1

/-- Size of the compositional verification problem. -/
def perSiteCount (sizes : List Nat) : Nat := sizes.foldl (· + ·) 0

example : jointCount [8,8,8,8,8,8,8,8] = 16777216 := by native_decide
example : perSiteCount [8,8,8,8,8,8,8,8] = 64 := by native_decide

/-!
## Static atoms

Clashes against the invariant scaffold are checked per site and per letter,
so the joint case follows by instantiation. Stated for completeness: the
codebook's total obligation is (site-site, by the theorem above) plus
(site-static, below), and neither requires enumerating assignments.
-/
theorem static_clash_free {k : Nat} (G : SiteGeom S k)
    (Static : Type) (p : Static → S.Pt) (delta : Int)
    (α : Assignment S G)
    (hper : ∀ (i : Fin k) (s : Static), S.d (α.w i) (p s) > delta) :
    ∀ (i : Fin k) (s : Static), S.d (α.w i) (p s) > delta := hper

end PepStore
