/-
Kernel checker for the certificate of the two-translate packing problem in R^5.
Core Lean only. Matrices are lists of rows of integers; `ent A i j` is the entry (i, j),
`0` outside the list. Every function returns a `Bool`; `D5Lean.Sound` proves what a
`true` result means.
-/
namespace D5

abbrev Mat := List (List Int)

/-- Entry `(i, j)` of a list matrix. -/
def ent (A : Mat) (i j : Nat) : Int := (A.getD i []).getD j 0

/-- Coordinate `i` of a list vector. -/
def vg (v : List Int) (i : Nat) : Int := v.getD i 0

/-- `sumN n f = f 0 + ⋯ + f (n - 1)`. -/
def sumN : Nat → (Nat → Int) → Int
  | 0, _ => 0
  | n + 1, f => sumN n f + f n

/-- `prodN n f = f 0 * ⋯ * f (n - 1)`. -/
def prodN : Nat → (Nat → Int) → Int
  | 0, _ => 1
  | n + 1, f => prodN n f * f n

/-- `allN n p` holds when `p i` holds for every `i < n`. -/
def allN : Nat → (Nat → Bool) → Bool
  | 0, _ => true
  | n + 1, p => allN n p && p n

/-- Dot product of two lists (truncated to the shorter one). -/
def dot : List Int → List Int → Int
  | a :: as, b :: bs => a * b + dot as bs
  | _, _ => 0

/-- The quadratic form `J[k]`, computed by dot products. -/
def qf (A : Mat) (k : List Int) : Int := dot k (A.map fun r => dot r k)

/-- `A` has `n` rows, each of length `n`. -/
def shape (A : Mat) (n : Nat) : Bool := A.length == n && A.all fun r => r.length == n

def symm (A : Mat) : Bool := allN 6 fun i => allN 6 fun j => ent A i j == ent A j i

/-- Entry `(i, j)` of the product `A * B` of `n × n` matrices. -/
def mulE (n : Nat) (A B : Mat) (i j : Nat) : Int := sumN n fun k => ent A i k * ent B k j

/-- Entry `(i, j)` of `Tᵀ * A * T`. -/
def congE (T A : Mat) (i j : Nat) : Int :=
  sumN 6 fun a => sumN 6 fun b => ent T a i * ent A a b * ent T b j

/-! ### Membership in `R` by a box enumeration -/

def rangeAll (p : Int → Bool) : Int → Nat → Bool
  | _, 0 => true
  | v, k + 1 => p v && rangeAll p (v + 1) k

/-- `boxAll [(lo₁, n₁), …] p`: `p` holds at every integer point with `lo_i ≤ x_i < lo_i + n_i`. -/
def boxAll : List (Int × Nat) → (List Int → Bool) → Bool
  | [], p => p []
  | (lo, n) :: rest, p => rangeAll (fun v => boxAll rest (fun xs => p (v :: xs))) lo n

def negV (k : List Int) : List Int := k.map fun x => -x

def isZero (k : List Int) : Bool := k.all fun x => x == 0

/-- At a lattice point of the box: `k = 0`, or `J[k] > 4`, or `J[k] = 4` and `±k` is listed. -/
def ptOK (J : Mat) (mins : List (List Int)) (k : List Int) : Bool :=
  isZero k || decide (4 < qf J k) || (qf J k == 4 && (mins.contains k || mins.contains (negV k)))

/-- A listed minimal vector: of length 6, last coordinate 0 or 1, nonzero, and `J[k] = 4`. -/
def minOK (J : Mat) (k : List Int) : Bool :=
  k.length == 6 && (vg k 5 == 0 || vg k 5 == 1) && !isZero k && qf J k == 4

/-- The bound certificate for coordinate `i`: `C · J' = ∑ e_r m_r m_rᵀ` with `e_r ≥ 0`,
where `J' = cd · J + cn · E₅₅`, the first term is `(e, unit vector i)` and
`e (b + 1)² > C (4 cd + cn)`. -/
def boundOK (J : Mat) (cd cn C : Int) (i : Nat) (b : Nat)
    (terms : List (Int × List Int)) : Bool :=
  match terms with
  | [] => false
  | (e, m) :: _ =>
    m.length == 6 && (allN 6 fun a => vg m a == if a == i then 1 else 0) &&
    decide (C * (4 * cd + cn) < e * ((b : Int) + 1) * ((b : Int) + 1)) &&
    (terms.all fun t => decide (0 ≤ t.1) && t.2.length == 6) &&
    (allN 6 fun a => allN 6 fun c =>
      C * (cd * ent J a c + (if a == 5 && c == 5 then cn else 0)) ==
        (terms.map fun t => t.1 * vg t.2 a * vg t.2 c).foldr (· + ·) 0)

/-- The box `[-b_i, b_i]` in the first five coordinates and `{0, 1}` in the last. -/
def boxOf (bs : List Nat) : List (Int × Nat) :=
  bs.map (fun (b : Nat) => ((-(b : Int)), 2 * b + 1)) ++ [((0 : Int), (2 : Nat))]

/-- `J ∈ R`, with `mins` the minimal vectors: the five bound certificates and the box. -/
def memOK (J : Mat) (mins : List (List Int)) (cd cn C : Int) (bs : List Nat)
    (certs : List (List (Int × List Int))) : Bool :=
  decide (0 < cd) && decide (0 ≤ cn) && decide (0 < C) && bs.length == 5 && certs.length == 5 &&
  (allN 5 fun i => boundOK J cd cn C i (bs.getD i 0) (certs.getD i [])) &&
  mins.all (minOK J) && boxAll (boxOf bs) (ptOK J mins)

/-! ### Determinants by triangular elimination -/

/-- `E` is lower triangular, `E * A` is upper triangular, `∏ diag(E A) = d ∏ diag E` and the same
for the leading `5 × 5` blocks with `dq`, and `∏ diag E ≠ 0`. Then `det A = d`, `det A₅ = dq`. -/
def detOK (A E : Mat) (d dq : Int) : Bool :=
  (allN 6 fun i => allN 6 fun j => decide (j ≤ i) || ent E i j == 0) &&
  (allN 6 fun i => allN 6 fun j => decide (i ≤ j) || mulE 6 E A i j == 0) &&
  prodN 6 (fun i => ent E i i) != 0 &&
  prodN 6 (fun i => mulE 6 E A i i) == d * prodN 6 (fun i => ent E i i) &&
  prodN 5 (fun i => mulE 6 E A i i) == dq * prodN 5 (fun i => ent E i i)

/-! ### Polynomials: coefficient lists, lowest degree first -/

def powI (w : Int) : Nat → Int
  | 0 => 1
  | n + 1 => powI w n * w

/-- `hom p y w n = ∑ p_i y^i w^(n-i)`, the value `w^n p(y / w)`. -/
def hom : List Int → Int → Int → Nat → Int
  | [], _, _, _ => 0
  | a :: r, y, w, n => a * powI w n + y * hom r y w (n - 1)

def addP : List Int → List Int → List Int
  | [], q => q
  | p, [] => p
  | a :: p, b :: q => (a + b) :: addP p q

/-- Coefficients of `(u + h x) c(x)`. -/
def mulLin (u h : Int) (c : List Int) : List Int := addP (c.map (u * ·)) (0 :: c.map (h * ·))

/-- Coefficients in `x` of `hom p (u + h x) w n`. -/
def comp : List Int → Int → Int → Int → Nat → List Int
  | [], _, _, _, _ => []
  | a :: r, u, h, w, n => addP [a * powI w n] (mulLin u h (comp r u h w (n - 1)))

def negPart : List Int → Int
  | [] => 0
  | a :: r => min a 0 + negPart r

/-- Positive on `(0, 1]`: after the leading zeros, `c_m + ∑_{i > m} min(c_i, 0) > 0`. -/
def posBound : List Int → Bool
  | [] => false
  | a :: r => if a == 0 then posBound r else decide (0 < a + negPart r)

/-- Positive on `(0, ∞)`: all coefficients `≥ 0` and the first nonzero one `> 0`. -/
def posInf : List Int → Bool
  | [] => false
  | a :: r => if a == 0 then posInf r else decide (0 < a) && r.all fun c => decide (0 ≤ c)

/-- `g = det Q - 512` from the coefficients of `det Q`. -/
def gOf : List Int → List Int
  | [] => [-512]
  | a :: r => (a - 512) :: r

/-! ### Edges -/

/-- The value of the polynomial `p` (degree bound `n`) or of `g` on a piece. -/
def pieceOK (pJ g : List Int) (wd : Int) (modeJ dirL : Bool) (u v : Int) : Bool :=
  let P := if modeJ then pJ else g
  let n := if modeJ then 6 else 5
  if dirL then posBound (comp P u (v - u) wd n) else posBound (comp P v (u - v) wd n)

/-- At the breakpoint `t = u / wd`: `det J ≤ 0 → det Q ≥ 512`, and `det J = 0 ∧ det Q = 512`
only at a listed point. -/
def pointOK (pJ g : List Int) (wd : Int) (eqs : List Int) (u : Int) : Bool :=
  let a := hom pJ u wd 6
  let b := hom g u wd 5
  (decide (0 < a) || decide (0 ≤ b)) && (!(a == 0 && b == 0) || eqs.contains u)

def sorted : List Int → Bool
  | a :: b :: r => decide (a < b) && sorted (b :: r)
  | _ => true

/-- Consecutive breakpoints and the pieces between them. -/
def piecesOK (pJ g : List Int) (wd : Int) : List Int → List (Bool × Bool) → Bool
  | u :: v :: r, (m, d) :: ps => pieceOK pJ g wd m d u v && piecesOK pJ g wd (v :: r) ps
  | [_], [] => true
  | _, _ => false

/-- The sample points `t_j = sn_j / sd` with their elimination certificates for `sd J + sn_j Y`. -/
def samplesOK (J Y : Mat) (pJ pQ : List Int) (sd : Int) (sn : List Int) (Es : List Mat) : Bool :=
  decide (0 < sd) && sn.length == 7 && Es.length == 7 && sorted sn &&
  (allN 7 fun j =>
    let s := vg sn j
    detOK ((List.range 6).map fun a => (List.range 6).map fun b => sd * ent J a b + s * ent Y a b)
      (Es.getD j []) (hom pJ s sd 6) (hom pQ s sd 5))

/-- The sign certificate along an edge: breakpoints `us / wd` from `0`, pieces between them. -/
def signOK (pJ pQ : List Int) (wd : Int) (us : List Int) (ps : List (Bool × Bool))
    (eqs : List Int) : Bool :=
  let g := gOf pQ
  decide (0 < wd) && pJ.length ≤ 7 && pQ.length ≤ 6 && us.head? == some 0 && sorted us &&
  us.all (pointOK pJ g wd eqs) && piecesOK pJ g wd us ps

/-- A bounded edge `J + tY`, `0 ≤ t ≤ tn / td`: the end point is `Tᵀ J_to T` with `T ∈ Γ`
(inverse `Tinv`), and `kmax` is a minimal vector of the end point with `Y[kmax] < 0`. -/
def boundedOK (J Y Jto T Tinv : Mat) (tn td : Int) (kmax : List Int) : Bool :=
  shape Y 6 && symm Y && shape T 6 && shape Tinv 6 && decide (0 < tn) && decide (0 < td) &&
  (allN 6 fun i => allN 6 fun j => td * ent J i j + tn * ent Y i j == td * congE T Jto i j) &&
  (allN 6 fun i => allN 6 fun j => mulE 6 T Tinv i j == if i == j then 1 else 0) &&
  (allN 5 fun j => ent T 5 j == 0) && (ent T 5 5 == 1 || ent T 5 5 == -1) &&
  kmax.length == 6 && (vg kmax 5 == 0 || vg kmax 5 == 1) && !isZero kmax &&
  td * qf J kmax + tn * qf Y kmax == 4 * td && decide (qf Y kmax < 0)

/-- An unbounded edge: `Y = [[λ w wᵀ, μ w], [μ wᵀ, s]]` with `λ n² + 2 μ n + s ≥ 0` on `ℤ`
(checked at `n₀ = ⌊-μ/λ⌋` and `n₀ + 1`), so `Y[k] ≥ 0` for every `k`. -/
def unboundedOK (Y : Mat) (lam mu s n0 : Int) (w : List Int) : Bool :=
  shape Y 6 && symm Y && w.length == 5 &&
  (allN 5 fun a => allN 5 fun b => ent Y a b == lam * vg w a * vg w b) &&
  (allN 5 fun a => ent Y a 5 == mu * vg w a) && ent Y 5 5 == s && decide (0 ≤ lam) &&
  (if lam == 0 then mu == 0 && decide (0 ≤ s) else
    decide (lam * n0 ≤ -mu) && decide (-mu < lam * (n0 + 1)) &&
    decide (0 ≤ lam * n0 * n0 + 2 * mu * n0 + s) &&
    decide (0 ≤ lam * (n0 + 1) * (n0 + 1) + 2 * mu * (n0 + 1) + s))

/-- The sign certificate for an unbounded edge ends with the piece `[u_last, ∞)`. -/
def tailOK (pJ pQ : List Int) (wd : Int) (us : List Int) (modeJ : Bool) : Bool :=
  match us.getLast? with
  | none => false
  | some u => if modeJ then posInf (comp pJ u 1 wd 6) else posInf (comp (gOf pQ) u 1 wd 5)

/-! ### The certificate records -/

structure Vtx where
  J : Mat
  mins : List (List Int)
  cd : Int
  cn : Int
  C : Int
  bs : List Nat
  certs : List (List (Int × List Int))
  E : Mat
  dJ : Int
  dQ : Int
  eq : Bool
deriving Inhabited

structure Edge where
  Y : Mat
  bounded : Bool
  tn : Int
  td : Int
  target : Nat
  T : Mat
  Tinv : Mat
  kmax : List Int
  lam : Int
  mu : Int
  s : Int
  n0 : Int
  w : List Int
  pJ : List Int
  pQ : List Int
  sd : Int
  sn : List Int
  Es : List Mat
  wd : Int
  us : List Int
  ps : List (Bool × Bool)
  tailJ : Bool
  eqs : List Int
deriving Inhabited

def Vtx.ok (v : Vtx) : Bool :=
  shape v.J 6 && symm v.J && memOK v.J v.mins v.cd v.cn v.C v.bs v.certs &&
  detOK v.J v.E v.dJ v.dQ &&
  (decide (0 < v.dJ) || decide (512 ≤ v.dQ)) && ((v.dJ == 0 && v.dQ == 512) == v.eq)

def Edge.ok (J Jto : Mat) (e : Edge) : Bool :=
  samplesOK J e.Y e.pJ e.pQ e.sd e.sn e.Es && signOK e.pJ e.pQ e.wd e.us e.ps e.eqs &&
  if e.bounded then
    boundedOK J e.Y Jto e.T e.Tinv e.tn e.td e.kmax && e.us.getLast? == some (e.tn * e.wd / e.td) &&
      e.tn * e.wd % e.td == 0
  else
    unboundedOK e.Y e.lam e.mu e.s e.n0 e.w && tailOK e.pJ e.pQ e.wd e.us e.tailJ

end D5
