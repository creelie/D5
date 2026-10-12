# Densest packings of two translates of a lattice in five dimensions

Paper, certificate, verification programs and Lean 4 formalization by
Deep Bhattacharjee.

**Theorem.** Let Λ ⊂ R^5 be a lattice and b ∈ R^5 \ Λ such that distinct
points of P = Λ ∪ (Λ + b) are at distance at least 2. Then
covol(Λ) ≥ 16√2. Equality holds for exactly five pairs (Λ, b) up to
congruence: three splittings of √2 D_5 into two cosets of a sublattice of
index two, and two non-lattice packings, the uniform packing with kissing
configuration Q_5 and Leech's 2-periodic layered packing with kissing
configuration L_5. So the packing of unit balls centred at P has density at
most π²/(15√2), the density of the D_5 lattice packing, and every Voronoi
cell of P has volume at least 8√2.

The bound was reported by Andreanov and Kallus (Discrete Comput. Geom. 63
(2020)) from an enumeration that partly used floating-point arithmetic. The
paper proves it and determines all equality cases.

## Method

A pair (Λ, b) is a point of a locally finite polyhedron R in the
21-dimensional space of symmetric 6 × 6 matrices. A concavity argument
(Section 3) reduces the bound to the vertices of R and the interior points of
its edges. Up to the integral symmetries, R has 34 vertices; their cones have
80,695 extreme rays in 480 orbits, 25 of them unbounded edges (Section 4). The
equality case is proved by hand through the deep holes of the lattice
(Section 5), and the five configurations are identified in Section 6.

## Files

| Path | Content |
|---|---|
| `paper/main.tex`, `paper/main.pdf` | the paper (amsart, 17 pages) |
| `dist/` | the PDF, the LaTeX source zip, the arXiv tarball and `MasterD5TwoTranslatesZip.zip` with every file of the repository |
| `journal/JOURNAL.md` | target journal (Discrete & Computational Geometry), format and AI policy |
| `compute/` | the enumeration (`enumerate.py`) and the analysis of the equality points (`equality.py`, `classify_eq.py`, `lattice_check.py`, `configurations.py`) |
| `data/` | the 34 vertices with their rays and edges (`vertices.json`), the five equality points and their classes |
| `cert/certificate.txt` | the certificate (4 MB): for each vertex its entries, minimal vectors, symmetry generators, all extreme rays, and one edge per orbit with its end point and a matrix in Γ, or its unbounded direction X_{u,a} |
| `checkers/` | four independent checkers of the certificate in C, C++ (GMP rationals), Julia and Python (exact rationals), and two double-description enumerations of the extreme rays from scratch in C and Julia |
| `lean/` | Lean 4 + Mathlib project; `D5.certificate` proves, in the kernel, every vertex and edge statement of the certificate |
| `run_all.sh` | runs every check |
| `build_paper.sh` | builds the paper and the files in `dist/` |

## Running the checks

```
./run_all.sh          # the four checkers, the two ray enumerations, cddlib against lrs,
                      # the equality points and the five configurations, and the Lean build
./run_all.sh --full   # also recomputes the enumeration and the certificate and compares
                      # them byte for byte with data/ and cert/
```

Needs gcc and g++ with GMP, Julia, Python 3 with pycddlib (GMP build), sympy
and numpy, `lrs` (lrslib), and elan for Lean (`lean/lean-toolchain`:
leanprover/lean4:v4.34.1).

## What Lean checks

`D5.certificate` (in `lean/D5Lean/Main.lean`) states, for the 34 vertices and
the 480 representative edges: each J_i lies in R, its minimal vectors are the
listed ones, det J_i and det Q(J_i) have the listed values; along every
representative edge, det J > 0 or det Q ≥ 512, with det J = 0 and
det Q = 512 only at the listed points; a bounded edge ends at T^T J_j T with
T ∈ Γ and leaves R after it; an unbounded direction Y has Y[k] ≥ 0 for all
k in M*. Each vertex and each edge is one theorem `checker = true`, proved by
`decide +kernel`, so the kernel evaluates the Boolean checker of
`lean/D5Lean/Check.lean`; the soundness lemmas in `lean/D5Lean/Sound/` derive
the statements above from it. `native_decide` is not used, and `#print axioms D5.certificate` gives propext, Classical.choice and
Quot.sound. The rank conditions that make the J_i vertices and the rays
extreme, the completeness of the ray lists and the orbit decompositions are
checked by the programs, not in Lean; the reduction and the equality case are
proved by hand in the paper.

## License

MIT (see `LICENSE`).
