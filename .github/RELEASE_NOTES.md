**Densest packings of two translates of a lattice in five dimensions**, by Deep Bhattacharjee.

**Theorem.** Let Λ ⊂ R^5 be a lattice and b ∈ R^5 \ Λ such that distinct points of P = Λ ∪ (Λ + b) are at distance at least 2. Then covol(Λ) ≥ 16√2, so the packing of unit balls centred at P has density at most π²/(15√2), the density of the D_5 lattice packing. Equality holds for exactly five pairs (Λ, b) up to congruence: three splittings of √2 D_5 and two non-lattice packings (the uniform Q_5 packing and Leech's 2-periodic L_5 packing).

| Step | Content |
|---|---|
| Section 2 | the polyhedron R of 6 × 6 Gram matrices, its vertices, edges and symmetries |
| Proposition 3.2 | every packing moves, without increasing det Q, to a vertex of R or an interior point of an edge (concavity of log det Q and of the Schur complement) |
| Proposition 4.3 | R has 34 vertices up to equivalence; 80,695 extreme rays in 480 orbits, 25 unbounded (Table 4.2) |
| Proposition 4.5, Corollary 4.6 | det Q ≥ 512 at every candidate point, with equality at exactly five points E_0, …, E_4 |
| Proposition 5.3 | equality case, by hand: a deep hole of a lattice is surrounded by lattice points spanning the space affinely |
| Section 6 | the five configurations, their bases, and the identification with Q_5, L_5 and the index-two sublattices of D_5 |

The computation is in exact arithmetic. Four independent programs (C, C++, Julia, Python) check the certificate; two double-description programs (C, Julia), cddlib and lrs recompute all extreme rays. Lean 4 with Mathlib proves `D5.certificate`: every vertex and edge statement of the certificate, each checked by `decide +kernel`, without `native_decide`, using only the axioms propext, Classical.choice and Quot.sound. The rank conditions, the completeness of the ray lists and the orbit decompositions are checked by the programs, and the reduction and the equality case are proved by hand in the paper.

Files:
- `dist/d5-two-translates.pdf`: the paper
- `dist/d5-two-translates-tex.zip`: LaTeX source (`main.tex`)
- `dist/d5-two-translates-arxiv.tar.gz`: the same source, ready for arXiv
- `dist/MasterD5TwoTranslatesZip.zip`: every file of the repository with the three files above

Run `run_all.sh` to repeat the checks and `build_paper.sh` to rebuild the files above.
