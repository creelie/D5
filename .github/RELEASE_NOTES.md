**Densest packings of two translates of a lattice in five dimensions**, by Deep Bhattacharjee (corresponding author) and Ushashi Bhattacharya.

**Theorem.** Let Λ ⊂ R^5 be a lattice and b ∈ R^5 \ Λ such that distinct points of P = Λ ∪ (Λ + b) are at distance at least 2. Then covol(Λ) ≥ 16√2, so the packing of unit balls centred at P has density at most π²/(15√2), the density of the D_5 lattice packing. Equality holds for exactly five pairs (Λ, b) up to congruence, giving three point sets: √2 D_5, the uniform packing with Szöllősi's kissing configuration Q_5, and Leech's 2-periodic layered packing with kissing configuration L_5. Andreanov and Kallus (Discrete Comput. Geom. 63 (2020)) reported this from a partly floating-point enumeration.

| Step | Content |
|---|---|
| Section 2 | the polyhedron R of 6 × 6 Gram matrices, its vertices, edges and symmetries |
| Proposition 3.2 | every packing moves, without increasing det Q, to a vertex of R or an interior point of an edge (concavity of log det Q and of the Schur complement) |
| Proposition 4.3 | R has 34 vertices up to equivalence; 80,695 extreme rays in 480 orbits, 25 unbounded (Table 4.2) |
| Proposition 4.5, Corollary 4.6 | det Q ≥ 512 at every candidate point, with equality at exactly five points E_0, …, E_4 |
| Proposition 5.3 | equality case, by hand: a deep hole of a lattice is surrounded by lattice points spanning the space affinely |
| Section 6 | the five configurations, their bases, and the identification with Q_5, L_5 and the index-two sublattices of D_5 |

The computation is in exact arithmetic. Four independent programs (C, C++, Julia, Python) check the certificate; two double-description programs (C, Julia), cddlib and lrs recompute all extreme rays. Lean 4 with Mathlib proves `D5.certificate`: the minimal vectors and determinants of every vertex and every statement of the certificate on the 480 representative edges, each checked by `decide +kernel`, without `native_decide`, using only the axioms propext, Classical.choice and Quot.sound. Not formalized: the rank conditions, the completeness of the ray lists and the orbit decompositions, which the programs check; the reduction and the equality case are proved by hand in the paper.

Version 2.0.0 adds the second author, three schematic figures, the declarations of the journal and the preprint folder.

Files:
- `AndreanovKallusTwoPeriodicPackingD5_final2.pdf`: the paper
- `AndreanovKallusTwoPeriodicPackingD5_final2-tex.zip`: LaTeX source (`main.tex`, figures as PNG)
- `AndreanovKallusTwoPeriodicPackingD5_final2-arxiv.tar.gz`: the same source, ready for arXiv
- `MasterAndreanovKallusTwoPeriodicPackingD5_masterzipfinal2.zip`: every file of the repository with the files above and `preprintAndreanovKallusTwoPeriodicPackingD5/AndreanovKallusTwoPeriodicPackingD5.pdf`

Run `run_all.sh` to repeat the checks and `build_paper.sh` to rebuild the files above.

Zenodo, all versions: [10.5281/zenodo.23323366](https://doi.org/10.5281/zenodo.23323366).
