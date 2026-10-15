# Target journal: Discrete & Computational Geometry

| Item | Details |
|---|---|
| Journal | Discrete & Computational Geometry (DCG) |
| Publisher | Springer Nature (Springer US) |
| ISSN | 0179-5376 (print), 1432-0444 (electronic) |
| Homepage | https://link.springer.com/journal/454 |
| Submission | Editorial Manager, https://www.editorialmanager.com/dcge |
| Manuscript | `dist/AndreanovKallusTwoPeriodicPackingD5_final2.pdf` (18 pages, `amsart`) |
| Source files | `dist/AndreanovKallusTwoPeriodicPackingD5_final2-tex.zip` (figures as PNG); or `build_submission.sh`, which writes `build/submission/` with the PDF and a source zip holding `main.tex` and the figures as Fig1 to Fig3 (PDF and EPS) |

## Scope fit

DCG publishes research in which geometry plays a central role, and sphere
packing is a standing topic of the journal. The two papers the article builds
on appeared there: A. Andreanov and Y. Kallus, *Locally optimal 2-periodic
sphere packings*, Discrete Comput. Geom. 63 (2020), 182-208, which reported
the five-dimensional bound from a partly floating-point enumeration, and
H. Cohn and I. Rajagopal, *Variations on five-dimensional sphere packings*,
Discrete Comput. Geom. 76 (2026), 1660-1678, which cites that bound as a
theorem. The article proves the bound and determines all equality cases, with
an exact certificate checked by four independent programs and, for its
vertex and edge statements, by the Lean 4 kernel.

MSC 2020: primary 52C17; secondary 11H31, 11H55.

## Length and format

- 18 pages; the journal sets no page limit.
- One source file `main.tex` with the references in `thebibliography`.
- Two authors: Deep Bhattacharjee (corresponding author) and Ushashi
  Bhattacharya, with affiliations, e-mail addresses and ORCID iDs. Six
  keywords (the guidelines ask for 4 to 6) and the MSC codes are in the
  first-page footnote.
- The abstract has two sentences (48 words). Springer's guidelines ask for
  150 to 250 words; the authors chose a short abstract, so lengthen it only if
  the editor asks.
- Three figures: two translates of a lattice, the descent to the vertices and
  edges of the polyhedron, and the deep hole lemma in the plane.
- The section "Statements and Declarations" has competing interests (none),
  funding (none) and data availability (who computed what, the Zenodo DOI of
  the Lean formalization and the link to the PDF in `preprintAndreanovKallusTwoPeriodicPackingD5/`), as
  the guidelines require.
- Supplementary material: the Zenodo record of the repository,
  doi:10.5281/zenodo.23387908 (release v2.0.0; all versions doi:10.5281/zenodo.23323366),
  cited in the paper.

## Policy on AI assistance

- Springer Nature: large language models do not satisfy its authorship
  criteria, and their use is documented in the Methods section or a suitable
  alternative part of the manuscript (DCG submission guidelines,
  https://link.springer.com/journal/454/submission-guidelines).
- The acknowledgements say that the authors used Claude (Anthropic),
  through Claude Code, for the proofs, the Lean formalization, the
  computations and the preparation of the article, and take full
  responsibility for it; the data availability statement repeats it.

## What the Lean formalization covers

Appendix B of the paper states it exactly: the minimal vectors and
determinants of the 34 vertices and every statement of the certificate on the
480 representative edges. The rank conditions, the completeness of the lists
of extreme rays and the orbit decompositions are checked by the programs in
exact arithmetic and are not formalized.
