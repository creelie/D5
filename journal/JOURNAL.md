# Target journal: Discrete & Computational Geometry

| Item | Details |
|---|---|
| Journal | Discrete & Computational Geometry (DCG) |
| Publisher | Springer Nature (Springer US) |
| ISSN | 0179-5376 (print), 1432-0444 (electronic) |
| Homepage | https://link.springer.com/journal/454 |
| Submission | Editorial Manager, https://www.editorialmanager.com/dcge; the form asks which Co-Editor-in-Chief should handle the paper (Janos Pach or Kenneth Clarkson) |
| Review | the guidelines ask for the authors' names, affiliations and the corresponding author's e-mail on the title page and say nothing about anonymization; submit `dist/d5-two-translates.pdf` |
| Source files | LaTeX; Springer recommends, and does not require, its Springer Nature LaTeX template. The paper uses `amsart`; editable source files must be uploaded at every submission and revision (`dist/d5-two-translates-tex.zip`) |

## Scope fit

DCG publishes research in mathematics and computer science in which geometry
plays a central role, and sphere packing is a standing topic of the journal.
The two papers the article builds on appeared there: A. Andreanov and
Y. Kallus, *Locally optimal 2-periodic sphere packings*, Discrete Comput.
Geom. 63 (2020), 182-208, which reported the five-dimensional bound from a
partly floating-point enumeration, and H. Cohn and I. Rajagopal, *Variations
on five-dimensional sphere packings*, Discrete Comput. Geom. (2026), which
cites that bound as a theorem. The article proves the bound and determines all
equality cases, with an exact certificate checked by four independent programs
and by the Lean 4 kernel.

MSC 2020: primary 52C17; secondary 11H31, 11H55.

## Length and format

- 17 pages in `amsart`; the journal sets no page limit.
- The source is a single file `main.tex` with the references in the source
  (`thebibliography`), so no `.bib` or `.bbl` file is needed;
  `dist/d5-two-translates-tex.zip` and the arXiv tarball contain this file
  only.
- Title page: title, author, affiliation, e-mail address of the corresponding
  author and ORCID, as the guidelines ask. Six keywords (the guidelines ask
  for 4 to 6) and the MSC codes are in the first-page footnote.
- The abstract has 76 words. Springer's general instructions ask for 150 to
  250 words; lengthen it at submission if the editorial office asks.
- The paper has the sections "Statements and Declarations" (funding,
  competing interests, data and code availability), which the guidelines
  require, and "Use of AI tools".
- Supplementary material: upload `dist/MasterD5TwoTranslatesZip.zip`, or cite
  the Zenodo record of the repository, as the data availability statement does.

## Policy on AI assistance

- Springer Nature: "Large Language Models (LLMs), such as ChatGPT, do not
  currently satisfy our authorship criteria. ... Use of an LLM should be
  properly documented in the Methods section (and if a Methods section is not
  available, in a suitable alternative part) of the manuscript." AI-assisted
  copy editing need not be declared. (DCG submission guidelines,
  https://link.springer.com/journal/454/submission-guidelines, section
  "Title Page"; read on 2026-10-12.)
- The paper has a section "Use of AI tools" that names Claude Code
  (Anthropic) and says how it was used; the author takes responsibility for
  the content. Add the version of the tool if the editors ask for it.
- Every step of the computation of Proposition 4.3 is checked by four
  independent programs and recomputed by two further programs and two
  standard packages; the Lean 4 kernel checks the certificate of every vertex
  and every edge (Appendix B states exactly what Lean covers). The reduction
  (Sections 2 and 3), the equality case (Section 5) and the identification of
  the configurations (Section 6) are proved by hand in the paper.
