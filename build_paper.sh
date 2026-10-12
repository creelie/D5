#!/usr/bin/env bash
# build_paper.sh -- build the paper and the release files into dist/.
#
#   dist/d5-two-translates.pdf             the compiled paper
#   dist/d5-two-translates-tex.zip         LaTeX source (main.tex)
#   dist/d5-two-translates-arxiv.tar.gz    the same source, ready for arXiv
#   dist/MasterD5TwoTranslatesZip.zip      every file of the repository with the three files above,
#                                          in the folder MasterD5TwoTranslatesZip/
#
# Needs pdflatex (amsart), git, zip and tar.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAPER="$ROOT/paper"
DIST="$ROOT/dist"
NAME=d5-two-translates
MASTER=MasterD5TwoTranslatesZip
for t in pdflatex git zip tar; do
  command -v "$t" >/dev/null || { echo "missing tool: $t" >&2; exit 2; }
done
mkdir -p "$DIST"

latex() { pdflatex -interaction=nonstopmode -halt-on-error "$@" > /dev/null; }
build() { # build <jobname> <input>
  rm -f "$1.aux" "$1.out"
  for run in 1 2 3 4; do
    latex -jobname="$1" "$2"
    grep -q "Rerun to get" "$1.log" || break
  done
  if grep -q "Overfull\|Underfull\|undefined\|Rerun to get" "$1.log"; then
    grep "Overfull\|Underfull\|undefined\|Rerun to get" "$1.log"
    echo "build_paper.sh: fix the warnings above" >&2
    exit 1
  fi
  rm -f "$1.aux" "$1.log" "$1.out"
}

echo "== paper"
cd "$PAPER"
build main main.tex
cp main.pdf "$DIST/$NAME.pdf"

STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

echo "== tex.zip"
mkdir -p "$STAGE/$NAME"
cp main.tex "$STAGE/$NAME/"
rm -f "$DIST/$NAME-tex.zip"
( cd "$STAGE" && zip -qrX "$DIST/$NAME-tex.zip" "$NAME" )

echo "== arXiv tarball"
mkdir -p "$STAGE/arxiv"
cp main.tex "$STAGE/arxiv/"
( cd "$STAGE/arxiv" && latex main.tex && latex main.tex && rm -f main.aux main.log main.out main.pdf )
tar -czf "$DIST/$NAME-arxiv.tar.gz" -C "$STAGE/arxiv" main.tex

echo "== $MASTER.zip"
cd "$ROOT"
mkdir -p "$STAGE/master/$MASTER"
git ls-files --cached --others --exclude-standard -z \
  | grep -zv "^dist/" \
  | while IFS= read -r -d '' f; do
      [ -e "$f" ] || continue
      mkdir -p "$STAGE/master/$MASTER/$(dirname "$f")"
      cp -p "$f" "$STAGE/master/$MASTER/$f"
    done
mkdir -p "$STAGE/master/$MASTER/dist"
cp "$DIST/$NAME.pdf" "$DIST/$NAME-tex.zip" "$DIST/$NAME-arxiv.tar.gz" "$STAGE/master/$MASTER/dist/"
rm -f "$DIST/$MASTER.zip"
( cd "$STAGE/master" && zip -qrX "$DIST/$MASTER.zip" "$MASTER" )

ls -l "$DIST"
