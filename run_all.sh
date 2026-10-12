#!/usr/bin/env bash
# Runs every check of the archive.
#   ./run_all.sh          the four independent checkers (C, C++, Julia, Python) on cert/certificate.txt,
#                         two from-scratch enumerations of the extreme rays (C, Julia), the cdd/lrs
#                         comparison, the equality classification, the checks of the five extremal
#                         configurations, and the Lean build
#   ./run_all.sh --full   also recomputes the enumeration and the certificate and compares them with the
#                         files in data/ and cert/
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FULL=0
[ "${1:-}" = "--full" ] && FULL=1
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

need() { command -v "$1" >/dev/null 2>&1 || { echo "missing tool: $1" >&2; exit 2; }; }
for t in gcc g++ julia python3 lrs; do need "$t"; done
[ -d "$ROOT/lean" ] && need lake
python3 -c "import cdd.gmp, sympy, numpy" 2>/dev/null || { echo "missing Python modules: pycddlib, sympy, numpy" >&2; exit 2; }

CERT="$ROOT/cert/certificate.txt"
step() { echo; echo "== $*"; }

step "C checker"
gcc -O2 -std=c11 -Wall -o "$WORK/check_c" "$ROOT/checkers/check.c" -lgmp
"$WORK/check_c" "$CERT"

step "C++ checker"
g++ -O2 -std=c++17 -o "$WORK/check_cpp" "$ROOT/checkers/check.cpp" -lgmpxx -lgmp
"$WORK/check_cpp" "$CERT"

step "Julia checker"
julia "$ROOT/checkers/check.jl" "$CERT"

step "Python checker"
python3 -I "$ROOT/checkers/check.py" "$CERT"

step "extreme rays: double description from scratch (C and Julia)"
gcc -O2 -std=c11 -Wall -o "$WORK/dd" "$ROOT/checkers/dd.c" -lgmp
"$WORK/dd" "$CERT"
julia "$ROOT/checkers/dd.jl" "$CERT"

step "extreme rays: cdd against lrs"
(cd "$ROOT/compute" && python3 rays_lrs.py "$ROOT/data/vertices.json")

step "equality points and their classes"
(cd "$ROOT/compute" && python3 equality.py "$ROOT/data/vertices.json" "$WORK/equality_points.json")
cmp <(python3 -c "import json,sys; print(json.dumps(json.load(open(sys.argv[1])), sort_keys=True))" "$WORK/equality_points.json") \
    <(python3 -c "import json,sys; print(json.dumps(json.load(open(sys.argv[1])), sort_keys=True))" "$ROOT/data/equality_points.json")
(cd "$ROOT/compute" && python3 classify_eq.py "$ROOT/data/equality_points.json" "$WORK/optimal_classes.json")
cmp <(python3 -c "import json,sys; print(json.dumps(json.load(open(sys.argv[1])), sort_keys=True))" "$WORK/optimal_classes.json") \
    <(python3 -c "import json,sys; print(json.dumps(json.load(open(sys.argv[1])), sort_keys=True))" "$ROOT/data/optimal_classes.json")
(cd "$ROOT/compute" && python3 lattice_check.py "$ROOT/data/optimal_classes.json")

step "the five configurations of the paper (Section 6)"
python3 -I "$ROOT/compute/configurations.py" "$ROOT/data/equality_points.json"

if [ "$FULL" = 1 ]; then
  step "recompute the enumeration and the certificate"
  (cd "$ROOT/compute" && python3 enumerate.py "$WORK/vertices.json")
  cmp <(python3 -c "import json,sys; print(json.dumps(json.load(open(sys.argv[1])), sort_keys=True))" "$WORK/vertices.json") \
      <(python3 -c "import json,sys; print(json.dumps(json.load(open(sys.argv[1])), sort_keys=True))" "$ROOT/data/vertices.json")
  python3 "$ROOT/cert/make_cert.py" "$WORK/vertices.json" "$ROOT/data/equality_points.json" "$WORK/certificate.txt"
  cmp "$WORK/certificate.txt" "$CERT"
fi

if [ -d "$ROOT/lean" ]; then
  step "Lean: the certificate files are generated from cert/certificate.txt"
  python3 -I "$ROOT/lean/generate/make_lean.py" "$CERT" "$ROOT/data/equality_points.json" "$WORK/lean" >/dev/null
  diff -r "$WORK/lean/D5Lean/Data" "$ROOT/lean/D5Lean/Data"
  diff -r "$WORK/lean/D5Lean/Cert" "$ROOT/lean/D5Lean/Cert"
  step "Lean: build and axioms"
  (cd "$ROOT/lean" && lake exe cache get && lake build && lake build D5Lean.Axioms) | tee "$WORK/lean.log"
  grep -q "'D5.certificate' depends on axioms: \[propext, Classical.choice, Quot.sound\]" "$WORK/lean.log" \
    || { echo "unexpected axioms" >&2; exit 1; }
  if grep -rnE "sorry|native_decide|admit|^axiom" "$ROOT/lean/D5Lean" --include=*.lean; then
    echo "forbidden construct in the Lean sources" >&2; exit 1
  fi
fi

echo
echo "all checks passed"
