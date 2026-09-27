#!/usr/bin/env bash
# Bake the Alcedo Studio LUT set (print-chain LMTs, short names).
#   scripts/aces_lmt/bake_alcedo_set.sh <out_dir> [extra bake_print_drt.py options]
# Examples:
#   scripts/aces_lmt/bake_alcedo_set.sh experiment_results/luts_alcedo/plain
#   scripts/aces_lmt/bake_alcedo_set.sh experiment_results/luts_alcedo/N  --neutralize
#   scripts/aces_lmt/bake_alcedo_set.sh experiment_results/luts_alcedo/NH --neutralize --hue-preserve 0.25
set -euo pipefail
OUT="$1"; shift
PY="${PYTHON:-.venv/Scripts/python.exe}"
[ -x "$PY" ] || PY=.venv/bin/python
BAKE=(scripts/aces_lmt/bake_print_drt.py --out "$OUT" "$@")
CINE="kodak_vision3_50d kodak_vision3_200t kodak_vision3_250d kodak_vision3_500t kodak_verita_200d"
KSTILL="kodak_portra_160 kodak_portra_400 kodak_portra_800 kodak_portra_800_push1 kodak_portra_800_push2 kodak_ektar_100 kodak_gold_200 kodak_ultramax_400"
FSTILL="fujifilm_c200 fujifilm_pro_400h fujifilm_xtra_400"
SLIDE="fujifilm_velvia_100 fujifilm_provia_100f kodak_ektachrome_100 kodak_kodachrome_64"
"$PY" "${BAKE[@]}" $CINE --print kodak_2383
"$PY" "${BAKE[@]}" $CINE --print kodak_2393
"$PY" "${BAKE[@]}" $KSTILL
"$PY" "${BAKE[@]}" $KSTILL --print fujifilm_crystal_archive_typeii
"$PY" "${BAKE[@]}" $FSTILL
"$PY" "${BAKE[@]}" $FSTILL --print kodak_portra_endura
"$PY" "${BAKE[@]}" $SLIDE
