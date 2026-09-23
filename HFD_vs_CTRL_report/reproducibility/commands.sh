#!/usr/bin/env bash
# Replay this ClawBio proteomics-de run
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="$(dirname "$SCRIPT_DIR")"
: "${CLAWBIO_ROOT:=C:\Users\Ibrah\.gemini\config\plugins\science}"
if [ ! -d "$CLAWBIO_ROOT" ]; then
  echo "Invalid CLAWBIO_ROOT: $CLAWBIO_ROOT" >&2
  exit 1
fi

python "$CLAWBIO_ROOT/skills\proteomics-de\proteomics_de.py" \
  --input \
  "diann_matrix.tsv" \
  --input-type \
  diann \
  --metadata \
  "metadata.csv" \
  --contrast \
  HFD,CTRL \
  --s0 \
  0.1 \
  --fdr \
  0.05 \
  --ttest-df \
  4 \
  --imputation-shift \
  1.8 \
  --imputation-scale \
  0.3 \
  --output \
  "$OUTPUT_DIR"
