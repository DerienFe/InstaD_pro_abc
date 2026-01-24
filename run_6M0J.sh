#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python "$ROOT_DIR/extract_fasta.py" \
  --structure "$ROOT_DIR/data/6M0J.cif" \
  --chain A \
  --label 6M0J \
  --output-dir "$ROOT_DIR/data"

python "$ROOT_DIR/extract_fasta.py" \
  --structure "$ROOT_DIR/data/6M0J.cif" \
  --chain E \
  --label 6M0J \
  --output-dir "$ROOT_DIR/data"
