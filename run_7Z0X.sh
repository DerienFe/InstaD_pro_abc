#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for chain_id in A B C; do
  python "$ROOT_DIR/extract_fasta.py" \
    --structure "$ROOT_DIR/data/7Z0X.cif" \
    --chain "$chain_id" \
    --label 7Z0X \
    --output-dir "$ROOT_DIR/data"
done
