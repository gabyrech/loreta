#!/usr/bin/env bash
# Usage: maybe_emit.sh <names.txt> <out.fastq.gz>
# Reads the embedded heredoc as the command to run IF names.txt is non-empty.
set -euo pipefail
NAMES="$1"
OUTFQ="$2"
shift 2
CMD=$(cat)
if [ -s "$NAMES" ]; then
bash -c "$CMD"
else
# emit empty gz so downstream steps don't fail
printf "" | gzip -c > "$OUTFQ"
fi