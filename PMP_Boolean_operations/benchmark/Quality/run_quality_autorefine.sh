#!/bin/bash

set -e

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_file> <timeout>"
    exit 1
fi

INPUT_FILE=$1
TIMEOUT=$2

TMP_LOG=$(mktemp)
timeout "$TIMEOUT"s quality_autorefine "$INPUT_FILE" > "$TMP_LOG"
cat "$TMP_LOG"
rm -f "$TMP_LOG"
