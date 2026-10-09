#!/bin/bash
set -e
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_file> <timeout>"
    exit 1
fi
INPUT_FILE=$1
TIMEOUT=$2
TMP_LOG=$(mktemp)
/usr/bin/time -f "TIME:%e\nMEM:%M" timeout "$TIMEOUT"s performance_boolop "$INPUT_FILE" 2> "$TMP_LOG"
SECONDS=$(grep "TIME" "$TMP_LOG" | cut -d':' -f2)
MEMORY=$(grep "MEM" "$TMP_LOG" | cut -d':' -f2)
rm -f "$TMP_LOG"
echo "{\"seconds\": \"$SECONDS\", \"memory_peaks\": \"$MEMORY\"}"
