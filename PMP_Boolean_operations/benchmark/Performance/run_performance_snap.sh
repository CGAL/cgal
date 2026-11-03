#!/bin/bash

set -e

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <input_file> <timeout> [component_params...]"
    exit 1
fi

INPUT_FILE=$1
TIMEOUT=$2
GRID_SIZE=$3
ERASE_ALL_DUPLICATE=$4

# Use /usr/bin/time for memory usage (maximum resident set size in KB)
TMP_LOG_BOOLOP=$(mktemp)
TMP_LOG_AUTOREFINE=$(mktemp)
TMP_LOG_SNAP=$(mktemp)


# Run the benchmarked command
/usr/bin/time -f "TIME:%e\nMEM:%M" timeout "$TIMEOUT"s performance_boolop "$INPUT_FILE" 2> "$TMP_LOG_BOOLOP"
/usr/bin/time -f "TIME:%e\nMEM:%M" timeout "$TIMEOUT"s performance_autorefine "$INPUT_FILE" 2> "$TMP_LOG_AUTOREFINE"
/usr/bin/time -f "TIME:%e\nMEM:%M" timeout "$TIMEOUT"s performance_snap_polygon_soup "$INPUT_FILE" "$GRID_SIZE" "$ERASE_ALL_DUPLICATE" 2> "$TMP_LOG_SNAP"

# Parse time and memory
SECONDS=$(grep "TIME" "$TMP_LOG" | cut -d':' -f2)
MEMORY=$(grep "MEM" "$TMP_LOG" | cut -d':' -f2)

rm -f "$TMP_LOG"

# Output JSON
echo "{\"seconds\": \"$SECONDS\", \"memory_peaks\": \"$MEMORY\"}"
