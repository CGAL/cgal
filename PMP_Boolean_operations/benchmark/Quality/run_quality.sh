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

TMP_LOG=$(mktemp)
timeout "$TIMEOUT"s quality_snap_polygon_soup "$INPUT_FILE" "$GRID_SIZE" "$ERASE_ALL_DUPLICATE"  > "$TMP_LOG"

cat $TMP_LOG
rm -f "$TMP_LOG"