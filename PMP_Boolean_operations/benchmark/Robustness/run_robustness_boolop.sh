#!/bin/bash
set -e
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_file> <timeout>"
    exit 1
fi
INPUT_FILE=$1
TIMEOUT=$2
timeout --foreground "$TIMEOUT"s robustness_boolop "$INPUT_FILE"
echo '{"TAG_NAME": "VALID_OUTPUT", "TAG": "OK"}'
