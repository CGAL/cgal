#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_file> <timeout>"
    exit 1
fi

INPUT_FILE=$1
TIMEOUT=$2

timeout --foreground "$TIMEOUT"s robustness_autorefine "$INPUT_FILE"
EXIT_CODE=$?

declare -A TAGS
TAGS[0]="VALID_OUTPUT"
TAGS[1]="INPUT_IS_INVALID"
TAGS[3]="SELF_INTERSECTING_OUTPUT"
TAGS[139]="SIGSEGV"
TAGS[11]="SIGSEGV"
TAGS[6]="SIGABRT"
TAGS[8]="SIGFPE"
TAGS[132]="SIGILL"
TAGS[124]="TIMEOUT"

TAG_NAME=${TAGS[$EXIT_CODE]:-UNKNOWN}
TAG_DESC=$([[ "$EXIT_CODE" -eq 0 ]] && echo "OK" || echo "Error")
echo "{\"TAG_NAME\": \"$TAG_NAME\", \"TAG\": \"$TAG_DESC\"}"
