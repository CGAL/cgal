#!/bin/bash

# Temp directory for individual result JSONs
TMP_RESULT_DIR=$(mktemp -d)

# Job control
JOBS=0
MAX_JOBS=$NUM_THREADS

# Function to process a single file
process_file() {
    INPUT_PATH="$1"
    INPUT_ID=$(basename "$INPUT_PATH" | cut -d. -f1)
    COMPONENT_NAME="$2"
    PROJECT_DIR="$3"
    TIMEOUT="$4"
    OUTPUT_DIR="$5"
    TMP_RESULT_FILE="$6"
    GRID_SIZE="$7"
    ERASE_ALL_DUPLICATE="$8"
    {
        echo "    \"$INPUT_ID\": {"
        echo "      \"path\": \"$INPUT_PATH\","

        PERF_OUTPUT=$(bash "$PROJECT_DIR/Performance/run_performance.sh" "$INPUT_PATH" "$TIMEOUT" "$GRID_SIZE" "$ERASE_ALL_DUPLICATE" 2>> "$OUTPUT_DIR/Logs/$COMPONENT_NAME/Performance/$INPUT_ID.log")
        echo "      \"Performance\": $PERF_OUTPUT,"

        QUALITY_OUTPUT=$(bash "$PROJECT_DIR/Quality/run_quality.sh" "$INPUT_PATH" "$TIMEOUT" "$GRID_SIZE" "$ERASE_ALL_DUPLICATE" 2>> "$OUTPUT_DIR/Logs/$COMPONENT_NAME/Quality/$INPUT_ID.log")
        echo "      \"Quality\": $QUALITY_OUTPUT,"

        ROBUST_OUTPUT=$(bash "$PROJECT_DIR/Robustness/run_robustness.sh" "$INPUT_PATH" "$TIMEOUT" "$GRID_SIZE" "$ERASE_ALL_DUPLICATE" 2>> "$OUTPUT_DIR/Logs/$COMPONENT_NAME/Robustness/$INPUT_ID.log")
        echo "      \"Robustness\": $ROBUST_OUTPUT"

        echo "    }"
    } > "$TMP_RESULT_FILE"
}
export -f process_file

# Usage function
usage() {
    echo "Usage: $0 <project_dir> <input_data_dir> <output_results_dir> <timeout> <num_threads> [component_params...]"
    exit 1
}

# Check parameters
if [ "$#" -lt 5 ]; then
    usage
fi

# Arguments
PROJECT_DIR=$1
INPUT_DIR=$2
OUTPUT_DIR=$3
TIMEOUT=$4
NUM_THREADS=$5
GRID_SIZE=$6
ERASE_ALL_DUPLICATE=$7

# Get component name from the project directory name
COMPONENT_NAME=$(basename "$PROJECT_DIR")
DATE_TAG=$(date +"%Y-%m-%d")
TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
RESULT_JSON="$OUTPUT_DIR/${COMPONENT_NAME}_results_${DATE_TAG}.json"

# Compile
# Do not forget to define CGAL_DIR
cmake "$PROJECT_DIR" "-DCMAKE_BUILD_TYPE=Release" "-DCMAKE_CXX_FLAGS=-O3"
make -j $NUM_THREADS

# Prepare log directories
mkdir -p "$OUTPUT_DIR/Logs/$COMPONENT_NAME/Performance"
mkdir -p "$OUTPUT_DIR/Logs/$COMPONENT_NAME/Quality"
mkdir -p "$OUTPUT_DIR/Logs/$COMPONENT_NAME/Robustness"

# Initialize JSON
echo "{" > "$RESULT_JSON"
echo "  \"$COMPONENT_NAME\": {" >> "$RESULT_JSON"
echo "    \"Thingi10K\": {" >> "$RESULT_JSON"

#process_file "$INPUT_DIR/100036.stl" "$COMPONENT_NAME" "$PROJECT_DIR" "$TIMEOUT" "$OUTPUT_DIR" "$TMP_RESULT_FILE" "$GRID_SIZE" "$ERASE_ALL_DUPLICATE"
# Loop input files and spawn parallel jobs
for INPUT_FILE in "$INPUT_DIR"/*; do
    INPUT_ID=$(basename "$INPUT_FILE" | cut -d. -f1)
    TMP_RESULT_FILE="$TMP_RESULT_DIR/$INPUT_ID.json"

    process_file "$INPUT_FILE" "$COMPONENT_NAME" "$PROJECT_DIR" "$TIMEOUT" "$OUTPUT_DIR" "$TMP_RESULT_FILE" "$GRID_SIZE" "$ERASE_ALL_DUPLICATE"

    ((JOBS+=1))
    if [ "$JOBS" -ge "$NUM_THREADS" ]; then
        wait
        JOBS=0
    fi
done

wait

# Merge all partial JSONs
echo "{" > "$RESULT_JSON"
echo "  \"$COMPONENT_NAME\": {" >> "$RESULT_JSON"
echo "    \"Thingi10K\": {" >> "$RESULT_JSON"

FIRST_ENTRY=true
for FILE in "$TMP_RESULT_DIR"/*.json; do
    if [ "$FIRST_ENTRY" = true ]; then
        FIRST_ENTRY=false
    else
        echo "," >> "$RESULT_JSON"
    fi
    cat "$FILE" >> "$RESULT_JSON"
done

echo "" >> "$RESULT_JSON"
echo "    }," >> "$RESULT_JSON"
echo "    \"finished_at\": \"$TIMESTAMP\"" >> "$RESULT_JSON"
echo "  }" >> "$RESULT_JSON"
echo "}" >> "$RESULT_JSON"