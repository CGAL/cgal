#!/bin/bash

SCRIPT_NAME="${0##*/}"
START_TIME=$(date +%s)

############## DO NOT COMMIT
CUSTOMER=""
DATA_PATH=/home/mrouxell/Data/Customers/${CUSTOMER}/offset_data_june_2024
# script_test_folder
# offset_data_june_2024
# offset_data_july_2024 # everything has some zero weight there
# offset_data_october_2026
# high_degree_data
############## DO NOT COMMIT

BUILD_DIR="build-release"
TIMEOUT_VALUE=120
MAX_ITEM_NUMBER=20
NUMBER_OF_SEEDS=50
NUMBER_OF_THREADS=2
OFFSET_DIRECTION="out" # in / out

cd "$BUILD_DIR" || exit 1
CONFIG_TEMPLATE="$PWD/StraightSkel.ini"
if [[ ! -r "$CONFIG_TEMPLATE" ]]; then
  echo "Not found or unreadable: $CONFIG_TEMPLATE" >&2
  exit 1
fi

PREPROCESS_FAILURES=()
OFFSET_FAILURES=()
POSTPROCESS_FAILURES=()
OTHER_FAILURES=()
TIMEOUTS=()
BAD_OUTPUTS=()
SUCCESSES=()

FILES=$(find ${DATA_PATH}/input_* -type f -name '*.ply' | sort | head -n "$MAX_ITEM_NUMBER")

CURRENT_DATE=$(date +%Y-%m-%d)
CURRENT_TIME=$(date +%H:%M)
TOP_OUTPUT_DIRECTORY=compare_results-${CURRENT_DATE}_${CURRENT_TIME}-${MAX_ITEM_NUMBER}i_${TIMEOUT_VALUE}s_${NUMBER_OF_SEEDS}seeds
mkdir -p "$TOP_OUTPUT_DIRECTORY"

GLOBAL_LOG=${TOP_OUTPUT_DIRECTORY}/log.txt # /dev/stdout
exec > ${GLOBAL_LOG} 2>&1

echo "Build DIR = ${BUILD_DIR}"
echo "Timeout = ${TIMEOUT_VALUE}"
echo "Max #inputs = ${MAX_ITEM_NUMBER}"
echo "Number of threads = ${NUMBER_OF_THREADS}"
echo "Inwards / Outwards = ${OFFSET_DIRECTION}"
echo ""

function process_file_single_function {
  if [ "${OFFSET_DIRECTION}" != "out" ]; then
    echo "====== [ERROR]: process_file_single_function() is currently only usable if OFFSET_DIRECTION is 'out' ======"
    return
  fi

  FILE=$1
  BASE_NAME=$(basename "$FILE" | sed 's/\.[^.]*$//')

  # Extract the number between "input_" and the next underscore
  FULL_ID=$(echo ${BASE_NAME} | sed -E 's/^[^_]*_(.*\..*|.*)/\1/') # input_nnnnnn_n.ply to nnnnnn_n
  ID=$(echo "${BASE_NAME}" | sed 's/^[^_]*_\([^_]*\)_.*$/\1/') # input_nnnnnn_n.ply to nnnnnn, for the offset file

  LOCAL_OUTPUT_DIRECTORY="${TOP_OUTPUT_DIRECTORY}/${FULL_ID}"
  mkdir -p "$LOCAL_OUTPUT_DIRECTORY"

  LOCAL_LOG="${LOCAL_OUTPUT_DIRECTORY}/local_log.txt"
  exec > "$LOCAL_LOG" 2>&1

  echo "FILE: $FILE" | tee -a "${LOCAL_LOG}"

  echo "BASE_NAME: ${BASE_NAME}"
  echo "Extracted FULL_ID: $FULL_ID"
  echo "Extracted ID: $ID"

  WEIGHT_FILE="${DATA_PATH}/offsets_${ID}.txt"
  echo "WEIGHT_FILE: ${WEIGHT_FILE}"

  for ((seed_index = 1; seed_index <= NUMBER_OF_SEEDS; seed_index++)); do
    SEED=$RANDOM
    OUTPUT_DIRECTORY="${LOCAL_OUTPUT_DIRECTORY}/seed_${seed_index}_${SEED}"
    mkdir -p "$OUTPUT_DIRECTORY"

    SEED_CONFIG="${OUTPUT_DIRECTORY}/StraightSkel.ini"
    cp "$CONFIG_TEMPLATE" "$SEED_CONFIG"
    sed -i -E "s/^([[:space:]]*seed[[:space:]]*=[[:space:]]*).*/\\1${SEED}/" "$SEED_CONFIG"

    RESULT_FILE="${OUTPUT_DIRECTORY}/result.txt"
    LOG_FILE="${OUTPUT_DIRECTORY}/log_offset.txt"
    RES_CODE_FILE="${OUTPUT_DIRECTORY}/res_code.txt"
    TIMING_FILE="${OUTPUT_DIRECTORY}/timing.txt"

    if [ ! -f "$WEIGHT_FILE" ]; then
      echo "====== [INFO]: missing offset file; using test_skeleton_3 ======"
      CMD="./test_skeleton_3 ${FILE} --save-times 1 --save-path ${OUTPUT_DIRECTORY} --config-path ${SEED_CONFIG}"
    else
      CMD="./offset_with_cardinal_weights ${FILE} --weight-path ${WEIGHT_FILE} --save-times 1 --save-path ${OUTPUT_DIRECTORY} --config-path ${SEED_CONFIG}"
    fi

    echo "Seed run ${seed_index}/${NUMBER_OF_SEEDS}: ${SEED}"
    echo "---- Calling:"
    echo "  $CMD"
    echo "  $LOG_FILE"

    time ( timeout --preserve-status $TIMEOUT_VALUE $CMD > "$LOG_FILE" 2>&1; echo $? > "$RES_CODE_FILE" ) 2> "$TIMING_FILE"

  # Extract runtime from the output of timeout
  runtime=$(cat "$TIMING_FILE" | grep "real" | awk '{print $2}')

  # Convert runtime to seconds
  if [[ $runtime =~ m ]]; then
    minutes=$(echo "$runtime" | sed 's/m.*//')
    seconds=$(echo "$runtime" | sed 's/.*m//;s/s//')
    runtime=$(echo "$minutes * 60 + $seconds" | bc)
  else
    runtime=$(echo "$runtime" | sed 's/s//')
  fi

  RUNTIME_FILE=${OUTPUT_DIRECTORY}/runtime.txt
  echo $runtime > ${RUNTIME_FILE}

  # Result of the main function
  RES=$(cat "$RES_CODE_FILE")
  echo "RES: $RES"
  if [ "$RES" -eq 0 ]; then
    # ----------------------------------------------------------------
    # Got a result, compare it

    OC_OUTPUT=${FILE/input/output}
    echo "OC_OUTPUT: ${OC_OUTPUT}"

    if [ ! -f "$OC_OUTPUT" ]; then
      echo "[ERROR]: missing output file?!"
      echo "$CMD" > $RESULT_FILE
      echo "PREPROCESS FAILURE" >> $RESULT_FILE
      return
    fi

    cp "${OC_OUTPUT}" "$OUTPUT_DIRECTORY"

    # ----------------------------------------------------------------
    # Compare

    OUTPUT=${OUTPUT_DIRECTORY}/result_1.obj
    CMD="./compare_outputs ${OUTPUT} ${OC_OUTPUT} ${OUTPUT_DIRECTORY}"
    LOG_FILE=${OUTPUT_DIRECTORY}/log_comparison.txt

    echo "---- Calling:"
    echo "  $CMD"
    echo "  $LOG_FILE"

    timeout --preserve-status $TIMEOUT_VALUE $CMD > "$LOG_FILE" 2>&1
    RES=$?
    if [ "$RES" -eq 0 ]; then
      echo "====== [SUCCESS] ======"
      echo "$CMD" > $RESULT_FILE
      echo "SUCCESS" >> $RESULT_FILE
    elif grep -q -e "Error: 'ours' is BAD" "$LOG_FILE"; then
      echo "$CMD" > $RESULT_FILE
      echo "OFFSET FAILURE" >> $RESULT_FILE
    else
      echo "====== [ERROR]: outputs differ ======"
      echo "$CMD" > $RESULT_FILE
      echo "BAD OUTPUT" >> $RESULT_FILE
    fi
  elif [ "$RES" -eq 143 ]; then
    echo "====== [ERROR]: time out! ======"
    echo "$CMD" > $RESULT_FILE
    echo "TIMEOUT" >> $RESULT_FILE
  elif [ "$RES" -eq 137 ]; then
    echo "====== [ERROR]: process killed! ======"
    echo "$CMD" > $RESULT_FILE
    echo "OTHER FAILURE" >> $RESULT_FILE
  else
    echo "$CMD" > $RESULT_FILE

    # parse the log to detect the error type
    if grep -q -e "Error: Config file not found" \
               -e "Error: failed to read input" \
               -e "Error: invalid input" \
               -e "Error: empty input" \
               -e "Error: open input" \
               -e "Error: broken input" \
               -e "Error: failed to read weights" \
               -e "Error: negative weights" \
               -e "Error: all weights are zero" \
               -e "Error: input has self intersections" \
               -e "Error: input has nested connected components" \
               "$LOG_FILE"; then
      echo "PREPROCESS FAILURE" >> $RESULT_FILE
    elif grep -q -e "Error: failed to assign weights" \
                 -e "Error: failed to preprocess" \
                 -e "Error: failed to add property map" \
                 "$LOG_FILE"; then
        echo "PREPROCESS FAILURE" >> $RESULT_FILE
    elif grep -q -e "Error: failed to initialize polyhedron" \
                 "$LOG_FILE"; then
        echo "OFFSET FAILURE" >> $RESULT_FILE
    elif grep -q -e "Error: run() returned 'false'" \
                 -e "Error: broken output" \
                 -e "Error: open output" \
                 -e "Error: non-triangulated output" \
                 -e "Error: degenerate faces in output" \
                 -e "Error: output self-intersects" \
                 -e "Error: failed to write result" \
                 "$LOG_FILE"; then
        echo "POSTPROCESS FAILURE" >> $RESULT_FILE
    else
        echo "OFFSET FAILURE" >> $RESULT_FILE
    fi
  fi
  echo ""
  done
}

export DATA_PATH
export TIMEOUT_VALUE
export OFFSET_DIRECTION
export TOP_OUTPUT_DIRECTORY
export CONFIG_TEMPLATE
export NUMBER_OF_SEEDS
export -f process_file_single_function

# main call
echo "$FILES" | parallel --wd $PWD --jobs "$NUMBER_OF_THREADS" --env DATA_PATH --env TIMEOUT_VALUE --env OFFSET_DIRECTION --env TOP_OUTPUT_DIRECTORY --env CONFIG_TEMPLATE --env NUMBER_OF_SEEDS process_file_single_function '{}'

exec >> ${GLOBAL_LOG} 2>&1

# concatenate all local logs into the single global log
find "$TOP_OUTPUT_DIRECTORY" -type f -name "local_log.txt" -exec sh -c 'for file do cat "$file"; echo ""; done' sh {} + >> "$GLOBAL_LOG"

# extract individual info into the arrays
while IFS= read -r FILE; do
  {
    read -r FIRST_LINE
    read -r SECOND_LINE
  } < "$FILE"

  case "$SECOND_LINE" in
    "PREPROCESS FAILURE")
      PREPROCESS_FAILURES+=("$FIRST_LINE")
      ;;
    "OFFSET FAILURE")
      OFFSET_FAILURES+=("$FIRST_LINE")
      ;;
    "POSTPROCESS FAILURE")
      POSTPROCESS_FAILURES+=("$FIRST_LINE")
      ;;
    "TIMEOUT")
      TIMEOUTS+=("$FIRST_LINE")
      ;;
    "BAD OUTPUT")
      BAD_OUTPUTS+=("$FIRST_LINE")
      ;;
    "OTHER FAILURE")
      OTHER_FAILURES+=("$FIRST_LINE")
      ;;
    "SUCCESS")
      SUCCESSES+=("$FIRST_LINE")
      ;;
    *)
      echo "Unknown category: $SECOND_LINE in file $FILE"
      ;;
  esac
done < <(find "$TOP_OUTPUT_DIRECTORY" -type f -name "result.txt")

echo ""
echo "${#SUCCESSES[@]} SUCCESSES:"
for BASE_NAME in "${SUCCESSES[@]}"; do
  echo "$BASE_NAME"
done
echo ""
echo "${#OFFSET_FAILURES[@]} OFFSET FAILURES:"
for BASE_NAME in "${OFFSET_FAILURES[@]}"; do
  echo "$BASE_NAME"
done
echo ""
echo "${#PREPROCESS_FAILURES[@]} PREPROCESS FAILURES:"
for BASE_NAME in "${PREPROCESS_FAILURES[@]}"; do
  echo "$BASE_NAME"
done
echo ""
echo "${#POSTPROCESS_FAILURES[@]} POSTPROCESS FAILURES:"
for BASE_NAME in "${POSTPROCESS_FAILURES[@]}"; do
  echo "$BASE_NAME"
done
echo ""
echo "${#OTHER_FAILURES[@]} OTHER FAILURES:"
for BASE_NAME in "${OTHER_FAILURES[@]}"; do
  echo "$BASE_NAME"
done
echo ""
echo "${#TIMEOUTS[@]} TIMEOUTS:"
for BASE_NAME in "${TIMEOUTS[@]}"; do
  echo "$BASE_NAME"
done
echo ""
echo "${#BAD_OUTPUTS[@]} OUTPUTS DIFFER:"
for BASE_NAME in "${BAD_OUTPUTS[@]}"; do
  echo "$BASE_NAME"
done

echo ""
echo "${#SUCCESSES[@]} SUCCESSES"
echo "${#OFFSET_FAILURES[@]} OFFSET FAILURES"
echo "${#PREPROCESS_FAILURES[@]} PREPROCESS FAILURES"
echo "${#POSTPROCESS_FAILURES[@]} POSTPROCESS FAILURES"
echo "${#OTHER_FAILURES[@]} KILLED"
echo "${#TIMEOUTS[@]} TIMEOUTS"
echo "${#BAD_OUTPUTS[@]} OUTPUTS DIFFER"

SUCCESSES_COUNT=${#SUCCESSES[@]}
FAILURES_COUNT=$((${#OFFSET_FAILURES[@]} + ${#PREPROCESS_FAILURES[@]} + ${#POSTPROCESS_FAILURES[@]} + ${#OTHER_FAILURES[@]} +${#TIMEOUTS[@]}))
TOTAL_COUNT=$((${SUCCESSES_COUNT} + ${FAILURES_COUNT} + ${#BAD_OUTPUTS[@]}))
echo "TOTAL: ${TOTAL_COUNT}"
echo "Success %: $((${#SUCCESSES[@]} * 100 / ${TOTAL_COUNT}))"

ELAPSED=$(($(date +%s) - START_TIME))
printf "elapsed: %s\n\n" "$(date -d@$ELAPSED -u +%H\ hours\ %M\ min\ %S\ sec)"