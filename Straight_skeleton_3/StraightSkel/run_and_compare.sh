#!/bin/bash

SCRIPT_NAME="${0##*/}"
START_TIME=$(date +%s)

CUSTOMER=""
DATA_PATH=

BUILD_DIR="build-release"
TIMEOUT_VALUE=10
MAX_ITEM_NUMBER=10
NUMBER_OF_THREADS=4
OFFSET_DIRECTION="out" # in / out

cd "$BUILD_DIR" || exit 1

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
OUTPUT_DIRECTORY=compare_results-${CURRENT_DATE}_${CURRENT_TIME}-${MAX_ITEM_NUMBER}i_${TIMEOUT_VALUE}s
mkdir -p $OUTPUT_DIRECTORY

GLOBAL_LOG=${OUTPUT_DIRECTORY}/log.txt # /dev/stdout
exec > ${GLOBAL_LOG} 2>&1

echo "Build DIR = ${BUILD_DIR}"
echo "Timeout = ${TIMEOUT_VALUE}"
echo "Max #inputs = ${MAX_ITEM_NUMBER}"
echo "Number of threads = ${NUMBER_OF_THREADS}"
echo "Inwards / Outwards = ${OFFSET_DIRECTION}"
echo ""

function process_file {
  FILE=$1
  BASE_NAME=$(basename "$FILE" | sed 's/\.[^.]*$//')

  # Extract the number between "input_" and the next underscore
  FULL_ID=$(echo ${BASE_NAME} | sed -E 's/^[^_]*_(.*\..*|.*)/\1/') # input_nnnnnn_n.ply to nnnnnn_n
  ID=$(echo "${BASE_NAME}" | sed 's/^[^_]*_\([^_]*\)_.*$/\1/') # input_nnnnnn_n.ply to nnnnnn, for the offset file

  mkdir -p $OUTPUT_DIRECTORY/${FULL_ID}

  # to write the result, concatenated later into the arrays
  RESULT_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/result.txt

  LOCAL_LOG=${OUTPUT_DIRECTORY}/${FULL_ID}/local_log.txt
  exec > ${LOCAL_LOG} 2>&1

  echo "FILE: $FILE" | tee -a "${LOCAL_LOG}"

  echo "BASE_NAME: ${BASE_NAME}"
  echo "Extracted FULL_ID: $FULL_ID"
  echo "Extracted ID: $ID"

  # ----------------------------------------------------------------
  # Convert to weighted PLY

  WEIGHT_FILE="${DATA_PATH}/offsets_${ID}.txt"
  echo "WEIGHT_FILE: ${WEIGHT_FILE}"

  WEIGHTED_INPUT="$OUTPUT_DIRECTORY/${FULL_ID}/input.ply"
  echo "WEIGHTED_INPUT: ${WEIGHTED_INPUT}"

  CMD="./convert_to_weighted_PLY ${FILE} ${WEIGHTED_INPUT} ${WEIGHT_FILE}"
  LOG_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/log_conversion.txt

  if [ ! -f $WEIGHT_FILE ]; then
    echo "====== [ERROR]: missing offset file?! ======"
    echo "$CMD" > $RESULT_FILE
    echo "PREPROCESS FAILURE" >> $RESULT_FILE
    return
  fi

  echo "---- Calling:"
  echo "  $CMD"
  echo "  $LOG_FILE"

  # Run conversion code
  timeout --preserve-status $TIMEOUT_VALUE $CMD > "$LOG_FILE" 2>&1
  RES=$?
  if [ ! "$RES" -eq 0 ]; then
    echo "====== [ERROR]: failed to create weighted PLY?! ======"
    echo "$CMD" > $RESULT_FILE
    echo "PREPROCESS FAILURE" >> $RESULT_FILE
    return
  fi

  # ----------------------------------------------------------------
  # If outward offset, add bbox and invert
  if [ "${OFFSET_DIRECTION}" == "out" ]; then
    WEIGHTED_INPUT_WITH_BBOX="$OUTPUT_DIRECTORY/${FULL_ID}/input_inverted_with_bbox.ply"
    echo "WEIGHTED_INPUT_WITH_BBOX: ${WEIGHTED_INPUT_WITH_BBOX}"

    CMD="./add_or_remove_bbox ${WEIGHTED_INPUT} add ${WEIGHTED_INPUT_WITH_BBOX}"
    LOG_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/log_add_bbox.txt

    echo "---- Calling:"
    echo "  $CMD"
    echo "  $LOG_FILE"

    # Run code
    timeout --preserve-status $TIMEOUT_VALUE $CMD > "$LOG_FILE" 2>&1
    RES=$?
    if [ ! "$RES" -eq 0 ]; then
      echo "====== [ERROR]: failed to invert and add bbox?! ======"
      echo "$CMD" > $RESULT_FILE
      echo "PREPROCESS FAILURE" >> $RESULT_FILE
      return
    fi

    WEIGHTED_INPUT=${WEIGHTED_INPUT_WITH_BBOX}
  fi

  # ----------------------------------------------------------------
  # Run skeleton and compare results

  # below will write offset_-1.obj and offset_-1_exact.obj into the item's folder
  CMD="./StraightSkel 3d load ${WEIGHTED_INPUT} --no-window --save-offsets -1 --save-path ${OUTPUT_DIRECTORY}/${FULL_ID}"
  LOG_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/log_offset.txt

  echo "---- Calling:"
  echo "  $CMD"
  echo "  $LOG_FILE"

  TIMING_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/timing.txt
  echo "$TIMING_FILE"

  # Run offset code
  RES_CODE_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/res_code.txt
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

  RUNTIME_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/runtime.txt
  echo $runtime > ${RUNTIME_FILE}

  RES=$(cat "$RES_CODE_FILE")
  if [ "$RES" -eq 0 ]; then
    # ----------------------------------------------------------------
    # Got a result, compare it

    OC_OUTPUT=${FILE/input/output}
    echo "OC_OUTPUT: ${OC_OUTPUT}"

    if [ ! -f $OC_OUTPUT ]; then
      echo "[ERROR]: missing OC output file?!"
      echo "$CMD" > $RESULT_FILE
      echo "PREPROCESS FAILURE" >> $RESULT_FILE
      return
    fi

    cp ${OC_OUTPUT} ${OUTPUT_DIRECTORY}/${FULL_ID}

    OUTPUT=${OUTPUT_DIRECTORY}/${FULL_ID}/result.obj

    # ----------------------------------------------------------------
    # If outward offset, add Bbox and invert

    if [ "${OFFSET_DIRECTION}" == "out" ]; then
      OUTPUT_WITH_BBOX=${OUTPUT_DIRECTORY}/${FULL_ID}/offset_-1_exact.obj

      CMD="./add_or_remove_bbox ${OUTPUT_WITH_BBOX} remove ${OUTPUT}"
      LOG_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/log_remove_bbox.txt

      echo "---- Calling:"
      echo "  $CMD"
      echo "  $LOG_FILE"

      # Run code
      timeout --preserve-status $TIMEOUT_VALUE $CMD > "$LOG_FILE" 2>&1
      RES=$?
      if [ ! "$RES" -eq 0 ]; then
        echo "====== [ERROR]: failed to remove bbox and invert?! ======"
        echo "$CMD" > $RESULT_FILE
        echo "POSTPROCESS FAILURE" >> $RESULT_FILE
        return
      fi
    fi

    # ----------------------------------------------------------------
    # Compare

    CMD="./compare_outputs ${OUTPUT} ${OC_OUTPUT}"
    LOG_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/log_comparison.txt

    echo "---- Calling:"
    echo "  $CMD"
    echo "  $LOG_FILE"

    timeout --preserve-status $TIMEOUT_VALUE $CMD > "$LOG_FILE" 2>&1
    RES=$?
    if [ "$RES" -eq 0 ]; then
      echo "====== [SUCCESS] ======"
      echo "$CMD" > $RESULT_FILE
      echo "SUCCESS" >> $RESULT_FILE
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
    echo "====== [ERROR]: offset failure! ======"
    echo "$CMD" > $RESULT_FILE
    echo "OFFSET FAILURE" >> $RESULT_FILE
  fi
  echo ""
}

export DATA_PATH
export TIMEOUT_VALUE
export OFFSET_DIRECTION
export OUTPUT_DIRECTORY
export -f process_file

# main call
echo "$FILES" | parallel --wd $PWD --jobs "$NUMBER_OF_THREADS" --env DATA_PATH --env TIMEOUT_VALUE --env OFFSET_DIRECTION --env OUTPUT_DIRECTORY process_file '{}'

exec >> ${GLOBAL_LOG} 2>&1

# concatenate all local logs into the single global log
find "$OUTPUT_DIRECTORY" -type f -name "local_log.txt" -exec sh -c 'for file do cat "$file"; echo ""; done' sh {} + >> "$GLOBAL_LOG"

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
done < <(find "$OUTPUT_DIRECTORY" -type f -name "result.txt")

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