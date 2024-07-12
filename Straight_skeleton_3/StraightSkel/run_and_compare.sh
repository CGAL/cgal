#!/bin/bash

SCRIPT_NAME="${0##*/}"

CUSTOMER=""
DATA_PATH=""

BUILD_DIR="build"
OFFSET_DIRECTION="out" # in / out
TIMEOUT_VALUE=10
MAX_ITEM_NUMBER=100

cd "$BUILD_DIR" || exit 1

SUCCESS=()
OFFSET_FAILURE=()
OTHER_FAILURE=()
TIMEOUT=()
BAD_OUTPUT=()

FILES=$(\
  find ${DATA_PATH}/input_* -type f -name '*.ply' | sort \
)

CURRENT_DATE=$(date +%Y-%m-%d)
CURRENT_TIME=$(date +%H:%M)
OUTPUT_DIRECTORY=compare_results_${CURRENT_DATE}_${CURRENT_TIME}
mkdir -p $OUTPUT_DIRECTORY

GLOBAL_LOG=${OUTPUT_DIRECTORY}/log.txt # /dev/stdout
exec > ${GLOBAL_LOG} 2>&1

COUNTER=0

for FILE in $FILES; do
  ((COUNTER++))

  echo "FILE: $FILE" | tee -a "${GLOBAL_LOG}"

  BASE_NAME=$(basename "$FILE" | sed 's/\.[^.]*$//')
  echo "BASE_NAME: ${BASE_NAME}"

  # Extract the number between "input_" and the next underscore
  FULL_ID=$(echo ${BASE_NAME} | sed -E 's/^[^_]*_(.*\..*|.*)/\1/') # input_nnnnnn_n.ply to nnnnnn_n
  ID=$(echo "${BASE_NAME}" | sed 's/^[^_]*_\([^_]*\)_.*$/\1/') # input_nnnnnn_n.ply to nnnnnn, for the offset file
  echo "Extracted FULL_ID: $FULL_ID"
  echo "Extracted ID: $ID"

  mkdir -p $OUTPUT_DIRECTORY/${FULL_ID}

  # ----------------------------------------------------------------
  # Convert to weighted PLY

  WEIGHT_FILE="${DATA_PATH}/offsets_${ID}.txt"
  echo "WEIGHT_FILE: ${WEIGHT_FILE}"

  if [ ! -f $WEIGHT_FILE ]; then
    echo "Missing offset file?!"
    continue
  fi

  WEIGHTED_INPUT="$OUTPUT_DIRECTORY/${FULL_ID}/input.ply"
  echo "WEIGHTED_INPUT: ${WEIGHTED_INPUT}"

  CMD="./convert_to_weighted_PLY ${FILE} ${WEIGHTED_INPUT} ${WEIGHT_FILE}"
  LOG_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/log_conversion.txt

  echo "---- Calling:"
  echo "  $CMD"
  echo "  $LOG_FILE"

  # Run conversion code
  timeout --preserve-status $TIMEOUT_VALUE $CMD > "$LOG_FILE" 2>&1
  RES=$?
  if [ ! "$RES" -eq 0 ]; then
    echo "[ERROR]: failed to create weighted PLY?!"
    OTHER_FAILURE+=("$CMD")
    continue;
  fi

  # ----------------------------------------------------------------
  # If outward offset, add Bbox and invert
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
      echo "[ERROR]: failed to invert and add bbox?!"
      OTHER_FAILURE+=("$CMD")
      continue;
    fi

    WEIGHTED_INPUT=${WEIGHTED_INPUT_WITH_BBOX}
  fi

  # ----------------------------------------------------------------
  # Run skeleton and compare results

  CMD="./StraightSkel 3d load ${WEIGHTED_INPUT} --no-window --save-offsets -1"
  LOG_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/log_offset.txt

  echo "---- Calling:"
  echo "  $CMD"
  echo "  $LOG_FILE"

  # Run offset code
  timeout --preserve-status $TIMEOUT_VALUE $CMD > "$LOG_FILE" 2>&1

  RES=$?
  if [ "$RES" -eq 0 ]; then
    # ----------------------------------------------------------------
    # Got a result, compare it

    OC_OUTPUT=${FILE/input/output}
    echo "OC_OUTPUT: ${OC_OUTPUT}"

    if [ ! -f $OC_OUTPUT ]; then
      echo "[ERROR]: missing OC output file?!"
      OTHER_FAILURE+=("$CMD")
      continue
    fi

    cp ${OC_OUTPUT} ${OUTPUT_DIRECTORY}/${FULL_ID}

    OUTPUT=${OUTPUT_DIRECTORY}/${FULL_ID}/result.ply

    # ----------------------------------------------------------------
    # If outward offset, add Bbox and invert

    if [ "${OFFSET_DIRECTION}" == "out" ]; then
      OUTPUT_WITH_BBOX=${OUTPUT_DIRECTORY}/${FULL_ID}/result_inverted_with_bbox.obj
      cp offset_-1.obj ${OUTPUT_WITH_BBOX}

      CMD="./add_or_remove_bbox ${OUTPUT_WITH_BBOX} remove ${OUTPUT}"
      LOG_FILE=${OUTPUT_DIRECTORY}/${FULL_ID}/log_remove_bbox.txt

      echo "---- Calling:"
      echo "  $CMD"
      echo "  $LOG_FILE"

      # Run code
      timeout --preserve-status $TIMEOUT_VALUE $CMD > "$LOG_FILE" 2>&1
      RES=$?
      if [ ! "$RES" -eq 0 ]; then
        echo "[ERROR]: failed to remove bbox and invert?!"
        OTHER_FAILURE+=("$CMD")
        continue;
      fi
    else
      cp offset_-1.obj ${OUTPUT}
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
      SUCCESS+=("$CMD")
    else
      echo "[ERROR]: output does not match?!"
      BAD_OUTPUT+=("$CMD")
    fi
  elif [ "$RES" -eq 143 ]; then
    echo "====== offset time out! ====== "
    TIMEOUT+=("$CMD")
  else
    echo "====== offset failure! ====== "
    OFFSET_FAILURE+=("$CMD")
  fi
  echo ""

  if [ $COUNTER -ge ${MAX_ITEM_NUMBER} ]; then
    echo "Breaking the loop as the counter has reached ${MAX_ITEM_NUMBER}."
    break
  fi
done

echo ""
echo "${#SUCCESS[@]} SUCCESSES:"
for BASE_NAME in "${SUCCESS[@]}"; do
  echo "$BASE_NAME"
done
echo ""
echo "${#OFFSET_FAILURE[@]} OFFSET FAILURES:"
for BASE_NAME in "${OFFSET_FAILURE[@]}"; do
  echo "$BASE_NAME"
done
echo ""
echo "${#OTHER_FAILURE[@]} OTHER FAILURES:"
for BASE_NAME in "${OTHER_FAILURE[@]}"; do
  echo "$BASE_NAME"
done
echo ""
echo "${#TIMEOUT[@]} TIMEOUTS:"
for BASE_NAME in "${TIMEOUT[@]}"; do
  echo "$BASE_NAME"
done
echo ""
echo "${#BAD_OUTPUT[@]} BAD OUTPUT:"
for BASE_NAME in "${BAD_OUTPUT[@]}"; do
  echo "$BASE_NAME"
done

echo ""
echo "${#SUCCESS[@]} SUCCESSES"
echo "${#OFFSET_FAILURE[@]} OFFSET FAILURES"
echo "${#OTHER_FAILURE[@]} OTHER FAILURES"
echo "${#TIMEOUT[@]} TIMEOUTS"
echo "${#BAD_OUTPUT[@]} BAD OUTPUTS"

SUCCESSES_COUNT=${#SUCCESS[@]}
FAILURES_COUNT=$((${#OFFSET_FAILURE[@]} + ${#OTHER_FAILURE[@]} +${#TIMEOUT[@]} + ${#BAD_OUTPUT[@]}))
TOTAL_COUNT=$((${SUCCESSES_COUNT} + ${FAILURES_COUNT}))
echo "TOTAL: ${TOTAL_COUNT}"
echo "Success %: $((${#SUCCESS[@]} * 100 / ${TOTAL_COUNT})) "
