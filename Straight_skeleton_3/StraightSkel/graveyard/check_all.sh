#!/bin/sh

cd build

if [ -z "$POLYHEDRONS" ]; then
  POLYHEDRONS="
    $(ls ../res/polyhedrons/events/*.obj)
    $(ls ../res/polyhedrons/*.obj)
    $(ls ../res/polyhedrons/wolfgang/*.obj)"
fi

for POLYHEDRON in $POLYHEDRONS; do
  OUTPUT=$(./StraightSkel 3d load $POLYHEDRON --no-window --check-graph \
    | grep 'GraphChecker.cpp' | sed 's#\[DEBUG\] /.*: ##')
  RESULT=$(echo "$OUTPUT" | grep 'result=')
  NUM_LAYERS_VANISH=$(echo "$OUTPUT" | grep 'num_layers_vanish=')
  NUM_LAYERS_CONTACT=$(echo "$OUTPUT" | grep 'num_layers_contact=')
  echo "$POLYHEDRON & $RESULT & $NUM_LAYERS_VANISH & $NUM_LAYERS_CONTACT \\\\"
done
