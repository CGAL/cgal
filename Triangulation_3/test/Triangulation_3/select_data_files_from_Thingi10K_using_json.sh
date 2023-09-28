#!/bin/sh

NUM_MAX=${1:-10000}


IDS=$(jq -c "[input_filename, .num_vertices < $NUM_MAX and .solid == 1]" json/* | awk -F'[/.,]' '/,true/ { print $2 }')
pushd raw_meshes >/dev/null
for i in $IDS; do
    echo $i.*
done
