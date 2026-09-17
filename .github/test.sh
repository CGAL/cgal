#!/bin/bash
/usr/local/bin/cmake --version
FACTOR=$1
set -ex
cd Lab/demo
/usr/local/bin/cmake -S Lab -B build -DCGAL_DIR="$2"
mapfile -t PLUGINS_ARRAY <<< "$(/usr/local/bin/cmake --build build -t help | grep -E 'plugin$' |& cut -d\  -f2)"
NB_OF_PLUGINS=${#PLUGINS_ARRAY[@]}
DEL=$(( NB_OF_PLUGINS / 4 + 1))
cd build
NUM_PROCS=$(nproc)
make "-j${NUM_PROCS}" "${PLUGINS_ARRAY[@]:$((FACTOR * DEL)):$DEL}"
