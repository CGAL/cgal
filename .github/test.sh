#!/bin/bash
/usr/local/bin/cmake --version
FACTOR=$1
set -ex
cd Lab/demo
/usr/local/bin/cmake -S Lab -B build -DCGAL_DIR=$2
LIST_OF_PLUGINS=$(/usr/local/bin/cmake --build build -t help | egrep 'plugin$' |& cut -d\  -f2)
PLUGINS_ARRAY=(${LIST_OF_PLUGINS});
NB_OF_PLUGINS=${#PLUGINS_ARRAY[@]}
DEL=$(($NB_OF_PLUGINS / 4))
cd build
NUM_PROCS=$(nproc)
make -j${NUM_PROCS} ${PLUGINS_ARRAY[@]:$(($FACTOR * $DEL)):$((($FACTOR + 1) * $DEL))}
