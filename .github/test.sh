#!/bin/bash
FACTOR=$1
set -ex
cd Polyhedron/demo
LIST_OF_PLUGINS=$(for f in $(find ./Polyhedron/Plugins -iname "CMakeLists.txt" |egrep -v Three_examples); do\
  egrep "polyhedron_demo_plugin" $f |cut -d\( -f2 | cut -d" " -f1|sort -u|egrep -v "partition_plugin"\
 |egrep -v "basic_item_plugin"|egrep -v "vtk_plugin" |egrep -v "register_point_sets_plugin"| egrep -v "classification_plugin" | egrep -v "io_image_plugin" ||true; done)
PLUGINS_ARRAY=(${LIST_OF_PLUGINS});
NB_OF_PLUGINS=${#PLUGINS_ARRAY[@]}
DEL=$(($NB_OF_PLUGINS / 4))
mkdir build
cd build
/usr/local/bin/cmake -DCGAL_DIR=$2 ../Polyhedron
make -j2 ${PLUGINS_ARRAY[@]:$(($FACTOR * $DEL)):$((($FACTOR + 1) * $DEL))}
