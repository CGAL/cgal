#!/bin/bash

BASE_DIR=data-coref/coplanar_triangles/all_cases

rm error.txt

for k in  $BASE_DIR \
          $BASE_DIR/bugreport \
          $BASE_DIR/deg \
          $BASE_DIR/deg/collinear_edges \
          $BASE_DIR/deg/vertex_on_edge; do
  nf=`ls $k/*.off | wc -l`
  nf=`expr $nf '/' 2`
  echo $k
  for i in `seq 1 $nf`; do
      if ./intersection_of_polyhedra_coplanar_tests $k/tr$i-1.off $k/tr$i-2.off ; then 
        echo $k/tr$i OK >> error.txt
      else
        echo $k/tr$i Error >> error.txt
      fi
      if ./test_intersection_of_polyhedra $k/tr$i-1.off $k/tr$i-2.off ; then 
        echo $k/tr$i AF OK >> error.txt
      else
        echo $k/tr$i AF Error >> error.txt
      fi
  done
done