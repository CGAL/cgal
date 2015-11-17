#!/bin/bash

rm error.txt

for exe in test_intersection_of_polyhedra test_corefinement; do
  echo  "-------------------------------"
  echo  "||  Runnning " $exe "  ||"
  echo  "-------------------------------"

  k=`wc -l coref_testsuite_run_list | awk '{print $1}'`

  for i in `seq 1 $k`; do
    files=`head -n $i testsuite_run_list | tail -n 1`
    f1=`echo $files | awk '{print $1}'`
    f2=`echo $files | awk '{print $2}'`
    echo "========== " $f1 $f2 " =========="
    if ./$exe $f1 $f2; then
      echo ./$exe $f1 $f2 OK >> error.txt
    else
      echo ./$exe $f1 $f2 Error >> error.txt
    fi
    echo "========== " $f2 $f1 " =========="
    if ./$exe $f2 $f1; then
      echo ./$exe $f2 $f1 OK >> error.txt
    else
      echo ./$exe $f2 $f1 Error >> error.txt
    fi
  done
  echo ""
done

