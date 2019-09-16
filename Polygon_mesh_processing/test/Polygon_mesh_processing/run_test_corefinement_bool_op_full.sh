#!/bin/bash

GREEN="\\033[1;32m"
NORMAL="\\033[0;39m"
RED="\\033[1;31m"

k=`wc -l test_corefinement_bool_op_full.cmd | awk '{print $1}'`

for i in `seq 1 $k`; do
  files=`head -n $i test_corefinement_bool_op_full.cmd | tail -n 1`
  f1=`echo $files | awk '{print $1}'`
  f2=`echo $files | awk '{print $2}'`
  ru=`echo $files | awk '{print $4}'`
  ri=`echo $files | awk '{print $5}'`
  rm=`echo $files | awk '{print $6}'`
  rmr=`echo $files | awk '{print $7}'`
  echo -n "==== " $f1 $f2 " "

  if (./test_corefinement_bool_op $f1 $f2 ALL $ru $ri $rm $rmr|| false ) > /dev/null 2>&1; then
    echo -e " ==> $GREEN SUCCEED"
    echo -e -n "$NORMAL"
  else
    echo -e " ==> $RED FAILED"
    echo -e -n "$NORMAL"
  fi
done
