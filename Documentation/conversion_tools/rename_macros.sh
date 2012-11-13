#!/bin/bash

nb_lines=`wc -l rename_macros.input | awk '{print $1}'`
sort rename_macros.input -r > /tmp/sorted_rename_macros.input

for i in `seq 1 $nb_lines`; do
  l=`head -n $i /tmp/sorted_rename_macros.input | tail -n 1`
  echo $l | awk '{print "s/\\\\" $1 "/\\\\" $2 "/g"}'
done > /tmp/rename_macros.sed

find ../.. -name '*.txt' -o -name '*.h' | xargs sed -i -f /tmp/rename_macros.sed


for i in `seq 1 $nb_lines`; do
  l=`head -n $i /tmp/sorted_rename_macros.input | tail -n 1`
  echo $l | awk '{print "s/\\\"" $1 "=/\\\"" $2 "=/g"}'
  echo $l | awk '{print "s/\\\"" $1 "{/\\\"" $2 "{/g"}'
done > /tmp/rename_macros.sed

sed -i -f /tmp/rename_macros.sed ../Doxyfile

echo "Do not forget to manually update ../pkglist_filter.py and ../doxyassist.xml"
