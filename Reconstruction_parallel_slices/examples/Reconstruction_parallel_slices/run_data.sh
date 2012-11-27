#!/bin/sh

mkdir -p res

for i in `cat $1.txt`; do
  prefix=`echo $i | sed 's/\.cgal//'`
  echo ./conforming-Delaunay_3 $1/$i
  if ./conforming-Delaunay_3 $1/$i; then
    mv graph.off res/$prefix.off
  else
    touch res/error-$prefix
  fi
done
