#!/bin/sh

echo "Input r1 r2 r3 r4 r5 (on same line):"
read r1 r2 r3 r4 r5
echo "Input size_bound (for r<r1):"
read size

filename_base="spheres-$r1-$r2-$r3-$r4-$r5-$size"

echo "$r1 $r2 $r3 $r4 $r5 $size ${filename_base}.mesh" \
  > ${filename_base}.cin

echo "./${filename_base}.cin created."
