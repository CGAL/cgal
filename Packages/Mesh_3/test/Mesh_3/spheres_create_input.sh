#!/bin/bash

echo "Input r1 r2 r3 r4 r5 (on same line):"
read r1 r2 r3 r4 r5

echo "$r1 $r2 $r3 $r4 $r5 spheres-$r1-$r2-$r3-$r4-$r5.mesh" \
  > spheres-$r1-$r2-$r3-$r4-$r5.cin

echo "./spheres-$r1-$r2-$r3-$r4-$r5.cin created."
