#!/bin/bash

if ! git diff --exit-code > /dev/null || ! git diff --staged --exit-code  > /dev/null ; then
    echo 'Your working directory contains local modifications!' >&2
    exit 1
fi

if [ ! -f Installation/include/CGAL/version.h ]; then
 echo "This script should be run at the root of a CGAL branch" >&2
 exit 1
fi

# replace tabs by two spaces in AABB-tree and Minkowski_sum_2 packages
find AABB_tree -name '*.h' -o -name '*.cpp' -o -name '*.txt' | xargs sed -i -E 's/\t/  /g'
find Minkowski_sum_2 -name '*.h' -o -name '*.cpp' -o -name '*.txt' | xargs sed -i -E 's/\t/  /g'

#replace tabs by 8 spaces for all other packages
find . -name '*.h' -o -name '*.cpp' -o -name '*.hpp' -o -name '*.tcc' -o -name '*.txt' | xargs sed -i -E 's/\t/        /g'

#remove trailing whitespace
find . -name '*.h' -o -name '*.cpp' -o -name '*.hpp' -o -name '*.tcc' -o -name '*.txt' | xargs sed -i -E 's/\s+$//'
