#!/bin/bash
# Requirement:
# * data: a link to the directory containing /big/*.off big-genus/*.off and big-genus-with-holes/*.off
# * quadrangulation_computation_benchmarks: a link to the exe of benchmark/quadrangulation_computation_benchmarks.cpp
# * path_homotopy: a link to the exe of examples/Surface_mesh_topology/path_homotopy.cpp 
# * path_homotopy_with_schema: a link to the exe of benchmark/path_homotopy_with_schema.cpp
#
# Output:
# 3 files for the 3 benchmarks:
# * res-quadrangulation-computation.txt for bench1
# * res-path-homotopy.txt for bench2
# * res-polygonal-schema.txt for bench3
# These files can be transform as data for gnuplot, using the script 

# You can run all the benchmarks or only one.
BENCH=2 # One work among ALL 1 2 3 

# Bench 1
if [ ${BENCH} = ALL -o ${BENCH} = 1 ]
then
    all_files=`ls data/big/*.off data/big-genus/*.off data/big-genus-with-holes/*.off`
    ./quadrangulation_computation_benchmarks ${all_files} > res-quadrangulation-computation.txt
fi

# Bench 2
if [ ${BENCH} = ALL -o ${BENCH} = 2 ]
then
    all_files="data/big/happy.off data/big-genus/obj9.off data/big-genus-with-holes/obj10-1.000000-2.off"
    rm -f res-path-homotopy.txt
    for file in ${all_files}
    do
        ./path_homotopy -L 100 10000 -D 100 10000 -N 100 $file -time >> res-path-homotopy.txt
    done
fi

# Bench 3
if [ ${BENCH} = ALL -o ${BENCH} = 3 ]
then
    rm -f res-polygonal-schema.txt
    for k in 10 20 30 100 200 300 1000 2000 3000 10000 20000 30000 100000 200000 300000 1000000 2000000 3000000 10000000 20000000 30000000
    do
        ./path_homotopy_with_schema "a b -a -b c d -c -d e f -e -f g h -g -h i j -i -j" -nbedges ${k} -nbdefo 100 -nbtests 10 -time >> res-polygonal-schema.txt
    done
fi
