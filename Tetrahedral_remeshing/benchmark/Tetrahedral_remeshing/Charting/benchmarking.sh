#!/bin/bash

# ###################################################################################################
# # INTRODUCTION
# ###################################################################################################

# This script is meant to benchmark Tetrahedral_remeshing over a directory of data files (e.g. Thingi10k/raw_meshes)
# with unique names. It checks a number of metrics in 3 categories: Robustness (did it finish without
# errors, and did it produce a valid result?), Performance (how long did it take to run?), and Quality
# (how good is the result mesh?).
#
# This script gives access to the remeshing criteria for convenience, but all the other parameters
# are fixed and chosen in two configuration files.
#
# If `benchmark_tetrahedral_remeshing.cpp` is in parallel mode, you should use a single thread.
# You also probably don't want to use "BENCHMARK_WITH_1_TO_MAX_THREADS", so be sure it's off.
#
# ###################################################################################################
# # QUICK START
# ###################################################################################################
#
# Set up your preferences in the files:
# - Tetrahedral_remeshing/benchmark/Tetrahedral_remeshing/benchmarking_config.h (remeshing options)
# - Tetrahedral_remeshing/benchmark/Tetrahedral_remeshing/concurrent_mesher_config.cfg (parallel options)
#
# Compile benchmark_tetrahedral_remeshing.cpp in a folder called `build-release`
# Go to Tetrahedral_remeshing/benchmark/Tetrahedral_remeshing/Charting
#
# Call:
#   sh benchmarking.sh $1 ... $12
# Arguments:
#   $1: directory containing the Tetrahedral_remeshing project (i.e. project root)
#   $2: directory containing the input data folder (e.g. "~/Data/Thingi10k/raw_meshes")
#   $3: directory containing the output results (recommended to be .../Tetrahedral_remeshing/benchmark/Tetrahedral_remeshing/Charting)
#   $4: target edge length
#   $5: number of iterations
#   $6: (unused, kept for compatibility)
#   $7: (unused, kept for compatibility)
#   $8: timeout value (in seconds)
#   $9: number of threads used (runs multiple data in parallel. Do NOT use if benching parallel mode!!)
#   $10: test identifier (e.g. hash of the last commit)
#   $11: test identifier of a previous test, to perform the difference with ${10} [optional]
#
# Find the result in the output_dir/charts.
#
# ###################################################################################################
# # SCRIPT DETAILS
# ###################################################################################################

# The script generates data by calling `run_benchmark.py`, which runs `benchmark_tetrahedral_remeshing.cpp`.
# The data is written in an XML file (see the macro CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA
# and the file benchmark.xml). One XML file per input data.
# The XML file is parsed by the scripts `generate_[robustness|performance|quality]_benchmark_charts.py`
# to generate individual charts for each metric, which are then merged together at the end
# of this script.
#
# To add a new metric to the benchmark, one must:
# - Add it into the benchmark.xml
# - Measure it (either in include, or in benchmark_tetrahedral_remeshing.cpp) and log it with the macro
#   CGAL_TETRAHEDRAL_REMESHING_SET_PERFORMANCE_DATA.
# - Parse it and plot it in the corresponding script (e.g. generate_robustness_benchmark_charts.py)
# - change the "tail -n" at the end of this file, and the "3x10" (layout of the final chart)

# ###################################################################################################
# # TODO
# ###################################################################################################

# - Make it work for other OS
# - Other metrics (edge radius ratio, Hausdorff distance)

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

# $1: directory containing the Tetrahedral_remeshing project
# $2: directory containing the output results 
# $3: target edge length
# $4: number of iterations
# $5: timeout value
# $6: test ID
# $7: input file path (passed by parallel)
set -x
function compute_benchmark_data() {
  filename=$(basename -- "$7")
  filename="${filename%.*}"

  echo "param #1 (directory containing the Tetrahedral_remeshing project): " $1
  echo "param #2 (directory containing the output results): " $2
  echo "param #3 (target edge length): " $3
  echo "param #4 (number of iterations): " $4 
  echo "param #5 (timeout value): " $5
  echo "param #6 (test ID): " $6
  echo "param #7 (the input file path): " $7  # Changed from ${9} to ${7}

  python3 "$1/benchmark/Tetrahedral_remeshing/Charting/run_benchmark.py" \
    --exec "$1/benchmark/Tetrahedral_remeshing/build-release/benchmark_tetrahedral_remeshing" \
    -i "$7" \
    --target_edge_length "$3" \
    --num_iterations "$4" \
    -t "$5" \
    --out "$2" \
    --test_ID "$6" #\
    #> "$2/logs/$6/$filename.log"
}

export -f compute_benchmark_data

echo "param #1 (directory containing the Tetrahedral_remeshing project): " $1
echo "param #2 (directory containing the input data folder): " $2
echo "param #3 (directory containing the output results): " $3
echo "param #4 (target edge length): " $4
echo "param #5 (number of iterations): " $5
echo "param #6 (timeout value): " $6
echo "param #7 (number of threads used): " $7
echo "param #8 (test ID): " $8
echo "param #9 (another test ID, to perform the difference with ${8}): " $9

echo "---------------------------------"

echo $# "argument(s) provided."

if [[ $# -lt 8 ]]; then
  echo "Expected at least 8 arguments, see list above."
  exit 1
fi

rm -rf $3/Robustness/logs/${8}
rm -rf $3/Robustness/results/${8}
rm -rf $3/Quality/logs/${8}
rm -rf $3/Quality/results/${8}
rm -rf $3/Performance/logs/${8}
rm -rf $3/Performance/results/${8}

mkdir -p $3/Robustness/logs/${8}
mkdir -p $3/Robustness/results/${8}
mkdir -p $3/Robustness/charts_data
mkdir -p $3/Robustness/charts
mkdir -p $3/Quality/logs/${8}
mkdir -p $3/Quality/results/${8}
mkdir -p $3/Quality/charts
mkdir -p $3/Performance/logs/${8}
mkdir -p $3/Performance/results/${8}
mkdir -p $3/Performance/charts
mkdir -p $3/logs/${8}
mkdir -p $3/charts

rm $3/Robustness/results/${8}/* # only for that one, others need to keep their results in case of diff

find $2 -mindepth 1 | parallel -j${7} compute_benchmark_data $1 $3 $4 $5 $6 $8 {}

python3 $1/benchmark/Tetrahedral_remeshing/Charting/Robustness/generate_robustness_benchmark_charts.py -i $3/Robustness/results/${8} -o $3/Robustness -c ${8}

if [ -z "${9}" ]; then
  python3 $1/benchmark/Tetrahedral_remeshing/Charting/Performance/generate_performance_benchmark_charts.py -i $3/Performance/results/${8} -o $3/Performance -c ${8};
else
  python3 $1/benchmark/Tetrahedral_remeshing/Charting/Performance/generate_performance_benchmark_charts.py -i $3/Performance/results/${8} -o $3/Performance -c ${8} -p $3/Performance/results/${9} -d ${9};
fi

if [ -z "${9}" ]; then
  python3 $1/benchmark/Tetrahedral_remeshing/Charting/Quality/generate_quality_benchmark_charts.py -i $3/Quality/results/${8} -o $3/Quality -c ${8};
else
  python3 $1/benchmark/Tetrahedral_remeshing/Charting/Quality/generate_quality_benchmark_charts.py -i $3/Quality/results/${8} -o $3/Quality -c ${8} -p $3/Quality/results/${9} -d ${9};
fi

charts_path="$(ls "$3/Robustness/charts"/* -dArt | tail -n 1) $(ls "$3/Performance/charts"/* -dArt | tail -n 10) $(ls "$3/Quality/charts"/* -dArt | tail -n 22)"

pdfjam --nup 3x11 $charts_path --outfile $3/charts/results_${8}_${9}_$(date '+%Y-%m-%d_%H:%M:%S').pdf
