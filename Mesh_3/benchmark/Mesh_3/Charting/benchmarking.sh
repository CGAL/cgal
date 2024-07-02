#!/bin/bash

# ###################################################################################################
# # INTRODUCTION
# ###################################################################################################

# This script is meant to benchmark Mesh_3 over a directory of data files (e.g. Thingi10k/raw_meshes)
# with unique names. It checks a number of metrics in 3 categories: Robustness (did it finish without
# errors, and did it produce a valid result?), Performance (how long did it take to run?), and Quality
# (how good is the result mesh?).
#
# This script gives access to the mesh criteria for convenience, but all the other parameters
# are fixed and chosen two configuration files.
#
# If `benchmark_mesh_3.cpp` is in parallel mode, you should use a single thread.
# You also probably don't want to use "BENCHMARK_WITH_1_TO_MAX_THREADS", so be sure it's off.
#
# ###################################################################################################
# # QUICK START
# ###################################################################################################
#
# Set up your preferences in the files:
# - CGAL_root/Mesh_3/benchmark/Mesh_3/benchmarking_config.h (Mesh_3 options)
# - CGAL_root/Mesh_3/benchmark/Mesh_3/concurrent_mesher_config.cfg (parallel options)
#
# Compile benchmark_mesh_3.cpp in a folder called `build-release`
# Go to Mesh_3/benchmark/Mesh_3/Charting
#
# Call:
#   sh benchmarking.sh $1 ... $12
# Arguments:
#   $1: directory containing the Mesh_3 project (i.e. CGAL_root)
#   $2: directory containing the input data folder (e.g. "~/Data/Thingi10k/raw_meshes")
#   $3: directory containing the output results (recommended to be CGAL_root/Mesh_3/benchmark/Mesh_3/Charting)
#   $4: facet size
#   $5: facet distance
#   $6: facet angle
#   $7: cell size
#   $8: cell radius-edge ratio
#   $9: timeout value (in seconds)
#   $10: number of threads used (runs multiple data in parallel. Do NOT use if benching parallel mode!!)
#   $11: test identifier (e.g. hash of the last commit)
#   $12: test identifier of a previous test, to perform the difference with ${11} [optional]
#
# Find the result in the output_dir/charts.
#
# ###################################################################################################
# # SCRIPT DETAILS
# ###################################################################################################

# The script generates data by calling `run_benchmark.py`, which runs `benchmark_mesh_3.cpp`.
# The data is written in an XML file (see the macro CGAL_MESH_3_SET_PERFORMANCE_DATA
# and the file benchmark.xml). One XML file per input data.
# The XML file is parsed by the scripts `generate_[robustness|performance|quality]_benchmark_charts.py`
# to generate individual charts for each metric, which are then merged together at the end
# of this script.
#
# To add a new metric to the benchmark, one must:
# - Add it into the benchmark.xml
# - Measure it (either in Mesh_3/include, or in benchmark_mesh_3.cpp) and log it with the macro
#   CGAL_MESH_3_SET_PERFORMANCE_DATA.
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

# $1: directory containing the Mesh_3 project
# $2: directory containing the output results
# $3: facet size
# $4: facet distance
# $5: facet angle
# $6: cell size
# $7: cell radius-edge ratio
# $8: timeout value (in seconds)
# $9: test identifier (e.g. hash of the last commit)
# $10: the input file path
function compute_benchmark_data() {
  filename=$(basename -- "${10}")
  filename="${filename%.*}"

  # print filename:
  echo "param #1 (directory containing the Mesh_3 project): " $1
  echo "param #2 (directory containing the output results): " $2
  echo "param #3 (facet size): " $3
  echo "param #4 (facet distance): " $4
  echo "param #5 (facet angle): " $5
  echo "param #6 (cell size): " $6
  echo "param #7 (cell radius-edge ratio): " $7
  echo "param #8 (timeout value): " $8
  echo "param #9 (test ID): " $9
  echo "param #10 (the input file path): " ${10}

  python3 $1/Mesh_3/benchmark/Mesh_3/Charting/run_benchmark.py \
   --exec $1/Mesh_3/benchmark/Mesh_3/build-release/benchmark_mesh_3 \
   -i ${10} \
   --facet_size $3 --facet_approx $4 --facet_angle $5 \
   --cell_size $6 --cell_shape $7 -t $8 \
   --out $2 \
   --test_ID $9 \
  > $2/logs/$9/$filename.log

}
export -f compute_benchmark_data

echo "param #1 (directory containing the Mesh_3 project): " $1
echo "param #2 (directory containing the input data folder): " $2
echo "param #3 (directory containing the output results): " $3
echo "param #4 (facet size): " $4
echo "param #5 (facet distance): " $5
echo "param #6 (facet angle): " $6
echo "param #7 (cell size): " $7
echo "param #8 (cell radius-edge ratio): " $8
echo "param #9 (timeout value): " $9
echo "param #10 (number of threads used): " ${10}
echo "param #11 (test ID): " ${11}
echo "param #12 (another test ID, to perform the difference with ${11}): " ${12}

echo "---------------------------------"

echo $# "argument(s) provided."

if [[ $# -lt 11 ]]; then
  echo "Expected 12 arguments, see list above."
  exit 1
fi

# just for convenience: remove everything
# rm $3/*.mesh
# rm -rf $3/Robustness/logs
# rm -rf $3/Robustness/results
# rm -rf $3/Robustness/charts_data
# rm -rf $3/Robustness/charts
# rm -rf $3/Quality/logs
# rm -rf $3/Quality/results
# rm -rf $3/Quality/charts
# rm -rf $3/Performance/logs
# rm -rf $3/Performance/results
# rm -rf $3/Performance/charts
# rm -rf $3/logs
# rm -rf $3/charts
# exit

rm -rf $3/Robustness/logs/${11}
rm -rf $3/Robustness/results/${11}
rm -rf $3/Quality/logs/${11}
rm -rf $3/Quality/results/${11}
rm -rf $3/Performance/logs/${11}
rm -rf $3/Performance/results/${11}

mkdir -p $3/Robustness/logs/${11}
mkdir -p $3/Robustness/results/${11}
mkdir -p $3/Robustness/charts_data
mkdir -p $3/Robustness/charts
mkdir -p $3/Quality/logs/${11}
mkdir -p $3/Quality/results/${11}
mkdir -p $3/Quality/charts
mkdir -p $3/Performance/logs/${11}
mkdir -p $3/Performance/results/${11}
mkdir -p $3/Performance/charts
mkdir -p $3/logs/${11}
mkdir -p $3/charts

rm $3/Robustness/results/${11}/* # only for that one, others need to keep their results in case of diff

find $2 -mindepth 1 | parallel -j${10} compute_benchmark_data $1 $3 $4 $5 $6 $7 $8 $9 ${11} :::

python3 $1/Mesh_3/benchmark/Mesh_3/Charting/Robustness/generate_robustness_benchmark_charts.py -i $3/Robustness/results/${11} -o $3/Robustness -c ${11}

if [ -z "${12}" ]; then
  python3 $1/Mesh_3/benchmark/Mesh_3/Charting/Performance/generate_performance_benchmark_charts.py -i $3/Performance/results/${11} -o $3/Performance -c ${11};
else
  python3 $1/Mesh_3/benchmark/Mesh_3/Charting/Performance/generate_performance_benchmark_charts.py -i $3/Performance/results/${11} -o $3/Performance -c ${11} -p $3/Performance/results/${12} -d ${12};
fi

if [ -z "${12}" ]; then
  python3 $1/Mesh_3/benchmark/Mesh_3/Charting/Quality/generate_quality_benchmark_charts.py -i $3/Quality/results/${11} -o $3/Quality -c ${11};
else
  python3 $1/Mesh_3/benchmark/Mesh_3/Charting/Quality/generate_quality_benchmark_charts.py -i $3/Quality/results/${11} -o $3/Quality -c ${11} -p $3/Quality/results/${12} -d ${12};
fi

charts_path="$(ls "$3/Robustness/charts"/* -dArt | tail -n 1) $(ls "$3/Performance/charts"/* -dArt | tail -n 10) $(ls "$3/Quality/charts"/* -dArt | tail -n 22)"

pdfjam --nup 3x11 $charts_path --outfile $3/charts/results_${11}_${12}_$(date '+%Y-%m-%d_%H:%M:%S').pdf
