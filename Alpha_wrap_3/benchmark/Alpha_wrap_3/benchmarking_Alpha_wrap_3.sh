#!/bin/bash
CGAL_directory=$1
input_path=$2
output_dir=$3
alpha_value=$4
timeout_value=$5
virtual_thread=$6

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR=${BUILD_DIR:-"$SCRIPT_DIR/build-release"}

mkdir -p $output_dir/Robustness/results
mkdir -p $output_dir/Performance/results
mkdir -p $output_dir/Quality/results

function compute_benchmark_data() {
    local CGAL_directory="$1"
    local output_dir="$2"
    local alpha_value="$3"
    local timeout_value="$4"
    local build_dir="$5"
    local input_entry="$6"
    filename=$(basename -- "$input_entry")
    filename="${filename%.*}"

    echo "Running benchmarks for: $filename"

    python3 $CGAL_directory/Alpha_wrap_3/benchmark/Alpha_wrap_3/Robustness/compute_robustness_benchmark_data.py \
        -e $build_dir/robustness_benchmark -i "$input_entry" -a $alpha_value -t $timeout_value \
        > $output_dir/Robustness/results/$filename.log

    python3 $CGAL_directory/Alpha_wrap_3/benchmark/Alpha_wrap_3/Performance/compute_performance_benchmark_data.py \
        -e $build_dir/performance_benchmark -i "$input_entry" -a $alpha_value \
        > $output_dir/Performance/results/$filename.log

    python3 $CGAL_directory/Alpha_wrap_3/benchmark/Alpha_wrap_3/Quality/compute_quality_benchmark_data.py \
        -e $build_dir/quality_benchmark -i "$input_entry" -a $alpha_value \
        > $output_dir/Quality/results/$filename.log

    echo "Benchmarks completed for: $filename"
}
export -f compute_benchmark_data

if [ -d "$input_path" ]; then
  find "$input_path" -type f | parallel -j"$virtual_thread" compute_benchmark_data "$CGAL_directory" "$output_dir" "$alpha_value" "$timeout_value" "$BUILD_DIR" {}
else
  echo "Error: $input_path is not a valid directory."
  exit 1
fi


exit 0