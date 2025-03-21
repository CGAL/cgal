#!/bin/bash
CGAL_directory=$1
input_path=$2
output_dir=$3
alpha_value=$4
timeout_value=$5
virtual_thread=$6
commit_hash=$7
mode=${8:-""}

SCRIPT_DIR="$CGAL_directory/Alpha_wrap_3/benchmark/Alpha_wrap_3"
BUILD_DIR=${BUILD_DIR:-"/app/build/Alpha_wrap_3"}

mkdir -p $output_dir/Robustness/results/$commit_hash
mkdir -p $output_dir/Performance/results/$commit_hash
mkdir -p $output_dir/Quality/results/$commit_hash

function compute_benchmark_data() {
    local input_file="$1"
    filename=$(basename -- "$input_file")
    filename="${filename%.*}"
    
    echo "Running benchmarks for: $filename"
    
    python3 $SCRIPT_DIR/Robustness/compute_robustness_benchmark_data.py \
    -e $BUILD_DIR/robustness_benchmark -i "$input_file" -a $alpha_value -t $timeout_value \
    > $output_dir/Robustness/results/$commit_hash/$filename.log
    
    python3 $SCRIPT_DIR/Performance/compute_performance_benchmark_data.py \
    -e $BUILD_DIR/performance_benchmark -i "$input_file" -a $alpha_value \
    > $output_dir/Performance/results/$commit_hash/$filename.log
    
    python3 $SCRIPT_DIR/Quality/compute_quality_benchmark_data.py \
    -e $BUILD_DIR/quality_benchmark -i "$input_file" -a $alpha_value \
    > $output_dir/Quality/results/$commit_hash/$filename.log
    
    echo "Benchmarks completed for: $filename"
}

case "$mode" in
  "--single-file")
    if [ -f "$input_path" ]; then
      compute_benchmark_data "$input_path"
    else
      echo "Error: $input_path is not a valid file."
      exit 1
    fi
    ;;
    
  *)
    if [ -d "$input_path" ]; then
      find "$input_path" -type f | parallel -j"$virtual_thread" compute_benchmark_data {}
    else
      echo "Error: $input_path is not a valid directory."
      exit 1
    fi
    ;;
esac

exit 0