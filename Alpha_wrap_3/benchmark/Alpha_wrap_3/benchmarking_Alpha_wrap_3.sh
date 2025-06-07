#!/bin/bash
CGAL_directory=$1
input_path=$2
output_dir=$3
timeout_value=$4
virtual_thread=$5
alpha_value=$6

Benchmark_output_dir="$output_dir/Logs/Alpha_wrap_3"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR=${BUILD_DIR:-"$SCRIPT_DIR/build-release"}

mkdir -p $Benchmark_output_dir/Robustness/results
mkdir -p $Benchmark_output_dir/Performance/results
mkdir -p $Benchmark_output_dir/Quality/results

function compute_benchmark_data() {
    local CGAL_directory="$1"
    local Benchmark_output_dir="$2"
    local alpha_value="$3"
    local timeout_value="$4"
    local build_dir="$5"
    local input_entry="$6"
    filename=$(basename -- "$input_entry")
    filename="${filename%.*}"

    echo "Running benchmarks for: $filename"

    python3 $CGAL_directory/Alpha_wrap_3/benchmark/Alpha_wrap_3/Robustness/compute_robustness_benchmark_data.py \
        -e $build_dir/robustness_benchmark -i "$input_entry" -a $alpha_value -t $timeout_value \
        > $Benchmark_output_dir/Robustness/results/$filename.log

    python3 $CGAL_directory/Alpha_wrap_3/benchmark/Alpha_wrap_3/Performance/compute_performance_benchmark_data.py \
        -e $build_dir/performance_benchmark -i "$input_entry" -a $alpha_value \
        > $Benchmark_output_dir/Performance/results/$filename.log

    python3 $CGAL_directory/Alpha_wrap_3/benchmark/Alpha_wrap_3/Quality/compute_quality_benchmark_data.py \
        -e $build_dir/quality_benchmark -i "$input_entry" -a $alpha_value \
        > $Benchmark_output_dir/Quality/results/$filename.log

    echo "Benchmarks completed for: $filename"
}
export -f compute_benchmark_data

if [ ! -f "$BUILD_DIR/robustness_benchmark" ]; then
    echo "-> Compiling benchmarks in $BUILD_DIR"
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"

    if [ -f "$CGAL_directory/Alpha_wrap_3/benchmark/Alpha_wrap_3/CMakeLists.txt" ]; then
        cmake "$CGAL_directory/Alpha_wrap_3/benchmark/Alpha_wrap_3" \
            -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CGAL_directory"
        make -j"$(nproc)"

        if [ $? -ne 0 ]; then
            echo "Compilation failed."
            exit 1
        fi
        echo "Compilation done."
    else
        echo "No CMakeLists.txt found. Skipping compilation."
    fi
    cd "$SCRIPT_DIR"
else
    echo "Benchmarks already compiled in $BUILD_DIR"
fi

valid_extensions=( -iname '*.off' -o -iname '*.obj' -o -iname '*.ply' -o -iname '*.stl' -o -iname '*.STL' -o -iname '*.ts' -o -iname '*.vtp' )

if [ -d "$input_path" ]; then
  find "$input_path" -type f \( "${valid_extensions[@]}" \) | \
    parallel -j"$virtual_thread" compute_benchmark_data "$CGAL_directory" "$Benchmark_output_dir" "$alpha_value" "$timeout_value" "$BUILD_DIR" {}
else
  echo "Error: $input_path is not a valid directory."
  exit 1
fi

python3 "$CGAL_directory/Alpha_wrap_3/benchmark/Alpha_wrap_3/alpha_wrap_results_parser.py" \
    --meshes-dir "$input_path" \
    --benchmark-results-dir "$Benchmark_output_dir" \
    --summary-json-path "$output_dir"

exit 0