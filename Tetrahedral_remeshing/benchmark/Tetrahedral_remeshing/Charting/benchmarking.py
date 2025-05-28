#!/usr/bin/env python3

import os
import sys
import glob
import subprocess
from datetime import datetime
from pathlib import Path
import multiprocessing
import argparse

def compute_benchmark_data(args):
    """
    Run benchmark for a single input file
    """
    project_dir, output_dir, target_edge_length, num_iterations, timeout, test_id, input_file = args
    filename = os.path.splitext(os.path.basename(input_file))[0]
    
    print(f"param #1 (project dir): {project_dir}")
    print(f"param #2 (output dir): {output_dir}")
    print(f"param #3 (target edge length): {target_edge_length}")
    print(f"param #4 (num iterations): {num_iterations}")
    print(f"param #5 (timeout): {timeout}")
    print(f"param #6 (test ID): {test_id}")
    print(f"param #7 (input file): {input_file}")
    
    benchmark_script = os.path.join(project_dir, "benchmark/Tetrahedral_remeshing/Charting/run_benchmark.py")
    benchmark_exec = os.path.join(project_dir, "benchmark/Tetrahedral_remeshing/build-release/benchmark_tetrahedral_remeshing")
    
    cmd = [
        "python3", benchmark_script,
        "--exec", benchmark_exec,
        "-i", input_file,
        "--target_edge_length", str(target_edge_length),
        "--num_iterations", str(num_iterations),
        "-t", str(timeout),
        "--out", output_dir,
        "--test_ID", test_id
    ]
    
    subprocess.run(cmd)

def main():
    parser = argparse.ArgumentParser(description="Benchmark Tetrahedral Remeshing")
    parser.add_argument("project_dir", help="Directory containing the Tetrahedral_remeshing project")
    parser.add_argument("input_dir", help="Directory containing input data folder")
    parser.add_argument("output_dir", help="Directory containing output results")
    parser.add_argument("target_edge_length", type=float, help="Target edge length")
    parser.add_argument("num_iterations", type=int, help="Number of iterations")
    parser.add_argument("timeout", type=int, help="Timeout value in seconds")
    parser.add_argument("num_threads", type=int, help="Number of parallel threads")
    parser.add_argument("test_id", help="Test identifier")
    parser.add_argument("diff_test_id", nargs='?', help="Previous test ID for comparison")
    
    args = parser.parse_args()
    
    # Create directory structure
    for category in ["Robustness", "Quality", "Performance"]:
        for subdir in ["logs", "results", "charts", "charts_data"]:
            path = os.path.join(args.output_dir, category, subdir)
            if subdir in ["logs", "results"]:
                path = os.path.join(path, args.test_id)
            os.makedirs(path, exist_ok=True)
    
    # Clean previous results
    results_dir = os.path.join(args.output_dir, "Robustness/results", args.test_id)
    for f in glob.glob(os.path.join(results_dir, "*")):
        os.remove(f)
    
    # Process input files in parallel
    input_files = glob.glob(os.path.join(args.input_dir, "*"))
    pool_args = [(args.project_dir, args.output_dir, args.target_edge_length, 
                 args.num_iterations, args.timeout, args.test_id, f) for f in input_files]
    
    with multiprocessing.Pool(args.num_threads) as pool:
        pool.map(compute_benchmark_data, pool_args)
    
    # Generate charts
    charts_script = {
        "Robustness": "Robustness/generate_robustness_benchmark_charts.py",
        "Performance": "Performance/generate_performance_benchmark_charts.py",
        "Quality": "Quality/generate_quality_benchmark_charts.py"
    }
    
    for category, script in charts_script.items():
        script_path = os.path.join(args.project_dir, "benchmark/Tetrahedral_remeshing/Charting", script)
        cmd = [
            "python3", script_path,
            "-i", os.path.join(args.output_dir, f"{category}/results", args.test_id),
            "-o", os.path.join(args.output_dir, category),
            "-c", args.test_id
        ]
        
        if args.diff_test_id:
            cmd.extend([
                "-p", os.path.join(args.output_dir, f"{category}/results", args.diff_test_id),
                "-d", args.diff_test_id
            ])
        
        subprocess.run(cmd)
    
    # Combine charts
    timestamp = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
    charts = []
    charts.extend(sorted(glob.glob(f"{args.output_dir}/Robustness/charts/*"))[-1:])
    charts.extend(sorted(glob.glob(f"{args.output_dir}/Performance/charts/*"))[-10:])
    charts.extend(sorted(glob.glob(f"{args.output_dir}/Quality/charts/*"))[-22:])
    
    output_pdf = f"{args.output_dir}/charts/results_{args.test_id}_{args.diff_test_id or ''}_{timestamp}.pdf"
    subprocess.run(["pdfjam", "--nup", "3x11"] + charts + ["--outfile", output_pdf])

if __name__ == "__main__":
    main()