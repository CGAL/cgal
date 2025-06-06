import yaml
import os
import glob
import subprocess
import argparse
import datetime
import sys
import json
import platform
import hashlib
from itertools import product
import pandas as pd
from generate_run_matrix import generate_run_matrix
from parse_results import validate_results_json

def get_latest_results_dir(base_dir):
    subdirs = [os.path.join(base_dir, d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    timestamped = [d for d in subdirs if os.path.basename(d)[:4].isdigit()]
    if not timestamped:
        raise FileNotFoundError("No timestamped results directories found.")
    return max(timestamped, key=os.path.getmtime)

parser = argparse.ArgumentParser(description="Run benchmarks as specified in a YAML config and run matrix.")
parser.add_argument('--config', type=str, default='config.yaml', help='Path to the YAML config file')
args = parser.parse_args()

with open(args.config) as f:
    config = yaml.safe_load(f)

pipeline_cfg = config.get('pipeline', {})
results_dir = pipeline_cfg['results_dir']
run_benchmarks = pipeline_cfg.get('run_benchmarks', True)
charts = pipeline_cfg.get('charts', True)
report = pipeline_cfg.get('report', True)

all_csvs = []
all_runs = []  # Collect all validated run dicts here

if run_benchmarks:
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    results_dir = os.path.join(results_dir, timestamp)
    os.makedirs(results_dir, exist_ok=True)

    # --- Save environment/version info to metadata.json ---
    metadata = {}
    metadata['timestamp'] = datetime.datetime.now().isoformat()
    metadata['python_version'] = platform.python_version()
    metadata['platform'] = platform.platform()
    metadata['config_path'] = os.path.abspath(args.config)
    with open(args.config, 'rb') as f:
        config_bytes = f.read()
        metadata['config_sha256'] = hashlib.sha256(config_bytes).hexdigest()
    try:
        repo_dir = os.path.dirname(os.path.abspath(__file__))
        git_hash = subprocess.check_output(['git', '-C', repo_dir, 'rev-parse', 'HEAD'], stderr=subprocess.DEVNULL).decode().strip()
        metadata['git_commit'] = git_hash
    except Exception:
        metadata['git_commit'] = None
    try:
        with open(os.path.join(results_dir, 'environment_metadata.json'), 'w') as f:
            json.dump(metadata, f, indent=2)
    except Exception as e:
        print(f"Warning: Could not write metadata.json: {e}")

# Get schema file from config or default to same folder as config.yaml
schema_file = config.get('schema_file', config.get('pipeline', {}).get('schema_file', None))
if not schema_file:
    schema_file = os.path.join(os.path.dirname(os.path.abspath(args.config)), 'results_schema.yaml')

for bench in config['benchmarks']:
    exec_path = bench['exec']
    input_dir = bench['input_dir']
    sweep = bench['sweep']
    benchmark_name = bench.get('name', 'benchmark')
    # Generate run matrix using shared function
    runs = generate_run_matrix(sweep)
    keys = list(sweep.keys())
    input_files = glob.glob(os.path.join(input_dir, '*'))
    bench_results_dir = os.path.join(results_dir, benchmark_name)
    os.makedirs(bench_results_dir, exist_ok=True)
    if run_benchmarks:
        for input_file in input_files:
            mesh_name = os.path.splitext(os.path.basename(input_file))[0]
            mesh_dir = os.path.join(bench_results_dir, mesh_name)
            os.makedirs(mesh_dir, exist_ok=True)
            for sweep_args in runs:
                out_suffix = "_".join(f"{k}{v}" for k, v in sweep_args.items())
                results_file = os.path.join(mesh_dir, f"{out_suffix}_results.json")
                # Create run_args as a dict to match cmd order
                run_args = {}
                run_args['input_mesh'] = input_file
                for k in keys:
                    run_args[k] = sweep_args[k]
                run_args['results_file'] = results_file
                cmd = [exec_path, input_file]
                for k in keys:
                    cmd.append(str(run_args[k]))
                cmd.append(results_file)
                result = subprocess.run(cmd)
                status = "success" if result.returncode == 0 else "failure"
                exec_metadata = None
                if 'exec_metadata' in bench:
                    exec_metadata_path = bench['exec_metadata']
                else:
                    exec_metadata_path = exec_path + '.metadata.json'
                if os.path.isfile(exec_metadata_path):
                    try:
                        with open(exec_metadata_path) as emf:
                            exec_metadata = json.load(emf)
                    except Exception as e:
                        exec_metadata = None
                try:
                    if os.path.isfile(results_file):
                        with open(results_file, 'r+') as f:
                            data = json.load(f)
                            if 'run_metadata' not in data:
                                data['run_metadata'] = {}
                            data['run_metadata']['input_arguments'] = run_args.copy()
                            data['run_metadata']['status'] = status
                            data['run_metadata']['timestamp'] = datetime.datetime.now().isoformat()
                            if exec_metadata:
                                data['run_metadata']['exec_metadata'] = exec_metadata
                            # Validate schema
                            try:
                                validate_results_json(data, results_file)
                            except Exception as ve:
                                raise ve
                            f.seek(0)
                            json.dump(data, f, indent=4)
                            f.truncate()
                            all_runs.append(data)  # Only add validated run
                    else:
                        raise FileNotFoundError(f"results.json not found: {results_file}")
                except Exception as e:
                    # Write fallback results.json
                    fallback = {
                        'run_metadata': {
                            'status': 'failure',
                            'input_arguments': run_args.copy(),
                            'timestamp': datetime.datetime.now().isoformat(),
                            'error': str(e)
                        }
                    }
                    if exec_metadata:
                        fallback['run_metadata']['exec_metadata'] = exec_metadata
                    with open(results_file, 'w') as f:
                        json.dump(fallback, f, indent=4)
                    print(f"[WARNING] Wrote fallback results.json for failed run: {results_file}\n  Reason: {e}")
        # After all runs, call parse_results.py to aggregate results for this benchmark
        csv_path = os.path.join(bench_results_dir, 'benchmark_results.csv')
        subprocess.run([
            sys.executable,
            os.path.join(os.path.dirname(__file__), 'parse_results.py'),
            bench_results_dir,
            '--output', csv_path,
            '--benchmark_name', benchmark_name,
            '--schema', schema_file
        ])
        all_csvs.append(csv_path)
    else:
        # If not running exec, just look for existing CSVs
        csv_path = os.path.join(bench_results_dir, 'benchmark_results.csv')
        if os.path.isfile(csv_path):
            all_csvs.append(csv_path)

# Aggregate all CSVs into one
if all_csvs:
    dfs = [pd.read_csv(csv, sep=';') for csv in all_csvs]
    agg_csv_path = os.path.join(results_dir, 'pipeline_results.csv')
    pd.concat(dfs, ignore_index=True).to_csv(agg_csv_path, index=False, sep=';')
    csv_path = agg_csv_path
else:
    csv_path = None

# Write the aggregated pipeline_results.json
agg_json_path = os.path.join(results_dir, 'pipeline_results.json')
with open(agg_json_path, 'w', encoding='utf-8') as f:
    json.dump(all_runs, f, indent=2)

# If 'charts' exists in config and charts is enabled, call generate_charts.py
if charts and 'charts' in config and config['charts'] and csv_path:
    subprocess.run([
        sys.executable,
        os.path.join(os.path.dirname(__file__), 'generate_charts.py'),
        args.config,
        csv_path
    ])

# If report is enabled, call generate_report.py
report_path = None
if report and 'charts' in config and config['charts'] and csv_path:
    subprocess.run([
        sys.executable,
        os.path.join(os.path.dirname(__file__), 'generate_report.py'),
        args.config,
        csv_path
    ])
    report_path = os.path.join(results_dir, 'report.html')

# --- Print summary at end ---
print("\n===== Pipeline Summary =====")
print(f"Results directory: {results_dir}")
print(f"Aggregated CSV: {csv_path if csv_path and os.path.isfile(csv_path) else 'Not generated'}")
charts_dir = os.path.join(results_dir, 'Charts')
if charts and 'charts' in config and config['charts']:
    print(f"Charts directory: {charts_dir if os.path.isdir(charts_dir) else 'Not generated'}")
else:
    print("Charts: Skipped")
if report and 'charts' in config and config['charts']:
    print(f"HTML report: {report_path if report_path and os.path.isfile(report_path) else 'Not generated'}")
else:
    print("Report: Skipped")
print("===========================\n") 