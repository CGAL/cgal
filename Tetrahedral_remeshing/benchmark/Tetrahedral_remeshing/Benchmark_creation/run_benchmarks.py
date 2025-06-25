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
import logging
import time
from pathlib import Path
from generate_run_matrix import generate_run_matrix
from parse_results import validate_results_json


def flatten_dict(d, parent_key='', sep='.'):
    """Recursively flattens a nested dictionary with dot notation."""
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def setup_logging(results_dir):
    """Setup logging to both console and file."""
    log_file = os.path.join(results_dir, 'benchmark_run.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

def get_latest_results_dir(base_dir):
    """Get the most recent timestamped results directory."""
    if not os.path.exists(base_dir):
        raise FileNotFoundError(f"Base results directory does not exist: {base_dir}")
    
    subdirs = [os.path.join(base_dir, d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    timestamped = [d for d in subdirs if os.path.basename(d)[:4].isdigit()]
    if not timestamped:
        raise FileNotFoundError("No timestamped results directories found.")
    return max(timestamped, key=os.path.getmtime)

def collect_environment_metadata(config_path):
    """Collect environment and version information."""
    env_metadata = {
        'python_version': platform.python_version(),
        'platform': platform.platform(),
        'python_executable': sys.executable,
        'working_directory': os.getcwd()
    }
    
    # Config file hash
    try:
        with open(config_path, 'rb') as f:
            config_bytes = f.read()
            env_metadata['config_sha256'] = hashlib.sha256(config_bytes).hexdigest()
    except Exception as e:
        logging.warning(f"Could not compute config hash: {e}")
        env_metadata['config_sha256'] = None
    
    # Git information
    try:
        repo_dir = os.path.dirname(os.path.abspath(__file__))
        git_hash = subprocess.check_output(['git', '-C', repo_dir, 'rev-parse', 'HEAD'], 
                                         stderr=subprocess.DEVNULL).decode().strip()
        env_metadata['git_commit'] = git_hash
        
        # Check if there are uncommitted changes
        git_status = subprocess.check_output(['git', '-C', repo_dir, 'status', '--porcelain'], 
                                           stderr=subprocess.DEVNULL).decode().strip()
        env_metadata['git_dirty'] = bool(git_status)
    except Exception:
        env_metadata['git_commit'] = None
        env_metadata['git_dirty'] = None
    
    return env_metadata

def run_single_benchmark(cmd, results_file, timeout_seconds=3600, capture_output=False):
    """Run a single benchmark with timeout and error handling."""
    import subprocess
    
    try:
        # Run command directly with maximum performance
        if capture_output:
            result = subprocess.run(cmd, timeout=timeout_seconds, 
                                  capture_output=True, text=True)
            stdout, stderr = result.stdout, result.stderr
        else:
            # Use subprocess.run for maximum performance
            result = subprocess.run(cmd)
            stdout, stderr = '', ''
        
        return {
            'returncode': result.returncode,
            'stdout': stdout,
            'stderr': stderr,
            'execution_time': 0,  # No timing tracking
            'timeout': False,
            'resource_usage': {}
        }
        
    except subprocess.TimeoutExpired:
        return {
            'returncode': -1,
            'stdout': '',
            'stderr': f'Process timed out after {timeout_seconds} seconds',
            'execution_time': timeout_seconds,
            'timeout': True,
            'resource_usage': {}
        }
    except Exception as e:
        return {
            'returncode': -1,
            'stdout': '',
            'stderr': str(e),
            'execution_time': 0,
            'timeout': False,
            'resource_usage': {}
        }

def validate_config(config):
    """Validate the configuration file structure."""
    required_keys = ['benchmarks', 'pipeline']
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required key '{key}' in config")
    
    if not isinstance(config['benchmarks'], list):
        raise ValueError("'benchmarks' must be a list")
    
    for i, bench in enumerate(config['benchmarks']):
        if 'exec' not in bench:
            raise ValueError(f"Benchmark {i} missing 'exec' field")
        if 'sweep' not in bench:
            raise ValueError(f"Benchmark {i} missing 'sweep' field")
        
        # Check if executable exists
        if not os.path.isfile(bench['exec']):
            logging.warning(f"Executable not found: {bench['exec']}")

def main():
    parser = argparse.ArgumentParser(description="Run benchmarks as specified in a YAML config and run matrix.")
    parser.add_argument('--config', type=str, default='config.yaml', help='Path to the YAML config file')
    parser.add_argument('--timeout', type=int, default=3600, help='Timeout for each benchmark run in seconds')
    parser.add_argument('--continue-on-error', action='store_true', help='Continue running other benchmarks if one fails')

    args = parser.parse_args()

    # Load and validate config
    try:
        with open(args.config) as f:
            config = yaml.safe_load(f)
        validate_config(config)
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)

    pipeline_cfg = config.get('pipeline', {})
    results_dir = pipeline_cfg['results_dir']
    run_benchmarks = pipeline_cfg.get('run_benchmarks', True)
    export_csv = pipeline_cfg.get('export_csv', True)
    # Legacy: remove old charts reference
    # charts = pipeline_cfg.get('charts', True)
    report = pipeline_cfg.get('report', True)
    


    all_csvs = []
    all_runs = []
    env_metadata = {}

    if run_benchmarks:
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        results_dir = os.path.join(results_dir, timestamp)
        os.makedirs(results_dir, exist_ok=True)
        
        # Setup logging
        logger = setup_logging(results_dir)
        logger.info(f"Starting benchmark run with config: {args.config}")
        
        # Collect environment metadata
        env_metadata = collect_environment_metadata(args.config)
        logger.info(f"Environment: Python {env_metadata['python_version']}, {env_metadata['platform']}")

    # No schema file needed - pipeline is data-agnostic

    total_runs = 0
    successful_runs = 0
    failed_runs = 0

    for bench_idx, bench in enumerate(config['benchmarks']):
        exec_path = bench['exec']
        sweep = bench['sweep']
        benchmark_name = bench.get('name', f'benchmark_{bench_idx}')
        
        logger.info(f"Processing benchmark '{benchmark_name}' ({bench_idx + 1}/{len(config['benchmarks'])})")
        
        # Generate run matrix
        try:
            runs = generate_run_matrix(sweep)
            logger.info(f"Generated {len(runs)} runs for benchmark '{benchmark_name}'")
        except Exception as e:
            logger.error(f"Failed to generate run matrix for benchmark '{benchmark_name}': {e}")
            if not args.continue_on_error:
                sys.exit(1)
            continue

        keys = list(sweep.keys())
        bench_results_dir = os.path.join(results_dir, benchmark_name)
        os.makedirs(bench_results_dir, exist_ok=True)
        runs_dir = os.path.join(bench_results_dir, 'Runs')
        os.makedirs(runs_dir, exist_ok=True)

        if run_benchmarks:
            for run_idx, sweep_args in enumerate(runs):
                total_runs += 1
                logger.info(f"Running {run_idx + 1}/{len(runs)}: {sweep_args}")
                
                # Create filename suffix with automatic file path detection
                def is_file_path(value):
                    """Automatically detect if a value is a file path."""
                    if not isinstance(value, str):
                        return False
                    
                    # Check if it looks like a file path
                    if os.path.sep in str(value) or '/' in str(value) or '\\' in str(value):
                        return True
                    
                    # Check if it has a file extension
                    if '.' in str(value) and len(str(value).split('.')[-1]) <= 5:
                        return True
                    
                    # Check if the file actually exists
                    if os.path.isfile(str(value)):
                        return True
                    
                    return False
                
                # Use a simple run counter to avoid long filenames
                # This ensures unique filenames while keeping them short
                out_suffix = f"run_{run_idx:03d}"
                results_file = os.path.join(runs_dir, f"{out_suffix}_results.json")
                # Normalize path to ensure consistent separators
                results_file = os.path.normpath(results_file)
                
                # Ensure the directory exists before running the executable
                os.makedirs(os.path.dirname(results_file), exist_ok=True)
                
                # Create run_args
                run_args = sweep_args.copy()
                run_args['results_file'] = results_file
                
                # Build command
                cmd = [exec_path]
                for k in keys:
                    value = str(run_args[k])
                    # Normalize file paths to ensure consistent separators
                    if is_file_path(value):
                        value = os.path.normpath(value)
                    cmd.append(value)
                cmd.append(results_file)
                
                # Debug: Log the command being executed
                logger.info(f"Executing command: {' '.join(cmd)}")
                
                # Run benchmark
                run_result = run_single_benchmark(cmd, results_file, args.timeout)
                status = "success" if run_result['returncode'] == 0 else "failure"
                
                if status == "success":
                    successful_runs += 1
                else:
                    failed_runs += 1
                    logger.warning(f"Run failed with return code {run_result['returncode']}")
                    
                    # For failed runs, capture error output separately if needed for debugging
                    if logger.isEnabledFor(logging.DEBUG):  # Only capture detailed errors if debug logging is enabled
                        logger.debug("Capturing error output for failed run...")
                        error_result = run_single_benchmark(cmd, results_file, args.timeout, capture_output=True)
                        if error_result['stderr']:
                            logger.debug(f"Error output: {error_result['stderr']}")

                # Load exec metadata
                exec_metadata = None
                exec_metadata_path = bench.get('exec_metadata', exec_path + '.metadata.json')
                if os.path.isfile(exec_metadata_path):
                    try:
                        with open(exec_metadata_path) as emf:
                            exec_metadata = json.load(emf)
                    except Exception as e:
                        logger.warning(f"Could not load exec metadata: {e}")

                # Process results
                try:
                    if os.path.isfile(results_file):
                        with open(results_file, 'r+') as f:
                            data = json.load(f)
                            if 'run_metadata' not in data:
                                data['run_metadata'] = {}
                            
                            data['run_metadata']['input_arguments'] = sweep_args.copy()
                            data['run_metadata']['status'] = status
                            data['run_metadata']['timestamp'] = datetime.datetime.now().isoformat()
                            
                            if exec_metadata:
                                data['run_metadata']['exec_metadata'] = exec_metadata
                            
                            # Add environment metadata
                            data['env_metadata'] = env_metadata
                            
                            # Add benchmark name for consistency (required by report generation)
                            data['benchmark_name'] = benchmark_name
                            

                            
                            # No schema validation needed - pipeline is data-agnostic
                            
                            f.seek(0)
                            json.dump(data, f, indent=4)
                            f.truncate()
                            all_runs.append(data)
                    else:
                        raise FileNotFoundError(f"results.json not found: {results_file}")
                        
                except Exception as e:
                    # Write fallback results.json
                    fallback = {
                        'run_metadata': {
                            'status': 'failure',
                            'input_arguments': run_args.copy(),
                            'timestamp': datetime.datetime.now().isoformat(),
                            'error': str(e),
                            'stderr': run_result.get('stderr', ''),
                            'timeout': run_result.get('timeout', False)
                        },
                        'env_metadata': env_metadata
                    }
                    if exec_metadata:
                        fallback['run_metadata']['exec_metadata'] = exec_metadata
                    
                    # Add benchmark name for consistency
                    fallback['benchmark_name'] = benchmark_name
                    

                    
                    with open(results_file, 'w') as f:
                        json.dump(fallback, f, indent=4)
                    logger.warning(f"Wrote fallback results.json for failed run: {results_file}")
                    logger.warning(f"Reason: {e}")

                    # Add fallback to all_runs so report can be generated even for failed runs
                    all_runs.append(fallback)

            # Parse results for this benchmark (only if CSV export is enabled)
            if export_csv:
                csv_path = os.path.join(bench_results_dir, 'benchmark_results.csv')
                try:
                    subprocess.run([
                        sys.executable,
                        os.path.join(os.path.dirname(__file__), 'parse_results.py'),
                        bench_results_dir,
                        '--output', csv_path,
                        '--benchmark_name', benchmark_name
                    ], check=True)
                    all_csvs.append(csv_path)
                except subprocess.CalledProcessError as e:
                    logger.error(f"Failed to parse results for benchmark '{benchmark_name}': {e}")
        else:
            # Look for existing CSVs (only if CSV export is enabled)
            if export_csv:
                csv_path = os.path.join(bench_results_dir, 'benchmark_results.csv')
                if os.path.isfile(csv_path):
                    all_csvs.append(csv_path)

    # Load existing results if not running benchmarks
    if not run_benchmarks:
        try:
            results_dir = get_latest_results_dir(pipeline_cfg['results_dir'])
            existing_pipeline_json = os.path.join(results_dir, 'pipeline_results.json')
            if os.path.isfile(existing_pipeline_json):
                with open(existing_pipeline_json, 'r', encoding='utf-8') as f:
                    existing_data = json.load(f)
                    if isinstance(existing_data, list):
                        all_runs = existing_data
        except Exception as e:
            print(f"Warning: Could not load existing pipeline_results.json: {e}")

    # Aggregate CSVs (only if CSV export is enabled)
    csv_path = None
    if export_csv and all_csvs:
        try:
            dfs = [pd.read_csv(csv, sep=';') for csv in all_csvs]
            agg_csv_path = os.path.join(results_dir, 'pipeline_results.csv')
            pd.concat(dfs, ignore_index=True).to_csv(agg_csv_path, index=False, sep=';')
            csv_path = agg_csv_path
        except Exception as e:
            logging.error(f"Failed to aggregate CSVs: {e}")

    # Write aggregated pipeline_results.json
    agg_json_path = os.path.join(results_dir, 'pipeline_results.json')
    try:
        with open(agg_json_path, 'w', encoding='utf-8') as f:
            json.dump(all_runs, f, indent=2)
    except Exception as e:
        logging.error(f"Failed to write pipeline_results.json: {e}")

    # Generate charts
    export_charts = pipeline_cfg.get('export_charts', pipeline_cfg.get('charts', False))  # Support both names
    charts_should_generate = 'charts' in config and config['charts'] and all_runs
    if charts_should_generate:
        # Always generate charts if defined in config, regardless of export_charts boolean
        try:
            subprocess.run([
                sys.executable,
                os.path.join(os.path.dirname(__file__), 'generate_charts.py'),
                args.config,
                agg_json_path
            ], check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to generate charts: {e}")
    elif export_charts:
        # If export_charts boolean is True but no charts defined, create empty Charts directory
        charts_dir = os.path.join(results_dir, 'Charts')
        os.makedirs(charts_dir, exist_ok=True)
        logging.info("Created Charts directory (export_charts enabled but no charts configured)")
    else:
        logging.info("Charts generation skipped")

    # Print status in summary
    charts_dir = os.path.join(results_dir, 'Charts')
    if charts_should_generate:
        print(f"Charts directory: {charts_dir if os.path.isdir(charts_dir) else 'Generation failed'}")
    elif export_charts:
        print(f"Charts directory: {charts_dir} (created but no charts configured)")
    else:
        print("Charts: Disabled")

    # Generate report from JSON data
    report_path = None
    if report and all_runs:
        try:
            subprocess.run([
                sys.executable,
                os.path.join(os.path.dirname(__file__), 'generate_report.py'),
                args.config,
                results_dir
            ], check=True)
            report_path = os.path.join(results_dir, 'report.html')
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to generate report: {e}")

    # Print summary
    print("\n===== Pipeline Summary =====")
    print(f"Results directory: {results_dir}")
    if run_benchmarks:
        print(f"Total runs: {total_runs}")
        print(f"Successful: {successful_runs}")
        print(f"Failed: {failed_runs}")
        print(f"Success rate: {successful_runs/total_runs*100:.1f}%" if total_runs > 0 else "No runs")
    if export_csv:
        print(f"Aggregated CSV: {csv_path if csv_path and os.path.isfile(csv_path) else 'Not generated'}")
    else:
        print("CSV export: Disabled")
    if report:
        print(f"HTML report: {report_path if report_path and os.path.isfile(report_path) else 'Not generated'}")
    else:
        print("Report: Skipped")
    print("===========================\n")

if __name__ == "__main__":
    main() 