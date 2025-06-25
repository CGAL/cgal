import os
import glob
import yaml
import json
import numpy as np
from jinja2 import Environment
import plotly.graph_objects as go
import plotly.express as px
from plotly.offline import plot
import plotly.io as pio
import hashlib
import io
import sys
import argparse
from datetime import datetime
import statistics

def validate_config(config):
    """Validate the configuration file."""
    if not isinstance(config, dict):
        raise ValueError("Config must be a dictionary")
    
    required_keys = ['results_root']
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required config key: {key}")
    
    return True

def load_metadata(report_dir):
    path = os.path.join(report_dir, "environment_metadata.json")
    if not os.path.isfile(path):
        return {}
    with open(path, "r") as f:
        return json.load(f)

def collect_metric_paths(d, prefix=()):
    """Recursively collect all metric paths and their values from a nested dict."""
    paths = {}
    for k, v in d.items():
        if isinstance(v, dict):
            paths.update(collect_metric_paths(v, prefix + (k,)))
        else:
            paths[prefix + (k,)] = v
    return paths

def load_pipeline_results(results_dir, group_by_executable=True):
    run_map = {}
    # First check for pipeline_results.json in the top level directory
    pipeline_results_file = os.path.join(results_dir, "pipeline_results.json")
    if os.path.isfile(pipeline_results_file):
        with open(pipeline_results_file, "r") as f:
            results_list = json.load(f)
            # Handle list of results
            for data in results_list:
                # ==================================================================
                # CONFIGURABLE: Include or exclude exec_metadata in similarity comparison
                # - group_by_executable=True: Only compare runs with same exec_metadata
                #   (same executable, compile flags, optimizations, etc.)
                # - group_by_executable=False: Compare runs across different executables
                #   (useful for comparing different implementations)
                # ==================================================================
                input_args = data.get("run_metadata", {}).get("input_arguments", {})
                exec_metadata = data.get("run_metadata", {}).get("exec_metadata", {})
                
                # Create comparable args from all input arguments
                if input_args:
                    processed_items = []
                    for key, value in input_args.items():
                        # Process values to extract filename stems for file paths
                        if isinstance(value, str) and ('/' in value or '\\' in value or '.' in value):
                            # Likely a file path - extract stem (filename without extension)
                            filename = os.path.basename(value)
                            stem = os.path.splitext(filename)[0]
                            processed_items.append((key, stem))
                        else:
                            processed_items.append((key, value))
                    
                    comparable_args = dict(processed_items)
                else:
                    comparable_args = {}
                
                # Create run key - conditionally include exec_metadata based on config
                if group_by_executable:
                    # Include exec_metadata - separate groups for different executables
                    run_key_data = {
                        'input_args': comparable_args,
                        'exec_metadata': exec_metadata
                    }
                else:
                    # Exclude exec_metadata - allows comparison across different executables
                    run_key_data = {
                        'input_args': comparable_args
                    }
                
                run_key_str = str(run_key_data)
                
                if run_key_str not in run_map:
                    run_map[run_key_str] = []
                
                # Add benchmark source directory information to the result
                data["_benchmark_source_dir"] = os.path.basename(results_dir)
                
                run_map[run_key_str].append(data)
    return run_map

def compare_metrics(run_metrics_list):
    """Given a list of result objects (one per report), compare all numeric leaves in their metrics and resource usage."""
    all_paths = set()
    for result in run_metrics_list:
        # Extract metrics
        metrics = result.get("metrics", {})
        all_paths.update(collect_metric_paths(metrics).keys())
        
        # Extract resource usage data
        resource_usage = result.get("resource_usage", {})
        if resource_usage:
            resource_paths = collect_metric_paths(resource_usage, prefix=("resource_usage",))
            all_paths.update(resource_paths.keys())
    
    all_paths = sorted(all_paths)
    results = []
    for path in all_paths:
        values = []
        for result in run_metrics_list:
            # Try metrics first
            metrics = result.get("metrics", {})
            val = collect_metric_paths(metrics).get(path, None)
            
            # If not found in metrics, try resource usage
            if val is None and path[0] == "resource_usage":
                resource_usage = result.get("resource_usage", {})
                if resource_usage:
                    resource_paths = collect_metric_paths(resource_usage, prefix=("resource_usage",))
                    val = resource_paths.get(path, None)
            
            values.append(val)
        results.append((path, values))
    return results

def safe_extract_value_by_path(data, path):
    """Safely extract a value from nested dictionary structure using dot notation path."""
    try:
        current = data
        for key in path.split('.'):
            if isinstance(current, dict) and key in current:
                current = current[key]
            else:
                return None
        return current
    except:
        return None

def normalize_path_value(value):
    """Normalize file paths to handle different path separators consistently."""
    if isinstance(value, str) and ('/' in value or '\\' in value):
        # This looks like a file path - normalize it
        return os.path.normpath(value).replace('\\', '/')
    return value

def get_benchmark_id_from_result(result):
    """Extract benchmark identifier from a result using exec_metadata primarily."""
    benchmark_id = "Unknown_Benchmark"
    
    # First try to use benchmark_name if available (most reliable)
    benchmark_name = result.get("benchmark_name")
    if benchmark_name:
        benchmark_id = benchmark_name
    else:
        # Try to extract from exec_metadata (different executables = different benchmarks)
        exec_metadata = result.get("run_metadata", {}).get("exec_metadata", {})
        if exec_metadata:
            # Look for executable path or name in exec_metadata
            exec_path = exec_metadata.get("executable_path", exec_metadata.get("exec_path"))
            if exec_path:
                # Extract executable name without extension
                exec_name = os.path.splitext(os.path.basename(exec_path))[0]
                benchmark_id = exec_name
            else:
                # Fallback to first available exec_metadata field
                for key, value in exec_metadata.items():
                    if isinstance(value, str) and value:
                        benchmark_id = f"{key}_{value}"[:20]  # Truncate for readability
                        break
    
    # Add git commit info for version differentiation
    env_metadata = result.get("env_metadata", {})
    git_commit = env_metadata.get("git_commit")
    git_dirty = env_metadata.get("git_dirty", False)
    source_dir = result.get("_benchmark_source_dir", "")
    
    if git_commit:
        # Use short git hash (first 8 characters)
        git_hash = git_commit[:8]
        if git_dirty:
            git_hash += "*"
        
        # Combine benchmark name with git hash
        if source_dir:
            benchmark_id = f"{benchmark_id}_{git_hash}_{source_dir}"
        else:
            benchmark_id = f"{benchmark_id}_{git_hash}"
    else:
        # Fallback to source directory name if available
        if source_dir and benchmark_id == "Unknown_Benchmark":
            benchmark_id = source_dir
        elif source_dir:
            benchmark_id = f"{benchmark_id}_{source_dir}"
        
        # If still no identifier, use timestamp
        if benchmark_id == "Unknown_Benchmark":
            benchmark_id = env_metadata.get("timestamp", "Unknown_Benchmark")
            
            # If no env_metadata timestamp, try run_metadata timestamp
            if benchmark_id == "Unknown_Benchmark":
                run_timestamp = result.get("run_metadata", {}).get("timestamp", "Unknown")
                if run_timestamp != "Unknown":
                    # Truncate to minute to group runs from same benchmark execution
                    try:
                        benchmark_id = run_timestamp[:16]  # YYYY-MM-DDTHH:MM
                    except:
                        benchmark_id = run_timestamp
    
    return benchmark_id

def apply_chart_filter(result, chart_filter):
    """Apply filter conditions to a result. Returns True if result passes all filters."""
    if not chart_filter:
        return True
    
    for filter_path, expected_value in chart_filter.items():
        actual_value = safe_extract_value_by_path(result, filter_path)
        
        # Handle different comparison types
        if actual_value is None:
            return False
        
        # Convert to string for comparison to handle different types
        if str(actual_value) != str(expected_value):
            return False
    
    return True



def generate_charts(comparison_results, config):
    """Generate unified charts that support both time-series and parameter-based analysis."""
    charts = config.get("charts", [])
    if not charts:
        print("No charts defined in config.")
        return {}
        
    chart_files = {}
    
    # Process each chart configuration
    for chart_idx, chart_config in enumerate(charts):
        y_metric = chart_config.get("y")
        x_param = chart_config.get("x")  # If not specified, we'll use timestamp
        chart_filter = chart_config.get("filter", {})
        title = chart_config.get("title", f"{y_metric} vs {x_param or 'Time'}")
        style = chart_config.get("style", "o-")
        
        if not y_metric:
            print(f"Warning: Chart {chart_idx} missing required 'y' field. Skipping.")
            continue
            
        # If no x parameter specified, use timestamp (time-series chart)
        if not x_param:
            x_param = "run_metadata.timestamp"
            
        print(f"\nGenerating chart: {title}")
        print(f"Y metric: {y_metric}")
        print(f"X parameter: {x_param}")
        if chart_filter:
            print(f"Filters: {chart_filter}")
        
        # Collect data points for this chart, grouped by benchmark execution
        benchmark_data = {}
        total_results = 0
        filtered_results = 0
        
        for run_key, results_list in comparison_results.items():
            for result in results_list:
                total_results += 1
                
                # Apply chart filter first
                if not apply_chart_filter(result, chart_filter):
                    continue
                
                filtered_results += 1
                
                # Extract x value using the unified approach
                x_value = safe_extract_value_by_path(result, x_param)
                if x_value is None:
                    continue
                
                # Extract y value from metrics or resource usage using full path
                y_value = safe_extract_value_by_path(result, y_metric)
                if y_value is None or not isinstance(y_value, (int, float)):
                    continue
                    
                # Get benchmark identifier
                benchmark_id = get_benchmark_id_from_result(result)
                
                # Universal grouping: use run_key as the base grouping mechanism
                # Results are already grouped by input_args + exec_metadata
                data_key = (run_key, benchmark_id)
                
                if data_key not in benchmark_data:
                    benchmark_data[data_key] = {
                        'points': [],
                        'run_key': run_key,
                        'benchmark_id': benchmark_id
                    }
                
                benchmark_data[data_key]['points'].append({
                    'x': x_value,
                    'y': y_value
                })
        
        if chart_filter:
            print(f"Filter applied: {filtered_results}/{total_results} results passed filter")
        
        if not benchmark_data:
            print(f"Warning: No data points found for chart '{title}' after filtering")
            continue
        
        # Universal grouping rule: group by input_arguments + exec_metadata, 
        # excluding any variables used as x or y in the chart
        chart_series = {}
        
        # Extract variable names from x and y paths that should be excluded from grouping
        excluded_vars = set()
        
        # Extract variable from x path (e.g., "run_metadata.input_arguments.threads" -> exclude "threads")
        if x_param and 'input_arguments' in x_param:
            x_var_name = x_param.split('.')[-1]  # Get the last part (e.g., "threads")
            excluded_vars.add(x_var_name)
        
        # Y is typically a metric path, not an input variable, so usually no exclusion needed
        # But if y_metric contains input_arguments, we should exclude that too
        if 'input_arguments' in y_metric:
            y_var_name = y_metric.split('.')[-1]
            excluded_vars.add(y_var_name)
        
        for data_key, data_info in benchmark_data.items():
            run_key = data_info['run_key']
            benchmark_id = data_info['benchmark_id']
            
            # Create series ID based on input_args + exec_metadata, excluding x/y variables
            try:
                run_key_data = eval(run_key)
                input_args = run_key_data.get('input_args', {})
                exec_metadata = run_key_data.get('exec_metadata', {})
                
                # Filter input_args to exclude x/y variables
                filtered_input_args = {k: v for k, v in input_args.items() if k not in excluded_vars}
                
                # Create series identifier from benchmark + filtered args + exec_metadata
                series_parts = [benchmark_id]
                
                # Add filtered input arguments
                if filtered_input_args:
                    args_str = "_".join(f"{k}={v}" for k, v in sorted(filtered_input_args.items()))
                    series_parts.append(args_str)
                
                # Add exec_metadata (usually distinguishes different executables)
                if exec_metadata:
                    exec_str = "_".join(f"{k}={v}" for k, v in sorted(exec_metadata.items()))
                    series_parts.append(f"exec_{exec_str}")
                
                series_id = "_".join(series_parts)
                
            except Exception as e:
                # Fallback: use benchmark_id + run_key
                series_id = f"{benchmark_id}_{run_key}"
            
            if series_id not in chart_series:
                chart_series[series_id] = {
                    'points': [],
                    'run_key': run_key,
                    'benchmark_id': benchmark_id,
                    'label': benchmark_id  # Will be cleaned up later
                }
            
            chart_series[series_id]['points'].extend(data_info['points'])
        
        print(f"Found {len(chart_series)} series for this chart")
        
        if not chart_series:
            print(f"Warning: No series found for chart '{title}'")
            continue
        
        # Set chart title with filters
        chart_title = title
        if chart_filter:
            filter_str = ", ".join(f"{k}={v}" for k, v in chart_filter.items())
            chart_title = f"{title}\nFiltered by: {filter_str}"
        
        # Store the chart 
        chart_key = f"chart_{chart_idx}"
        if chart_key not in chart_files:
            chart_files[chart_key] = []
        
        # Create interactive Plotly chart
        print(f"Creating interactive Plotly chart for '{title}'")
        
        fig = go.Figure()
        
        # Add traces for each series
        for series_id, series_data in chart_series.items():
            points = series_data['points']
            
            # Sort points by x value for proper line connection
            sorted_points = sorted(points, key=lambda p: p['x'])
            x_vals = [p['x'] for p in sorted_points]
            y_vals = [p['y'] for p in sorted_points]
            
            # Create more compact hover text using plotly's customdata approach
            run_key = series_data['run_key']
            try:
                run_key_data = eval(run_key)
                input_args = run_key_data.get('input_args', {})
                exec_metadata = run_key_data.get('exec_metadata', {})
                
                customdata = []
                for i, (x_val, y_val) in enumerate(zip(x_vals, y_vals)):
                    # Get all input arguments, formatted nicely
                    config_parts = []
                    for key, value in input_args.items():
                        # For file paths, show only the filename or last folder
                        if isinstance(value, str) and ('/' in value or '\\' in value):
                            display_value = os.path.basename(value)
                        else:
                            display_value = str(value)
                        
                        config_parts.append(f"{key}:{display_value}")
                    
                    # Create config text with intelligent wrapping
                    if config_parts:
                        # If too many parameters, use line breaks for wrapping
                        if len(config_parts) > 3:
                            # Split into groups of 2-3 for better wrapping
                            config_text = ", ".join(config_parts[:3])
                            if len(config_parts) > 3:
                                config_text += "<br>     " + ", ".join(config_parts[3:])  # Add spacing for alignment
                        else:
                            config_text = ", ".join(config_parts)
                    else:
                        config_text = "Default"
                    
                    # Store data as list for custom hover template (benchmark, config)
                    customdata.append([
                        series_data['benchmark_id'],  # Already contains git hash info
                        config_text
                    ])
                
            except Exception as e:
                customdata = [[series_data['benchmark_id'], 'N/A'] for _ in x_vals]
            
            # Determine marker and line style
            if 'o' in style:
                mode = 'lines+markers'
                marker_size = 8
            elif '-' in style:
                mode = 'lines'
                marker_size = 6
            else:
                mode = 'markers'
                marker_size = 8
            
            # Use different markers for different series
            marker_symbol = 'circle'
            if 'd' in style:
                marker_symbol = 'diamond'
            elif 's' in style:
                marker_symbol = 'square'
            
            # Create a clean series name for the legend (shorter version)
            legend_name = series_id.replace('_', ' ').title()
            if len(legend_name) > 50:
                legend_name = legend_name[:47] + "..."
            
            fig.add_trace(go.Scatter(
                x=x_vals,
                y=y_vals,
                mode=mode,
                name=legend_name,
                customdata=customdata,
                hovertemplate=
                    f"<b>X:</b> %{{x}}<br>" +
                    f"<b>Y:</b> %{{y}}<br>" +
                    "<b>Benchmark:</b> %{customdata[0]}<br>" +
                    "<b>Input Arguments:</b> %{customdata[1]}" +
                    "<extra></extra>",
                marker=dict(
                    size=marker_size,
                    symbol=marker_symbol
                ),
                line=dict(width=2)
            ))
        
        # Update layout for better appearance
        fig.update_layout(
            title=chart_title,
            xaxis_title=x_param,
            yaxis_title=y_metric,
            hovermode='closest',
            height=600,
            legend=dict(
                orientation="v",
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.02,
                font=dict(size=10)  # Smaller font for multi-column layout
            ),
            margin=dict(r=180, l=60, t=80, b=60),  # Adjusted margins for multi-column
            template="plotly_white",  # Clean white background
            hoverlabel=dict(
                bgcolor="rgba(255,255,255,0.95)",
                bordercolor="rgba(0,0,0,0.3)",
                font_size=12,
                font_family="Arial",
                namelength=-1,  # Show full hover text
                align="left"
            )
        )
        
        # Convert to embeddable HTML div (not full HTML document)
        div_id = f"plotly-div-{hash(chart_title)}"
        
        # Create embeddable HTML with plotly CDN - use full width
        plotly_chart_html = f"""
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <div id="{div_id}" style="width:100%;height:550px;position:relative;margin:15px 0;"></div>
        <script>
            Plotly.newPlot('{div_id}', {fig.to_json()}, {{displayModeBar: true, responsive: true, displaylogo: false}});
        </script>
        """
        
        # Debug: Show what the plotly HTML looks like
        print(f"Plotly HTML starts with: {plotly_chart_html[:100]}...")
        
        # Store Plotly HTML directly
        chart_files[chart_key].append(plotly_chart_html)
        print(f"Successfully created interactive Plotly chart for '{title}' - HTML length: {len(plotly_chart_html)}")
        
        # Handle file saving if requested
        save_as = chart_config.get("save_as")
        if save_as:
            # Save as complete HTML file
            html_filename = save_as.replace('.png', '.html').replace('.jpg', '.html').replace('.jpeg', '.html')
            full_html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>{chart_title}</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
    <h1>{chart_title}</h1>
    <div id="{div_id}" style="width:100%;height:600px;"></div>
    <script>
        Plotly.newPlot('{div_id}', {fig.to_json()});
    </script>
</body>
</html>
"""
            with open(html_filename, 'w', encoding='utf-8') as f:
                f.write(full_html)
            print(f"Saved interactive chart: {html_filename}")
    
    return chart_files

def calculate_outliers(comparison_data, outlier_threshold):
    """
    Pre-calculate statistical information for each metric row including outliers, variance, etc.
    Returns a dict with statistical analysis for each value.
    """
    outlier_info = {}
    
    for run_key, metrics in comparison_data.items():
        outlier_info[run_key] = {}
        
        for path, values in metrics:
            # Extract numeric values only
            numeric_values = [v for v in values if isinstance(v, (int, float)) and v is not None]
            

            
            # Calculate statistics if we have enough data points
            row_outliers = [False] * len(values)  # Default: no outliers
            variance = 0
            coefficient_of_variation = 0
            percentiles = {}
            
            if len(numeric_values) > 0:
                try:
                    row_mean = statistics.mean(numeric_values)
                    
                    if len(numeric_values) > 1:
                        row_std = statistics.stdev(numeric_values)
                        variance = statistics.variance(numeric_values)
                        
                        # Coefficient of variation (CV) - relative variability
                        if row_mean != 0:
                            coefficient_of_variation = (row_std / abs(row_mean)) * 100
                        
                        # Calculate percentiles if we have enough data
                        if len(numeric_values) >= 3:
                            sorted_values = sorted(numeric_values)
                            percentiles = {
                                'p25': statistics.quantiles(sorted_values, n=4)[0] if len(sorted_values) >= 4 else sorted_values[0],
                                'p50': statistics.median(sorted_values),
                                'p75': statistics.quantiles(sorted_values, n=4)[2] if len(sorted_values) >= 4 else sorted_values[-1],
                                'min': min(sorted_values),
                                'max': max(sorted_values)
                            }
                        
                        # Detect outliers
                        if row_std > 0:  # Avoid division by zero
                            for i, value in enumerate(values):
                                if isinstance(value, (int, float)) and value is not None:
                                    z_score = abs(value - row_mean) / row_std
                                    row_outliers[i] = z_score > outlier_threshold
                    else:
                        row_std = 0
                        row_mean = numeric_values[0]
                        
                except Exception as e:
                    row_mean = statistics.mean(numeric_values) if numeric_values else 0
                    row_std = 0
            else:
                row_mean = 0
                row_std = 0
            
            # Store comprehensive statistical info for this path
            path_key = '.'.join(str(p) for p in path)
            outlier_info[run_key][path_key] = {
                'outliers': row_outliers,
                'mean': row_mean,
                'std': row_std,
                'variance': variance,
                'coefficient_of_variation': coefficient_of_variation,
                'percentiles': percentiles,
                'count': len(numeric_values),
                'outlier_count': sum(row_outliers)
            }
    
    return outlier_info

def main():
    parser = argparse.ArgumentParser(description="Compare benchmark results as specified in a YAML config.")
    parser.add_argument('--config', type=str, default='config.yaml', help='Path to the YAML config file')
    args = parser.parse_args()

    # Load and validate config
    try:
        with open(args.config) as f:
            config = yaml.safe_load(f)
        validate_config(config)
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)
    
    results_root = config.get("results_root", "")
    if not results_root:
        print("Error: 'results_root' not specified in config file")
        sys.exit(1)
    
    if not os.path.exists(results_root):
        print(f"Error: Results root directory not found: {results_root}")
        print("Please ensure the results_root path in the config is correct.")
        sys.exit(1)
    
    include_reports = config.get("include_reports", [])
    
    # Get analysis configuration (backward compatibility)
    analysis_config = config.get("analysis", {})
    significant_change_threshold = analysis_config.get("significant_change_threshold", 
                                                     config.get("significant_change_threshold", 0.0))
    outlier_threshold = analysis_config.get("outlier_threshold", 
                                           config.get("outlier_threshold", 3.0))
    group_by_executable = analysis_config.get("group_by_executable", 
                                             config.get("group_by_executable", True))
    
    # Get report configuration (backward compatibility)
    report_config = config.get("report", {})
    chart_grid_cols = report_config.get("chart_grid_cols", 
                                       config.get("report_chart_grid_cols", 2))
    
    charts = config.get("charts", [])

    # Find all report directories
    if include_reports:
        report_dirs = [os.path.join(results_root, d) for d in include_reports if os.path.isdir(os.path.join(results_root, d))]
        if not report_dirs:
            print(f"Warning: None of the specified report directories found in {results_root}")
            print(f"Specified directories: {include_reports}")
    else:
        try:
            all_dirs = [d for d in os.listdir(results_root) if os.path.isdir(os.path.join(results_root, d))]
            report_dirs = [os.path.join(results_root, d) for d in all_dirs]
            if not report_dirs:
                print(f"Warning: No subdirectories found in {results_root}")
        except OSError as e:
            print(f"Error accessing results directory: {e}")
            sys.exit(1)

    print(f"Found {len(report_dirs)} report directories to process:")
    for report_dir in report_dirs:
        print(f"  - {os.path.basename(report_dir)}")

    print(f"\nComparison mode: {'Same executable only' if group_by_executable else 'Cross-executable comparison'}")
    
    # Load pipeline results and metadata for each report
    report_metrics = {}
    report_timestamps = {}
    report_input_args = {}
    
    for report_dir in report_dirs:
        print(f"\nProcessing: {os.path.basename(report_dir)}")
        run_key_to_metrics_map = load_pipeline_results(report_dir, group_by_executable)
        if run_key_to_metrics_map:
            report_metrics[report_dir] = run_key_to_metrics_map
            metadata = load_metadata(report_dir)
            timestamp = metadata.get("timestamp", "unknown")
            report_timestamps[report_dir] = timestamp
            for run_key, data_list in run_key_to_metrics_map.items():
                if run_key not in report_input_args:
                    # Extract just the input_args part from the run_key for display
                    try:
                        run_key_data = eval(run_key)  # Convert string back to dict
                        input_args = run_key_data.get('input_args', {})
                        report_input_args[run_key] = input_args
                        print(f"  Debug: Extracted input_args from run_key: {input_args}")
                    except Exception as e:
                        # Fallback: extract from first result's run_metadata
                        fallback_args = data_list[0].get("run_metadata", {}).get("input_arguments", {})
                        report_input_args[run_key] = fallback_args
                        print(f"  Debug: Using fallback input_arguments: {fallback_args}")
                        print(f"  Debug: run_key parsing failed: {e}")
                        print(f"  Debug: run_key content: {run_key[:200]}...")
            print(f"  Loaded {len(run_key_to_metrics_map)} run configurations")
        else:
            print(f"  No results found in {report_dir}")

    if not report_metrics:
        print("\nError: No benchmark results found in any of the report directories.")
        print("Please ensure that:")
        print("1. The results_root path is correct")
        print("2. The directories contain pipeline_results.json files")
        print("3. The results were generated by the benchmark pipeline")
        sys.exit(1)

    # Compare metrics across reports
    comparison_results = {}
    for report_dir, run_key_to_metrics_map in report_metrics.items():
        for run_key, metrics_list in run_key_to_metrics_map.items():
            if run_key not in comparison_results:
                comparison_results[run_key] = []
            # Collect all results that have the same input arguments AND exec_metadata
            comparison_results[run_key].extend(metrics_list)

    # Generate comparison data - only for run_keys that have multiple results to compare
    comparison_data = {}
    for run_key, results_list in comparison_results.items():
        if len(results_list) > 1:
            # Only compare if we have multiple results with the same input arguments and exec_metadata
            comparison_data[run_key] = compare_metrics(results_list)
            try:
                run_key_data = eval(run_key)
                input_args = run_key_data.get('input_args', {})
                exec_metadata = run_key_data.get('exec_metadata', {})
                print(f"Comparing {len(results_list)} runs with input arguments: {input_args}")
                if exec_metadata:
                    print(f"  Exec metadata: {exec_metadata}")
                
                # Debug outlier detection
                metrics_comparison = compare_metrics(results_list)
                if metrics_comparison:
                    print(f"  Outlier detection debug (threshold={outlier_threshold}):")
                    for path, values in metrics_comparison[:2]:  # Show first 2 metrics
                        numeric_values = [v for v in values if isinstance(v, (int, float))]
                        if len(numeric_values) > 1:
                            mean_val = sum(numeric_values) / len(numeric_values)
                            variance = sum((v - mean_val) ** 2 for v in numeric_values) / (len(numeric_values) - 1)
                            std_dev = variance ** 0.5 if variance > 0 else 0
                            outliers = [v for v in numeric_values if std_dev > 0 and abs(v - mean_val) / std_dev > outlier_threshold]
                            print(f"    {'.'.join(path)}: {len(outliers)}/{len(numeric_values)} outliers detected")
                
                # Debug: Print some metric values for outlier detection testing
                metrics_comparison = compare_metrics(results_list)
                if metrics_comparison:
                    print(f"  Sample metrics for outlier detection:")
                    for path, values in metrics_comparison[:3]:  # Show first 3 metrics
                        numeric_values = [v for v in values if isinstance(v, (int, float))]
                        if len(numeric_values) > 1:
                            mean_val = sum(numeric_values) / len(numeric_values)
                            variance = sum((v - mean_val) ** 2 for v in numeric_values) / (len(numeric_values) - 1)
                            std_dev = variance ** 0.5 if variance > 0 else 0
                            max_z = max(abs(v - mean_val) / std_dev if std_dev > 0 else 0 for v in numeric_values)
                            print(f"    {'.'.join(path)}: values={numeric_values}, mean={mean_val:.3f}, std={std_dev:.3f}, max_z={max_z:.2f}")
            except:
                print(f"Comparing {len(results_list)} runs with run key: {run_key[:100]}...")
        else:
            print(f"Skipping comparison for run key - only {len(results_list)} result(s) found")

    # Generate charts
    chart_files = generate_charts(comparison_results, config)
    
    # Print chart generation info
    print("\nChart generation completed:")
    print(f"Number of charts generated: {len(chart_files)}")
    print("Using interactive Plotly charts with hover information and legends")
    print("Chart file keys:", list(chart_files.keys()) if chart_files else "None")
    
    # Generate HTML report
    output_file = os.path.join(results_root, "benchmark_comparison_report.html")
    try:
        env = Environment(
            autoescape=True,
            enable_async=False,
            keep_trailing_newline=False,
            extensions=['jinja2.ext.do', 'jinja2.ext.loopcontrols']
        )
        template = env.from_string("""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Benchmark Comparison</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; }
                table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
                th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
                th { background-color: #f2f2f2; }
                tr:nth-child(even) { background-color: #f9f9f9; }
                h2 { margin-top: 20px; }
                        .charts-grid { 
                            display: grid; 
                            grid-template-columns: repeat({{ chart_grid_cols }}, 1fr); 
                            gap: 20px; 
                            margin: 20px 0; 
                            width: 100%;
                        }
                        .chart-container { 
                            text-align: center;
                            overflow: visible;
                            position: relative;
                            padding: 20px;
                            width: 100%;
                            max-width: none;
                        }
                        .chart-grid { 
                            max-width: 100%; 
                            height: auto; 
                            border: 1px solid #ddd; 
                            border-radius: 8px; 
                        }
                        /* Plotly chart styling */
                        .chart-container .plotly-graph-div {
                            border: 1px solid #ddd;
                            border-radius: 8px;
                            margin: 0 auto;
                        }
                .run-section { border-bottom: 2px solid #ccc; padding-bottom: 30px; margin-bottom: 30px; }
                        .charts-section { margin-top: 40px; }
                        .filter-info { background-color: #e8f4fd; padding: 10px; margin: 10px 0; border-left: 4px solid #2196F3; }
                        .outlier { background-color: #ffebee !important; color: #c62828; font-weight: bold; }
                        .outlier-tooltip { position: relative; cursor: help; }
                        .outlier-tooltip:hover::after { 
                            content: 'Outlier (exceeds threshold standard deviations from mean)'; 
                            position: absolute; 
                            background: #333; 
                            color: white; 
                            padding: 5px; 
                            border-radius: 3px; 
                            font-size: 12px; 
                            white-space: nowrap; 
                            z-index: 1000; 
                            bottom: 100%; 
                            left: 50%; 
                            transform: translateX(-50%); 
                        }
                        /* Statistical styling */
                        .variance-col { background-color: #f0f8ff; font-weight: bold; }
                        .stats-tooltip { position: relative; cursor: help; }
                        .stats-tooltip:hover::after { 
                            content: attr(data-stats); 
                            position: absolute; 
                            background: #333; 
                            color: white; 
                            padding: 8px; 
                            border-radius: 3px; 
                            font-size: 11px; 
                            white-space: pre-line; 
                            z-index: 1000; 
                            bottom: 100%; 
                            left: 50%; 
                            transform: translateX(-50%); 
                            max-width: 300px;
                        }
                        .high-variance { background-color: #fff3cd; }
                        .low-variance { background-color: #d1f2eb; }
                        /* Test styles for debugging */
                        .test-outlier { background-color: #ff0000 !important; color: white !important; font-weight: bold !important; border: 3px solid #ff0000 !important; }
                        
                        /* Responsive design for smaller screens */
                        @media (max-width: 768px) {
                            .charts-grid { 
                                grid-template-columns: 1fr; 
                            }
                        }
            </style>
        </head>
        <body>
            <h1>Benchmark Comparison</h1>
            
                <!-- Configuration Summary -->
                <div style="background-color: #f8f9fa; padding: 15px; margin: 15px 0; border-radius: 5px;">
                    <h3>Configuration Summary</h3>
                    <p><strong>Results Root:</strong> {{ config.get('results_root', 'Not specified') }}</p>
                    <p><strong>Analysis:</strong> Significant change threshold: {{ config.get('analysis', {}).get('significant_change_threshold', config.get('significant_change_threshold', 0.0)) }}%, 
                       Outlier threshold: {{ config.get('analysis', {}).get('outlier_threshold', config.get('outlier_threshold', 3.0)) }} sigma</p>
                    <p><strong>Comparison Mode:</strong> {{ 'Same executable only' if config.get('analysis', {}).get('group_by_executable', config.get('group_by_executable', True)) else 'Cross-executable comparison' }}</p>
                    <p><strong>Charts:</strong> {{ config.get('charts', [])|length }} charts configured, {{ chart_grid_cols }} columns</p>
            </div>
            

                
                <!-- Charts Section -->
                {% if chart_files %}
                    <div class="charts-section">
                        <h2>Charts</h2>
                        <div class="charts-grid">
                            {% for key in chart_files.keys() %}
                            <div class="chart-container">
                                {% for chart in chart_files[key] %}
                                    {{ chart|safe }}
                                {% endfor %}
                            </div>
                        {% endfor %}
                        </div>
                    </div>
            {% endif %}
            
                <!-- Detailed Comparison Section -->
            {% if comparison_data %}
                    <h2>Detailed Comparisons</h2>
            {% for run_key, metrics in comparison_data.items() %}
                {% set run_results = comparison_results[run_key] %}
                <div class="run-section">
                            <h3>Run Configuration {{ loop.index }}</h3>
                            <h4>Input Arguments</h4>
                    <table>
                        <tr>
                            <th>Parameter</th>
                            <th>Value</th>
                        </tr>
                        {% if run_key in report_input_args %}
                            {% for key, value in report_input_args[run_key].items() %}
                                <tr>
                                    <td>{{ key }}</td>
                                    <td>{{ value }}</td>
                                </tr>
                            {% endfor %}
                        {% endif %}
                    </table>
                    
                            <!-- Show exec metadata -->
                            {% if not group_by_executable %}
                                <!-- When comparing across executables, show all unique exec_metadata -->
                                <h4>Executable Metadata (All Implementations)</h4>
                                {% set all_exec_metadata = {} %}
                                {% for result in run_results %}
                                    {% set exec_metadata = result.get('run_metadata', {}).get('exec_metadata', {}) %}
                                    {% set benchmark_id = get_benchmark_id_from_result(result) %}
                                    {% if exec_metadata and benchmark_id not in all_exec_metadata %}
                                        {% do all_exec_metadata.update({benchmark_id: exec_metadata}) %}
                                    {% endif %}
                                {% endfor %}
                                
                                {% if all_exec_metadata %}
                                    {% for benchmark_id, exec_metadata in all_exec_metadata.items() %}
                                        <h5>{{ benchmark_id }}</h5>
                                        <table style="margin-bottom: 15px;">
                                            <tr>
                                                <th>Property</th>
                                                <th>Value</th>
                                            </tr>
                                            {% for key, value in exec_metadata.items() %}
                                            <tr>
                                                <td>{{ key }}</td>
                                                <td>{{ value }}</td>
                                            </tr>
                                            {% endfor %}
                                        </table>
                                    {% endfor %}
                                {% endif %}
                            {% else %}
                                <!-- When comparing within same executable, show single exec_metadata -->
                                {% set first_result = run_results[0] if run_results else {} %}
                                {% set exec_metadata = first_result.get('run_metadata', {}).get('exec_metadata', {}) %}
                                {% if exec_metadata %}
                                <h4>Executable Metadata</h4>
                                <table>
                                    <tr>
                                        <th>Property</th>
                                        <th>Value</th>
                                    </tr>
                                    {% for key, value in exec_metadata.items() %}
                                    <tr>
                                        <td>{{ key }}</td>
                                        <td>{{ value }}</td>
                                    </tr>
                                    {% endfor %}
                                </table>
                                {% endif %}
                            {% endif %}
                    
                    {% set categories = {} %}
                    {% for path, values in metrics %}
                        {% set category = path[0] %}
                        {% if category not in categories %}
                            {% do categories.update({category: []}) %}
                        {% endif %}
                        {% do categories[category].append((path, values)) %}
                    {% endfor %}
                    {% for category, category_metrics in categories.items() %}
                            <h4>{{ category }} Metrics</h4>
                        <table>
                            <tr>
                                <th>Metric</th>
                                    <th>Variance (CV%)</th>
                                {% for result in run_results %}
                                    <th>{{ get_benchmark_id_from_result(result) }}</th>
                                {% endfor %}
                                    <th>Statistics</th>
                            </tr>
                            {% for path, values in category_metrics %}
                                    {% set path_key = '.'.join(path) %}
                                    {% set outlier_data = outlier_info.get(run_key, {}).get(path_key, {}) %}
                                    {% set variance = outlier_data.get('variance', 0) %}
                                    {% set cv = outlier_data.get('coefficient_of_variation', 0) %}
                                    {% set percentiles = outlier_data.get('percentiles', {}) %}
                                    {% set outlier_count = outlier_data.get('outlier_count', 0) %}
                                    
                                <tr>
                                    <td>{{ '.'.join(path[1:]) }}</td>
                                        <td class="variance-col stats-tooltip {% if cv > 20 %}high-variance{% elif cv < 5 %}low-variance{% endif %}"
                                            data-stats="Variance: {{ "%.6f"|format(variance) }}&#10;CV: {{ "%.1f"|format(cv) }}%&#10;{% if percentiles %}Min: {{ "%.6f"|format(percentiles.min) }}&#10;Q1: {{ "%.6f"|format(percentiles.p25) }}&#10;Median: {{ "%.6f"|format(percentiles.p50) }}&#10;Q3: {{ "%.6f"|format(percentiles.p75) }}&#10;Max: {{ "%.6f"|format(percentiles.max) }}{% endif %}">
                                            {{ "%.2f"|format(cv) }}%
                                        </td>
                                    {% for value in values %}
                                            {# Use pre-calculated outlier information #}
                                            {% set outliers_list = outlier_data.get('outliers', []) %}
                                            {% set is_outlier = outliers_list[loop.index0] if loop.index0 < outliers_list|length else False %}
                                            {% set row_mean = outlier_data.get('mean', 0) %}
                                            {% set row_std = outlier_data.get('std', 0) %}
                                            {% set row_count = outlier_data.get('count', 0) %}
                                            
                                            {# Calculate z-score for display #}
                                            {% set z_score = 0 %}
                                            {% if value is not none and value != "" and not (value is string) and row_std > 0 %}
                                                {% set z_score = ((value|float - row_mean)|abs / row_std) %}
                                            {% endif %}
                                            
                                            <td class="{% if is_outlier %}outlier outlier-tooltip{% endif %}" 
                                                title="Raw: {{ value }}, Mean: {{ "%.6f"|format(row_mean) }}, Std: {{ "%.6f"|format(row_std) }}, Z: {{ "%.2f"|format(z_score) }}, Threshold: {{ outlier_threshold }}, Count: {{ row_count }}, Outlier: {{ is_outlier }}">
                                                {% if value is not none and value != "" %}
                                                    {% if value is string %}
                                                        {{ value }}
                                                    {% else %}
                                                        {% set float_val = value|float %}
                                                        {% if float_val != 0 and (float_val|abs < 0.000001 or float_val|abs >= 1000000) %}
                                                            {{ "%.3e"|format(float_val) }}
                                                        {% else %}
                                                            {{ "%.6f"|format(float_val) }}
                                                        {% endif %}
                                                    {% endif %}
                                                {% else %}
                                                    N/A
                                                {% endif %}
                                            </td>
                                    {% endfor %}
                                        <td class="stats-tooltip" 
                                            data-stats="Mean: {{ "%.6f"|format(outlier_data.get('mean', 0)) }}&#10;Std Dev: {{ "%.6f"|format(outlier_data.get('std', 0)) }}&#10;Outliers: {{ outlier_count }}/{{ outlier_data.get('count', 0) }}">
                                            μ={{ "%.3f"|format(outlier_data.get('mean', 0)) }}
                                        </td>
                                </tr>
                            {% endfor %}
                        </table>
                    {% endfor %}
                </div>
            {% endfor %}
        {% else %}
                    <p>No comparison data available. This could mean:</p>
                    <ul>
                        <li>No runs with identical input arguments and exec_metadata were found across different benchmark executions</li>
                        <li>Only single runs were found for each configuration</li>
                        <li>The results contain different exec_metadata (different executables/compile options)</li>
                    </ul>
        {% endif %}
    </body>
    </html>
    """)
        
        # Calculate outlier information before rendering
        outlier_info = calculate_outliers(comparison_data, outlier_threshold)
        
        html = template.render(
            comparison_data=comparison_data, 
            comparison_results=comparison_results, 
            chart_files=chart_files, 
            report_input_args=report_input_args,
            config=config,
            outlier_threshold=outlier_threshold,
            outlier_info=outlier_info,
            chart_grid_cols=chart_grid_cols,
            group_by_executable=group_by_executable,
            get_benchmark_id_from_result=get_benchmark_id_from_result
        )
        
        with open(output_file, "w", encoding='utf-8') as f:
            f.write(html)
        
        print(f"\nHTML report generated successfully: {output_file}")
        
    except Exception as e:
        print(f"Error generating HTML report: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 