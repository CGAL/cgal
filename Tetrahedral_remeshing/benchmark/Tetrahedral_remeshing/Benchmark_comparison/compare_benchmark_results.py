import os
import glob
import yaml
import json
import numpy as np
from jinja2 import Environment
import matplotlib.pyplot as plt
import seaborn as sns
import hashlib
import base64
import io

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

def load_pipeline_results(results_dir):
    run_map = {}
    # First check for pipeline_results.json in the top level directory
    pipeline_results_file = os.path.join(results_dir, "pipeline_results.json")
    if os.path.isfile(pipeline_results_file):
        with open(pipeline_results_file, "r") as f:
            results_list = json.load(f)
            # Handle list of results
            for data in results_list:
                # ==================================================================
                # IMPORTANT ASSUMPTION: We assume that the LAST input argument is 
                # the JSON file path where results are saved. Since the output file
                # path is irrelevant for determining if two runs are "similar" 
                # (and thus comparable), we exclude it from the run_key computation.
                # This allows us to compare runs with identical parameters that 
                # saved their results to different files.
                # ==================================================================
                input_args = data.get("run_metadata", {}).get("input_arguments", {})
                if input_args:
                    # Convert to list of items, exclude the last one, then back to dict
                    items = list(input_args.items())
                    if len(items) > 1:
                        comparable_items = items[:-1]  # Exclude last argument
                    else:
                        comparable_items = list(items)  # Keep as-is if only one argument
                    
                    # Process values to extract filename stems for file paths
                    processed_items = []
                    for key, value in comparable_items:
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
                
                run_key_str = str(comparable_args)
                if run_key_str not in run_map:
                    run_map[run_key_str] = []
                run_map[run_key_str].append(data)
    return run_map

def compare_metrics(run_metrics_list):
    """Given a list of result objects (one per report), compare all numeric leaves in their metrics."""
    all_paths = set()
    for result in run_metrics_list:
        # Extract just the metrics part from each result
        metrics = result.get("metrics", {})
        all_paths.update(collect_metric_paths(metrics).keys())
    all_paths = sorted(all_paths)
    results = []
    for path in all_paths:
        values = []
        for result in run_metrics_list:
            # Extract just the metrics part from each result
            metrics = result.get("metrics", {})
            val = collect_metric_paths(metrics).get(path, None)
            values.append(val)
        results.append((path, values))
    return results

def generate_charts(comparison_data, comparison_results, config):
    charts = config.get("charts", [])
    if not charts:
        print("No charts defined in config.")
        return
    
    # Get grid configuration
    grid_cols = config.get("report_chart_grid_cols", 2)
    chart_files = {}
    
    for run_key, metrics in comparison_data.items():
        # Get the raw results for this run_key to access timestamps
        raw_results = comparison_results.get(run_key, [])
        
        # Create a grid of subplots for this run_key
        num_charts = len(charts)
        if num_charts == 0:
            continue
            
        # Calculate grid dimensions
        grid_rows = (num_charts + grid_cols - 1) // grid_cols  # Ceiling division
        
        # Create figure with subplots
        fig, axes = plt.subplots(grid_rows, grid_cols, figsize=(grid_cols * 6, grid_rows * 4))
        fig.suptitle(f'Benchmark Comparison - Run Configuration {list(comparison_data.keys()).index(run_key) + 1}', fontsize=14)
        
        # Handle single subplot case
        if num_charts == 1:
            axes = [axes]
        elif grid_rows == 1:
            axes = [axes] if grid_cols == 1 else axes
        else:
            axes = axes.flatten()
        
        chart_idx = 0
        
        for chart in charts:
            y_metric = chart.get("y")
            if not y_metric:
                print(f"Warning: Chart missing y component. Skipping.")
                continue
                
            if chart_idx >= len(axes):
                break
                
            ax = axes[chart_idx]
            x_metric = chart.get("x")
            title = chart.get("title", f"{y_metric}")
            
            x_values = []
            y_values = []
            
            # Extract data from each result
            for i, result in enumerate(raw_results):
                # Find x value
                if x_metric:
                    # User specified a custom x metric
                    x_val = None
                    for path, values in metrics:
                        if '.'.join(path) == x_metric:
                            x_val = values[i] if i < len(values) else None
                            break
                else:
                    # Default to timestamp from run metadata
                    timestamp = result.get("run_metadata", {}).get("timestamp", f"Run {i+1}")
                    x_val = timestamp
                
                # Find y value
                y_val = None
                for path, values in metrics:
                    if '.'.join(path) == y_metric:
                        y_val = values[i] if i < len(values) else None
                        break
                        
                if x_val is not None and y_val is not None:
                    x_values.append(x_val)
                    y_values.append(y_val)
            
            if x_values and y_values:
                style = chart.get("style", "o-")
                ax.plot(x_values, y_values, style)
            
            ax.set_title(title, fontsize=10)
            ax.set_xlabel(x_metric if x_metric else "Timestamp", fontsize=9)
            ax.set_ylabel(y_metric.split('.')[-1], fontsize=9)  # Use last part of metric name
            ax.tick_params(axis='both', which='major', labelsize=8)
            
            # Rotate x-axis labels if they are timestamps
            if not x_metric:
                ax.tick_params(axis='x', rotation=45)
            
            chart_idx += 1
        
        # Hide unused subplots
        for i in range(chart_idx, len(axes)):
            axes[i].set_visible(False)
        
        # Save chart to memory and encode as base64
        plt.tight_layout()
        buffer = io.BytesIO()
        plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
        buffer.seek(0)
        
        # Encode image as base64
        image_base64 = base64.b64encode(buffer.read()).decode('utf-8')
        plt.close()
        
        # Create data URI
        data_uri = f"data:image/png;base64,{image_base64}"
        
        # Also save to file for debugging/reference
        sanitized_run_key = str(run_key).replace('{', '').replace('}', '').replace(':', '_').replace(',', '_').replace(' ', '_')
        grid_filename = f"{sanitized_run_key}_chart_grid.png"
        
        # Save the chart grid again to file for reference
        fig_copy, axes_copy = plt.subplots(grid_rows, grid_cols, figsize=(grid_cols * 6, grid_rows * 4))
        fig_copy.suptitle(f'Benchmark Comparison - Run Configuration {list(comparison_data.keys()).index(run_key) + 1}', fontsize=14)
        
        # Recreate the plot for file saving (since we closed the original)
        if num_charts == 1:
            axes_copy = [axes_copy]
        elif grid_rows == 1:
            axes_copy = [axes_copy] if grid_cols == 1 else axes_copy
        else:
            axes_copy = axes_copy.flatten()
        
        chart_idx = 0
        for chart in charts:
            y_metric = chart.get("y")
            if not y_metric or chart_idx >= len(axes_copy):
                continue
                
            ax = axes_copy[chart_idx]
            x_metric = chart.get("x")
            title = chart.get("title", f"{y_metric}")
            
            x_values = []
            y_values = []
            for i, result in enumerate(raw_results):
                if x_metric:
                    x_val = None
                    for path, values in metrics:
                        if '.'.join(path) == x_metric:
                            x_val = values[i] if i < len(values) else None
                            break
                else:
                    timestamp = result.get("run_metadata", {}).get("timestamp", f"Run {i+1}")
                    x_val = timestamp
                y_val = None
                for path, values in metrics:
                    if '.'.join(path) == y_metric:
                        y_val = values[i] if i < len(values) else None
                        break
                if x_val is not None and y_val is not None:
                    x_values.append(x_val)
                    y_values.append(y_val)
            
            if x_values and y_values:
                style = chart.get("style", "o-")
                ax.plot(x_values, y_values, style)
            
            ax.set_title(title, fontsize=10)
            ax.set_xlabel(x_metric if x_metric else "Timestamp", fontsize=9)
            ax.set_ylabel(y_metric.split('.')[-1], fontsize=9)
            ax.tick_params(axis='both', which='major', labelsize=8)
            if not x_metric:
                ax.tick_params(axis='x', rotation=45)
            chart_idx += 1
        
        for i in range(chart_idx, len(axes_copy)):
            axes_copy[i].set_visible(False)
            
        plt.tight_layout()
        plt.savefig(grid_filename, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Chart grid saved as {grid_filename} and embedded in HTML")
        chart_files[run_key] = [data_uri]  # Store the data URI instead of filename
    
    return chart_files

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(script_dir, "benchmark_comparison_config.yaml")
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    results_root = config["results_root"]
    include_reports = config.get("include_reports", [])
    significant_change_threshold = config.get("significant_change_threshold", 0.0)
    outlier_threshold = config.get("outlier_threshold", 3.0)
    charts = config.get("charts", [])

    # Find all report directories
    if include_reports:
        report_dirs = [os.path.join(results_root, d) for d in include_reports if os.path.isdir(os.path.join(results_root, d))]
    else:
        report_dirs = [os.path.join(results_root, d) for d in os.listdir(results_root) if os.path.isdir(os.path.join(results_root, d))]

    # Load pipeline results and metadata for each report
    report_metrics = {}
    report_timestamps = {}
    report_input_args = {}
    for report_dir in report_dirs:
        run_key_to_metrics_map = load_pipeline_results(report_dir)
        if run_key_to_metrics_map:
            report_metrics[report_dir] = run_key_to_metrics_map
            metadata = load_metadata(report_dir)
            timestamp = metadata.get("timestamp", "unknown")
            report_timestamps[report_dir] = timestamp
            for run_key, data_list in run_key_to_metrics_map.items():
                if run_key not in report_input_args:
                    report_input_args[run_key] = data_list[0].get("run_metadata", {}).get("input_arguments", {})

    # Compare metrics across reports
    comparison_results = {}
    for report_dir, run_key_to_metrics_map in report_metrics.items():
        for run_key, metrics_list in run_key_to_metrics_map.items():
            if run_key not in comparison_results:
                comparison_results[run_key] = []
            # Collect all results that have the same input arguments
            comparison_results[run_key].extend(metrics_list)

    # Generate comparison data - only for run_keys that have multiple results to compare
    comparison_data = {}
    for run_key, results_list in comparison_results.items():
        if len(results_list) > 1:
            # Only compare if we have multiple results with the same input arguments
            comparison_data[run_key] = compare_metrics(results_list)
            print(f"Comparing {len(results_list)} runs with input arguments: {report_input_args.get(run_key, {})}")
        else:
            print(f"Skipping comparison for run key {run_key} - only {len(results_list)} result(s) found")

    # Generate charts
    chart_files = generate_charts(comparison_data, comparison_results, config)

    # Generate HTML report
    env = Environment()
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
            .chart-container { margin: 20px 0; text-align: center; }
            .chart-grid { max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 8px; }
            .run-section { border-bottom: 2px solid #ccc; padding-bottom: 30px; margin-bottom: 30px; }
        </style>
    </head>
    <body>
        <h1>Benchmark Comparison</h1>
        {% for run_key, metrics in comparison_data.items() %}
            {% set run_results = comparison_results[run_key] %}
            <div class="run-section">
                <h2>Run Configuration {{ loop.index }}</h2>
                <h3>Input Arguments</h3>
                <table>
                    <tr>
                        <th>Parameter</th>
                        <th>Value</th>
                    </tr>
                    {% for key, value in report_input_args[run_key].items() %}
                        <tr>
                            <td>{{ key }}</td>
                            <td>{{ value }}</td>
                        </tr>
                    {% endfor %}
                </table>
                
                {% if run_key in chart_files %}
                    <h3>Charts</h3>
                    <div class="chart-container">
                        {% for chart_file in chart_files[run_key] %}
                            <img src="{{ chart_file }}" alt="Benchmark Charts Grid" class="chart-grid">
                        {% endfor %}
                    </div>
                {% endif %}
                
                {% set categories = {} %}
                {% for path, values in metrics %}
                    {% set category = path[0] %}
                    {% if category not in categories %}
                        {% set _ = categories.update({category: []}) %}
                    {% endif %}
                    {% set _ = categories[category].append((path, values)) %}
                {% endfor %}
                {% for category, category_metrics in categories.items() %}
                    <h3>{{ category }} Metrics</h3>
                    <table>
                        <tr>
                            <th>Metric</th>
                            {% for result in run_results %}
                                <th>{{ result.get('run_metadata', {}).get('timestamp', 'Unknown') }}</th>
                            {% endfor %}
                        </tr>
                        {% for path, values in category_metrics %}
                            <tr>
                                <td>{{ '.'.join(path[1:]) }}</td>
                                {% for value in values %}
                                    <td>{{ "%.6f"|format(value) if value is number else value }}</td>
                                {% endfor %}
                            </tr>
                        {% endfor %}
                    </table>
                {% endfor %}
            </div>
        {% endfor %}
    </body>
    </html>
    """)
    html = template.render(comparison_data=comparison_data, comparison_results=comparison_results, chart_files=chart_files, report_input_args=report_input_args)
    with open(os.path.join(results_root, "benchmark_comparison_report.html"), "w") as f:
        f.write(html)

if __name__ == "__main__":
    main() 