import os
import sys
import yaml
import pandas as pd
import datetime
from jinja2 import Environment
import base64
import json

def image_to_base64(path):
    with open(path, 'rb') as img_f:
        encoded = base64.b64encode(img_f.read()).decode('utf-8')
    ext = os.path.splitext(path)[1].lower()
    mime = 'image/png' if ext == '.png' else 'image/jpeg'
    return f"data:{mime};base64,{encoded}"

def run(config, aggregated_data, charts_dir='Charts', report_dir=None):
    # Prepare data for the report
    date = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')

    # Grouped by benchmark: parameter sweep and raw results
    report_metrics = config.get('report_metrics', None)
    benchmark_sections = []
    if 'benchmark_name' in aggregated_data.columns:
        benchmarks = aggregated_data['benchmark_name'].unique()
    else:
        benchmarks = ['all']
    for benchmark in benchmarks:
        bench_df = aggregated_data[aggregated_data['benchmark_name'] == benchmark] if benchmark != 'all' else aggregated_data
        # Get sweep parameter names from config for this benchmark
        sweep_params = []
        for bench_cfg in config.get('benchmarks', []):
            if bench_cfg.get('name', 'benchmark') == benchmark:
                sweep_params = list(bench_cfg.get('sweep', {}).keys())
                break
        # Find the corresponding columns in the DataFrame (look for run_metadata.input_arguments.<param>)
        param_cols = [f'run_metadata.input_arguments.{p}' for p in sweep_params if f'run_metadata.input_arguments.{p}' in bench_df.columns]
        if param_cols:
            param_table = bench_df[param_cols].drop_duplicates().to_html(index=False)
        else:
            param_table = '<i>No sweep parameters found in results for this benchmark.</i>'
        # Raw results table
        if report_metrics:
            entry = next((rm for rm in report_metrics if rm['benchmark'] == benchmark), None)
            if entry and 'columns' in entry:
                cols = [col for col in entry['columns'] if col in bench_df.columns]
                raw_table = bench_df[cols].to_html(index=False)
            else:
                raw_table = bench_df.to_html(index=False)
        else:
            raw_table = bench_df.to_html(index=False)
        benchmark_sections.append({'benchmark': benchmark, 'param_table': param_table, 'raw_table': raw_table})

    chart_paths = [os.path.join(charts_dir, f) for f in os.listdir(charts_dir) if f.endswith('.png') or f.endswith('.jpg')]
    chart_paths.sort()
    charts = [image_to_base64(path) for path in chart_paths]

    # --- Load metadata.json if present ---
    metadata = None
    if report_dir is None:
        report_dir = os.path.dirname(os.path.abspath(charts_dir))
    metadata_path = os.path.join(report_dir, 'environment_metadata.json')
    if os.path.isfile(metadata_path):
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)

    # Get chart grid columns from config
    report_chart_grid_cols = config.get('pipeline', {}).get('report_chart_grid_cols', 2)

    # Jinja2 template
    template_str = '''
    <html>
    <head><title>Remeshing Benchmark Report</title></head>
    <body>
    <h1>Remeshing Benchmark Report</h1>
    <p><b>Date:</b> {{ date }}</p>
    {% for section in benchmark_sections %}
      <h2>Benchmark: {{ section.benchmark }}</h2>
      <h3>Parameter Sweep</h3>
      {{ section.param_table | safe }}
      <h3>Raw Results</h3>
      {{ section.raw_table | safe }}
    {% endfor %}
    <h2>Charts</h2>
    <table>
      {% for row in charts|batch(report_chart_grid_cols, '') %}
        <tr>
          {% for chart in row %}
            <td>{% if chart %}<img src="{{ chart }}" width="600"><br>{% endif %}</td>
          {% endfor %}
        </tr>
      {% endfor %}
    </table>
    {% if metadata %}
    <h2>Run Metadata</h2>
    <table border="1" cellpadding="4" cellspacing="0">
      {% for k, v in metadata.items() %}
        <tr><td><b>{{ k }}</b></td><td><code>{{ v }}</code></td></tr>
      {% endfor %}
    </table>
    {% endif %}
    </body>
    </html>
    '''
    env = Environment()
    template = env.from_string(template_str)
    html = template.render(date=date, benchmark_sections=benchmark_sections, charts=charts, metadata=metadata, report_chart_grid_cols=report_chart_grid_cols)
    os.makedirs(report_dir, exist_ok=True)
    report_path = os.path.join(report_dir, 'report.html')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"Wrote HTML report to {report_path}")
    return report_path

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python generate_report.py <config.yaml> <results.csv>")
        sys.exit(1)
    config_path = sys.argv[1]
    csv_path = sys.argv[2]
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    # Auto-detect delimiter
    def read_csv_with_auto_delimiter(csv_path):
        return pd.read_csv(csv_path, delimiter=';')
    df = read_csv_with_auto_delimiter(csv_path)
    charts_dir = os.path.join(os.path.dirname(os.path.abspath(csv_path)), 'Charts')
    report_dir = os.path.dirname(os.path.abspath(csv_path))
    run(config, df, charts_dir=charts_dir, report_dir=report_dir) 