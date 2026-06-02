# Generic Benchmarking Pipeline

This repository provides a **generic, schema-driven benchmarking pipeline** for scientific and engineering workflows. It is designed to be easily adapted to any domain by simply editing a YAML schema file and configuration.

---

## Features

- **Schema-driven validation**: Define your required results structure in a simple YAML file.
- **Configurable parameter sweeps**: Run benchmarks over any parameter grid.
- **Automatic results aggregation**: Collects and flattens results into a CSV and a full-fidelity JSON.
- **Interactive visualization**: Generate dynamic HTML reports with embedded charts and data tables.
- **Cross-run comparison**: Compare results across multiple benchmark runs with automatic chart generation.
- **Extensible and domain-agnostic**: No hardcoded field names or domain logic.
- **Rich metadata tracking**: Automatically captures environment info, git commit, and run configuration.

---

## Quickstart

1. **Prepare your benchmark executables**  
   Your executables should output a JSON file with a structure matching your schema (see below).

2. **Write your config file**  
   Example: `config.yaml`
   ```yaml
   pipeline:
     results_dir: ./Results
     run_benchmarks: true
     charts: true
     report: true
     report_chart_grid_cols: 2  # Configure chart grid layout

   benchmarks:
     - name: my_benchmark
       exec: ./my_benchmark_exec
       input_dir: ./input_data
       macros_config_file: ./my_benchmark_macros.h  # Optional
       sweep:
         param1: [1, 2]
         param2: [A, B]

   charts:
     - y: "metrics.Performance.Total_Time.Value"
       x: "run_metadata.run_info.Mesh.Cells"  # Optional, defaults to timestamp
       title: "Performance vs Mesh Size"
       style: "o-"  # Optional plotting style
   ```

3. **Write your results schema**  
   Example: `results_schema.yaml`
   ```yaml
   - metrics:
       - Performance
       - Quality
   - run_metadata:
       - status
       - run_info
       - input_arguments
   ```

4. **Run the pipeline**
   ```sh
   python run_benchmarks.py --config config.yaml
   ```

5. **Compare results across runs**
   ```sh
   python Benchmark_comparison/compare_benchmark_results.py
   ```

---

## Directory Structure

```
.
├── config.yaml
├── results_schema.yaml
├── run_benchmarks.py
├── parse_results.py
├── generate_charts.py
├── generate_report.py
├── generate_run_matrix.py
├── Benchmark_comparison/
│   └── compare_benchmark_results.py
├── Results/
│   ├── <timestamped_run>/
│   │   ├── my_benchmark/
│   │   │   ├── <input_name>/
│   │   │   │   ├── ..._results.json
│   │   │   ├── benchmark_results.csv
│   │   ├── pipeline_results.csv
│   │   ├── pipeline_results.json
│   │   ├── environment_metadata.json
│   │   ├── Charts/
│   │   │   ├── ...
│   │   ├── report.html
│   │   └── benchmark_comparison_report.html
```

---

## Main Scripts

- **run_benchmarks.py**: Orchestrates the pipeline, runs benchmarks, patches results, triggers aggregation, charting, and reporting.
- **parse_results.py**: Flattens and validates results JSONs into a CSV and a full-fidelity JSON, using your schema.
- **generate_charts.py**: Creates charts from the aggregated CSV, driven by your config.
- **generate_report.py**: Produces a self-contained HTML report with embedded charts, parameter sweeps, and raw results.
- **generate_run_matrix.py**: Generates the parameter sweep grid.
- **compare_benchmark_results.py**: Compares results across multiple runs, generating interactive HTML reports with embedded charts and data tables.

---

## Results Schema

- The schema is a YAML file describing the required structure of your results JSON.
- Each dash/indentation level corresponds to a level in the JSON.
- Only the presence of keys is checked; extra keys are ignored.

Example:
```yaml
- metrics:
    - Performance
    - Quality
- run_metadata:
    - status
    - run_info
    - input_arguments
```

---

## Preprocessor Macros Configuration

The pipeline supports tracking preprocessor macros that are defined when building your benchmarks. This allows you to:

- **Tag results** with the compile-time configuration used
- **Compare runs** with different macro configurations
- **Group results** by macro combinations for analysis

### Configuration

Add a `macros_config_file` field to your benchmark configuration:

```yaml
benchmarks:
  - name: my_benchmark
    exec: ./my_benchmark_exec
    macros_config_file: ./my_benchmark_macros.h  # Path to header file with macros
    sweep:
      param1: [1, 2]
```

### Macros Config File Format

The macros config file should be a simple C/C++ header file containing `#define` statements:

```c
#define FEATURE_A_ENABLED
#define FEATURE_B_ENABLED
// #define FEATURE_C_ENABLED  // Commented out = not defined
#define OPTIMIZATION_LEVEL 2
```

### Important Notes

- **Assumption**: The pipeline assumes that macros defined in the config file are actually used in your executable. **We do not validate that the macros are actually used in the code**.
- **Purpose**: Macros are parsed and used to **tag results** for comparison and grouping purposes only.
- **Parsing**: Only active `#define` statements are parsed; commented out defines (e.g., `// #define MACRO`) are ignored.
- **Comparison**: Runs are only compared if they have the same input arguments, exec metadata, AND the same set of defined macros.

### Results Structure

Parsed macros appear in the results JSON as:

```json
{
  "run_metadata": {
    "preprocessor_macros": [
      "FEATURE_A_ENABLED",
      "FEATURE_B_ENABLED",
      "OPTIMIZATION_LEVEL"
    ]
  }
}
```

---

## Data Aggregation and Comparison

- After each pipeline run, both a CSV (`pipeline_results.csv`) and a JSON (`pipeline_results.json`) are exported.
- **pipeline_results.json** is the canonical, full-fidelity data source for all downstream analysis and comparison.
- The CSV is provided for spreadsheet compatibility and quick tabular inspection.
- Environment metadata is captured in `environment_metadata.json`.
- Use `compare_benchmark_results.py` to generate interactive comparison reports across multiple runs.

---

## Extending and Customizing

- To add or remove required fields, just edit `results_schema.yaml`.
- To change parameter sweeps, edit the `sweep` section in your config.
- To add new charts or change chart types, edit the `charts` section in your config.
- To adjust chart grid layout, set `report_chart_grid_cols` in the pipeline config.

---

## Troubleshooting

- If a run fails schema validation, you'll get a clear error message showing the missing key and its path.
- You can specify a different schema file with `--schema` if needed.
- The pipeline automatically captures environment metadata to help reproduce and debug issues.

---

## More Details

See [README_DEEPDIVE.md](README_DEEPDIVE.md) for a detailed explanation of each script, the schema mechanism, and advanced usage. 