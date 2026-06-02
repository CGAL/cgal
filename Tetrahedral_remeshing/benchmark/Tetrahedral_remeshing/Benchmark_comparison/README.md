# Benchmark Comparison Tool

This tool compares benchmark results across multiple benchmark runs, providing statistical analysis and visualization of performance metrics, quality measurements, and resource usage data.

## Features

- **Flexible Comparison Modes**: Compare within same executable or across different implementations
- **Statistical Analysis**: Automatic outlier detection and variance analysis
- **Interactive Charts**: Plotly-based charts with hover information and filtering
- **Comprehensive Reports**: HTML reports with detailed metric comparisons
- **Resource Usage Tracking**: Memory and CPU usage analysis (when available)

## Configuration Options

### Comparison Modes

The tool supports two distinct comparison scenarios controlled by the `group_by_executable` setting:

#### 1. Cross-Executable Comparison (`group_by_executable: false`)
**Use Case**: Compare different implementations that perform the same operation
- Compares `benchmark_tetrahedral_remeshing.exe` vs `benchmark_refactored_tetrahedral_remeshing.exe`
- Groups results by input arguments only (ignoring executable metadata)
- Perfect for comparing current vs refactored implementations
- Shows performance differences between algorithm variants

**Example**: Compare the current CGAL implementation against your refactored elementary operations framework

#### 2. Same-Executable Comparison (`group_by_executable: true`)
**Use Case**: Fine-tune and analyze variations within a single implementation
- Compares only runs from the same executable
- Groups results by input arguments AND executable metadata
- Perfect for performance tuning and consistency analysis
- Shows performance variations across different runs of the same code

**Example**: Analyze performance consistency of the refactored implementation across multiple runs

### Basic Configuration

```yaml
results_root: /path/to/Benchmark_results  # Root directory containing benchmark runs
include_reports: []                       # Specific subdirs to include (empty = all)

analysis:
  significant_change_threshold: 0.0       # Show metrics with >X% change (0 = all)
  outlier_threshold: 2.0                  # Highlight outliers >X std deviations
  group_by_executable: false              # false = cross-executable, true = same-executable

report:
  chart_grid_cols: 2                      # Chart grid layout columns

charts:
  - kind: line
    y: metrics.Performance.Total_Time.Value
    x: run_metadata.input_arguments.target_edge_factor
    title: "Performance vs Edge Factor"
    style: "o-"
```

## Usage Examples

### Comparing Current vs Refactored Implementation

1. **Run both benchmarks with identical parameters**:
   ```yaml
   # In Benchmark_creation/config.yaml
   benchmarks:
     - name: current
       exec: benchmark_tetrahedral_remeshing.exe
       sweep:
         input_mesh: [elephant.mesh]
         num_iterations: [50]
         target_edge_factor: [0.25, 0.5, 1.0, 2.0]
         threads: [1]
     - name: refactored  
       exec: benchmark_refactored_tetrahedral_remeshing.exe
       sweep:
         input_mesh: [elephant.mesh]
         num_iterations: [50]
         target_edge_factor: [0.25, 0.5, 1.0, 2.0]
         threads: [1]
   ```

2. **Configure comparison for cross-executable analysis**:
   ```yaml
   # In Benchmark_comparison/config.yaml
   analysis:
     group_by_executable: false  # Enable cross-executable comparison
   
   charts:
     - kind: line
       y: metrics.Performance.Total_Time.Value
       x: run_metadata.input_arguments.target_edge_factor
       title: "Current vs Refactored: Performance by Edge Factor"
       style: "o-"
   ```

3. **Run comparison**:
   ```bash
   python compare_benchmark_results.py --config config.yaml
   ```

### Fine-Tuning Single Implementation

1. **Configure for same-executable analysis**:
   ```yaml
   analysis:
     group_by_executable: true   # Compare only within same executable
     outlier_threshold: 1.0      # Sensitive outlier detection for consistency analysis
   ```

2. **Run multiple times with same parameters** to analyze consistency
3. **Analyze variance and outliers** in the generated report

## Chart Types and Filtering

### Time-Series Charts
```yaml
- kind: line
  y: metrics.Performance.Total_Time.Value
  title: "Performance Over Time"
  # No x specified = uses timestamp
```

### Parameter-Based Charts
```yaml
- kind: line
  y: metrics.Performance.Total_Time.Value
x: run_metadata.input_arguments.threads
  title: "Performance vs Thread Count"
```

### Filtered Charts
```yaml
- kind: line
  y: metrics.Performance.Total_Time.Value
  x: run_metadata.input_arguments.target_edge_factor
  title: "Performance for Large Meshes Only"
    filter:
    run_metadata.run_info.Mesh.Name: elephant
```

## Understanding Results

### Metric Comparison Tables
- **Green cells**: Low variance (CV < 5%)
- **Yellow cells**: Medium variance (5% ≤ CV ≤ 20%)
- **Red cells**: Outliers (beyond threshold standard deviations)
- **Hover tooltips**: Detailed statistics (mean, std dev, percentiles)

### Interactive Charts
- **Series**: Automatically grouped by benchmark + input parameters
- **Hover information**: Shows benchmark name, input arguments, and metric values
- **Legend**: Distinguishes between different implementations or parameter sets
- **Responsive**: Adapts to different screen sizes

### Statistical Analysis
- **Coefficient of Variation (CV%)**: Relative variability measure
- **Outlier Detection**: Z-score based identification
- **Percentile Analysis**: Min, Q1, Median, Q3, Max values
- **Cross-Implementation Comparison**: When enabled, shows both implementations in same tables

## File Structure

```
Benchmark_results/
├── 2024-01-15_10-30-00/    # Timestamp-based run directories
│   └── pipeline_results.json
├── 2024-01-15_11-45-00/
│   └── pipeline_results.json
└── benchmark_comparison_report.html  # Generated report
```

## Troubleshooting

### No Comparison Data Available
**Problem**: Report shows "No comparison data available"
**Solutions**:
- Ensure both executables use identical input parameters
- Check that `group_by_executable` setting matches your comparison intent
- Verify that `pipeline_results.json` files contain both implementations

### Charts Not Showing Expected Comparisons
**Problem**: Charts don't show both implementations
**Solutions**:
- Set `group_by_executable: false` for cross-executable comparison
- Ensure both executables have runs with identical input arguments
- Check chart filters aren't excluding one implementation

### High Variance Warnings
**Problem**: Many metrics marked as high variance
**Solutions**:
- Run more iterations to reduce natural variance
- Check for system interference during benchmarking
- Consider if differences are due to algorithmic improvements vs measurement noise

## Advanced Features

### Custom Grouping
The tool automatically handles different executable metadata when `group_by_executable: false`, allowing meaningful comparison between:
- Different compilation flags
- Different algorithm implementations
- Different versions of the same code

### Resource Usage Integration
When resource monitoring is enabled, the tool automatically includes:
- Peak memory usage
- Average CPU utilization
- Memory allocation patterns
- All integrated into the same comparison framework

### Extensible Chart System
Charts support various customizations:
- Multiple metric types (performance, quality, resource usage)
- Flexible filtering by any result field
- Multiple chart types (line, scatter, etc.)
- Custom styling and layout options