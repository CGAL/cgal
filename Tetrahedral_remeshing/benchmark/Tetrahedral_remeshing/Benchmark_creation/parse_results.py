import os
import glob
import json
import pandas as pd
import argparse
import yaml

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

def load_schema(schema_path):
    with open(schema_path, 'r') as f:
        return yaml.safe_load(f)

def check_schema_recursive(schema, data, path=None):
    """
    Recursively check that all keys in schema exist at the correct level in data.
    schema: list or dict (from YAML)
    data: dict (from JSON)
    path: list of keys (for error messages)
    """
    if isinstance(schema, list):
        for item in schema:
            if isinstance(item, dict):
                for k, v in item.items():
                    if k not in data:
                        raise ValueError(f"{'/'.join(path or [])}: Missing key '{k}' in results.json")
                    check_schema_recursive(v, data[k], (path or []) + [k])
            else:
                if item not in data:
                    raise ValueError(f"{'/'.join(path or [])}: Missing key '{item}' in results.json")
    elif isinstance(schema, dict):
        for k, v in schema.items():
            if k not in data:
                raise ValueError(f"{'/'.join(path or [])}: Missing key '{k}' in results.json")
            check_schema_recursive(v, data[k], (path or []) + [k])
    else:
        # schema is a string (leaf key), already checked
        pass

def validate_results_json(data, schema, path=None):
    if not isinstance(data, dict):
        raise ValueError(f"{path or 'results.json'}: Not a dict at top level.")
    check_schema_recursive(schema, data)

def parse_results(results_root, output_csv='pipeline_results.csv', benchmark_name=None, schema_path='results_schema.yaml'):
    schema = load_schema(schema_path)
    rows = []
    for mesh_dir in glob.glob(os.path.join(results_root, '*')):
        if not os.path.isdir(mesh_dir):
            continue
        mesh_name = os.path.basename(mesh_dir)
        for results_file in glob.glob(os.path.join(mesh_dir, '*_results.json')):
            with open(results_file) as f:
                data = json.load(f)
            validate_results_json(data, schema, results_file)
            # Extract parameter set from filename (remove _results.json)
            params = os.path.basename(results_file).replace('_results.json', '')
            row = {'mesh': mesh_name, 'params': params}
            if benchmark_name:
                row['benchmark_name'] = benchmark_name
            # Only support new nested structure: flatten 'metrics' and 'run_metadata'
            if 'metrics' in data:
                row.update(flatten_dict(data['metrics'], parent_key='metrics'))
            if 'run_metadata' in data:
                row.update(flatten_dict(data['run_metadata'], parent_key='run_metadata'))
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False, sep=';')
    # Also save as JSON for robust comparison and provenance
    output_json = os.path.splitext(output_csv)[0] + '.json'
    df.to_dict(orient='records')
    with open(output_json, 'w', encoding='utf-8') as f:
        json.dump(df.to_dict(orient='records'), f, indent=2)
    print(f"Wrote benchmark results to {output_csv} and {output_json}")
    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Aggregate remeshing results into a benchmark results CSV.')
    parser.add_argument('results_root', type=str, help='Root directory containing mesh subfolders with *_results.json files')
    parser.add_argument('--output', type=str, default='pipeline_results.csv', help='Output CSV file')
    parser.add_argument('--benchmark_name', type=str, default=None, help='Name of the benchmark (for distinguishing in multi-benchmark runs)')
    parser.add_argument('--schema', type=str, default='results_schema.yaml', help='Path to the schema YAML file')
    args = parser.parse_args()
    parse_results(args.results_root, args.output, args.benchmark_name, schema_path=args.schema) 