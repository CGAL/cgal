from itertools import product
import yaml
import argparse
import os

"""
Usage:
    python generate_run_matrix.py --config config.yaml
If --config is not provided, defaults to 'config.yaml'.
"""

def generate_run_matrix(sweep):
    """
    Generate the run matrix (list of parameter dicts) from a sweep dictionary.
    Each key in sweep maps to a list of values; the result is the Cartesian product.
    """
    keys = list(sweep.keys())
    values = [sweep[k] for k in keys]
    return [dict(zip(keys, v)) for v in product(*values)]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate run matrix from sweep in config file.")
    parser.add_argument('--config', type=str, default='config.yaml', help='Path to the YAML config file')
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    sweep = config['sweep']

    runs = generate_run_matrix(sweep)

    config_dir = os.path.dirname(os.path.abspath(args.config))
    out_path = os.path.join(config_dir, 'run_matrix.yaml')

    with open(out_path, 'w') as f:
        yaml.dump({'runs': runs}, f) 