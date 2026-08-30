from itertools import product
import os
import logging
import yaml
import argparse
import json
import hashlib
from pathlib import Path
import glob

def setup_logging():
    """Setup logging for run matrix generation."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)

def validate_path(path):
    """Validate that a path exists and is accessible."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Path does not exist: {path}")
    if not os.access(path, os.R_OK):
        raise PermissionError(f"Path is not readable: {path}")
    return True

def expand_sweep_value(val, base_dir=None):
    """
    Expand sweep values with improved error handling and more options.
    
    Args:
        val: The value to expand (can be string, list, dict, or other)
        base_dir: Base directory for relative paths
    
    Returns:
        List of expanded values
    """
    logger = logging.getLogger(__name__)
    
    # If it's a list, expand each item recursively and flatten
    if isinstance(val, list):
        expanded = []
        for item in val:
            try:
                expanded.extend(expand_sweep_value(item, base_dir))
            except Exception as e:
                logger.warning(f"Failed to expand list item {item}: {e}")
                continue
        return expanded
    
    # If it's a dict, handle special expansion patterns
    elif isinstance(val, dict):
        if 'range' in val:
            # Handle numeric ranges: {'range': [start, stop, step]}
            range_params = val['range']
            if len(range_params) == 2:
                start, stop = range_params
                step = 1
            elif len(range_params) == 3:
                start, stop, step = range_params
            else:
                raise ValueError(f"Range must have 2 or 3 parameters, got {len(range_params)}")
            
            if isinstance(start, float) or isinstance(stop, float) or isinstance(step, float):
                # Handle float ranges
                import numpy as np
                return np.arange(start, stop, step).tolist()
            else:
                # Handle integer ranges
                return list(range(start, stop, step))
        
        elif 'glob' in val:
            # Handle glob patterns: {'glob': 'pattern', 'base_dir': 'optional'}
            pattern = val['glob']
            pattern_base_dir = val.get('base_dir', base_dir or '.')
            full_pattern = os.path.join(pattern_base_dir, pattern)
            matches = glob.glob(full_pattern)
            if not matches:
                logger.warning(f"Glob pattern '{full_pattern}' matched no files")
                return []
            return sorted(matches)
        
        elif 'files' in val:
            # Handle explicit file lists: {'files': [...], 'base_dir': 'optional'}
            files = val['files']
            file_base_dir = val.get('base_dir', base_dir or '.')
            expanded = []
            for file in files:
                if os.path.isabs(file):
                    full_path = file
                else:
                    full_path = os.path.join(file_base_dir, file)
                
                if os.path.exists(full_path):
                    expanded.append(full_path)
                else:
                    logger.warning(f"File not found: {full_path}")
            return expanded
        
        else:
            # Unknown dict format, treat as single value
            logger.warning(f"Unknown dict format for sweep value: {val}")
            return [val]
    
    # If it's a string, check if it's a file or directory
    elif isinstance(val, str):
        # Handle relative paths
        if base_dir and not os.path.isabs(val):
            full_path = os.path.join(base_dir, val)
        else:
            full_path = val
        
        if os.path.isdir(full_path):
            try:
                validate_path(full_path)
                files = [os.path.join(full_path, f) for f in os.listdir(full_path) 
                        if os.path.isfile(os.path.join(full_path, f))]
                if not files:
                    raise FileNotFoundError(f"Sweep directory {full_path} is empty.")
                return sorted(files)  # Sort for reproducible results
            except Exception as e:
                logger.error(f"Failed to expand directory {full_path}: {e}")
                raise
                
        elif os.path.isfile(full_path):
            try:
                validate_path(full_path)
                return [full_path]
            except Exception as e:
                logger.error(f"Failed to access file {full_path}: {e}")
                raise
        else:
            # Not a file or directory, treat as a string value (e.g., a label)
            return [val]
    
    # If it's a number or other type, just wrap in a list
    else:
        return [val]

def validate_sweep_config(sweep):
    """Validate the sweep configuration."""
    if not isinstance(sweep, dict):
        raise ValueError("Sweep configuration must be a dictionary")
    
    if not sweep:
        raise ValueError("Sweep configuration cannot be empty")
    
    for key, value in sweep.items():
        if not isinstance(key, str):
            raise ValueError(f"Sweep key must be string, got {type(key)}")
        
        if key.strip() != key:
            raise ValueError(f"Sweep key '{key}' contains leading/trailing whitespace")
        
        if not key:
            raise ValueError("Sweep key cannot be empty")

def generate_run_matrix(sweep, base_dir=None, max_combinations=None):
    """
    Generate the run matrix (list of parameter dicts) from a sweep dictionary.
    Each key in sweep maps to a list of values; the result is the Cartesian product.
    Handles numeric, file, and directory sweep values.
    
    Args:
        sweep: Dictionary of parameter names to sweep values
        base_dir: Base directory for relative paths
        max_combinations: Maximum number of combinations to generate (safety limit)
    
    Returns:
        List of parameter dictionaries
    """
    logger = logging.getLogger(__name__)
    
    validate_sweep_config(sweep)
    
    logger.info(f"Generating run matrix for {len(sweep)} parameters")
    
    keys = list(sweep.keys())
    expanded_values = []
    
    for key in keys:
        try:
            values = expand_sweep_value(sweep[key], base_dir)
            if not values:
                raise ValueError(f"No values generated for parameter '{key}'")
            expanded_values.append(values)
            logger.info(f"Parameter '{key}': {len(values)} values")
        except Exception as e:
            logger.error(f"Failed to expand parameter '{key}': {e}")
            raise
    
    # Calculate total combinations
    total_combinations = 1
    for values in expanded_values:
        total_combinations *= len(values)
    
    logger.info(f"Total combinations: {total_combinations}")
    
    # Safety check
    if max_combinations and total_combinations > max_combinations:
        raise ValueError(f"Too many combinations ({total_combinations}), maximum allowed: {max_combinations}")
    
    # Generate the Cartesian product
    run_matrix = [dict(zip(keys, combination)) for combination in product(*expanded_values)]
    
    logger.info(f"Generated run matrix with {len(run_matrix)} runs")
    return run_matrix

def save_run_matrix(run_matrix, output_path, include_metadata=True):
    """Save run matrix to file with optional metadata."""
    logger = logging.getLogger(__name__)
    
    output_data = {'runs': run_matrix}
    
    if include_metadata:
        output_data['metadata'] = {
            'total_runs': len(run_matrix),
            'generated_at': __import__('datetime').datetime.now().isoformat(),
            'parameters': list(run_matrix[0].keys()) if run_matrix else []
        }
    
    try:
        with open(output_path, 'w') as f:
            yaml.dump(output_data, f, default_flow_style=False, sort_keys=False)
        logger.info(f"Saved run matrix to {output_path}")
    except Exception as e:
        logger.error(f"Failed to save run matrix: {e}")
        raise

def load_cached_matrix(cache_path, config_hash):
    """Load cached run matrix if it exists and is valid."""
    if not os.path.exists(cache_path):
        return None
    
    try:
        with open(cache_path, 'r') as f:
            cached_data = yaml.safe_load(f)
        
        if (cached_data.get('config_hash') == config_hash and 
            'runs' in cached_data):
            logging.info(f"Using cached run matrix from {cache_path}")
            return cached_data['runs']
    except Exception as e:
        logging.warning(f"Failed to load cached matrix: {e}")
    
    return None

def save_cached_matrix(run_matrix, cache_path, config_hash):
    """Save run matrix to cache with config hash."""
    try:
        cache_data = {
            'config_hash': config_hash,
            'runs': run_matrix,
            'cached_at': __import__('datetime').datetime.now().isoformat()
        }
        with open(cache_path, 'w') as f:
            yaml.dump(cache_data, f)
        logging.info(f"Cached run matrix to {cache_path}")
    except Exception as e:
        logging.warning(f"Failed to cache run matrix: {e}")

def main():
    """Main entry point when run as script."""
    parser = argparse.ArgumentParser(description="Generate run matrix from sweep in config file.")
    parser.add_argument('--config', type=str, default='config.yaml', 
                       help='Path to the YAML config file')
    parser.add_argument('--output', type=str, 
                       help='Output file path (default: run_matrix.yaml in config directory)')
    parser.add_argument('--max-combinations', type=int, default=10000,
                       help='Maximum number of combinations to generate (safety limit)')
    parser.add_argument('--no-cache', action='store_true',
                       help='Disable caching of run matrix')
    parser.add_argument('--base-dir', type=str,
                       help='Base directory for relative paths in sweep values')
    args = parser.parse_args()

    logger = setup_logging()
    
    try:
        # Load configuration
        with open(args.config) as f:
            config = yaml.safe_load(f)
        logger.info(f"Loaded configuration from {args.config}")
        
        if 'sweep' not in config:
            raise ValueError("Configuration file must contain 'sweep' section")
        
        sweep = config['sweep']

        # Determine base directory
        base_dir = args.base_dir or os.path.dirname(os.path.abspath(args.config))
        
        # Generate config hash for caching
        config_str = yaml.dump(sweep, sort_keys=True)
        config_hash = hashlib.md5(config_str.encode()).hexdigest()
        
        # Check cache if enabled
        run_matrix = None
        if not args.no_cache:
            cache_path = os.path.join(base_dir, '.run_matrix_cache.yaml')
            run_matrix = load_cached_matrix(cache_path, config_hash)
        
        # Generate run matrix if not cached
        if run_matrix is None:
            run_matrix = generate_run_matrix(sweep, base_dir, args.max_combinations)
            
            # Save to cache if enabled
            if not args.no_cache:
                save_cached_matrix(run_matrix, cache_path, config_hash)
        
        # Determine output path
        if args.output:
            output_path = args.output
        else:
            config_dir = os.path.dirname(os.path.abspath(args.config))
            output_path = os.path.join(config_dir, 'run_matrix.yaml')
        
        # Save run matrix
        save_run_matrix(run_matrix, output_path)
        
        logger.info(f"Successfully generated run matrix with {len(run_matrix)} runs")
        
    except Exception as e:
        logger.error(f"Run matrix generation failed: {e}")
        return 1
    
    return 0

if __name__ == '__main__':
    exit(main()) 