import os
import json
import argparse
from datetime import datetime

def parse_file(filepath, num_lines):
    try:
        with open(filepath, 'r', encoding='utf-8') as file:
            return [file.readline().strip() for _ in range(num_lines)]
    except FileNotFoundError:
        print(f"Warning: File not found - {filepath}")
        return ["N/A"] * num_lines
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
        return ["ERROR"] * num_lines

def get_performance(file_name, benchmark_results_dir):
    filepath = os.path.join(benchmark_results_dir, "Performance", "results", f"{file_name}.log")
    seconds, memory_peaks = parse_file(filepath, 2)
    return {
        "seconds": seconds,
        "memory_peaks": memory_peaks
    }

def get_quality(file_name, benchmark_results_dir):
    filepath = os.path.join(benchmark_results_dir, "Quality", "results", f"{file_name}.log")
    (
        mean_min_angle,
        mean_max_angle,
        mean_radius_ratio,
        mean_edge_ratio,
        mean_aspect_ratio,
        complexity,
        almost_degenerate_triangles,
        hausdorff_distance
    ) = parse_file(filepath, 8)
    return {
        "Mean_Min_Angle_(degree)": mean_min_angle,
        "Mean_Max_Angle_(degree)": mean_max_angle,
        "Mean_Radius_Ratio": mean_radius_ratio,
        "Mean_Edge_Ratio": mean_edge_ratio,
        "Mean_Aspect_Ratio": mean_aspect_ratio,
        "Complexity_(#_of_triangle)": complexity,
        "#_of_almost_degenerate_triangle": almost_degenerate_triangles,
        "Hausdorff_distance_output_to_input_(%_of_bbox_diag)": hausdorff_distance
    }

def get_robustness(file_name, benchmark_results_dir):
    filepath = os.path.join(benchmark_results_dir, "Robustness", "results", f"{file_name}.log")
    robustness_flag = parse_file(filepath, 1)[0]

    tag_mapping = {
        "VALID_SOLID_OUTPUT": "OK",
        "INPUT_IS_INVALID": "Error",
        "OUTPUT_IS_NOT_TRIANGLE_MESH": "Error",
        "OUTPUT_IS_COMBINATORIAL_NON_MANIFOLD": "Error",
        "OUTPUT_HAS_BORDERS": "Error",
        "OUTPUT_HAS_DEGENERATED_FACES": "Error",
        "OUTPUT_HAS_GEOMETRIC_SELF_INTERSECTIONS": "Error",
        "OUTPUT_DOES_NOT_BOUND_VOLUME": "Error",
        "OUTPUT_DOES_NOT_CONTAIN_INPUT": "Error",
        "OUTPUT_DISTANCE_IS_TOO_LARGE": "Error",
        "TIMEOUT": "Timeout"
    }

    tag_name = robustness_flag if robustness_flag in tag_mapping else "UNKNOWN"
    tag = tag_mapping.get(robustness_flag, "UNKNOWN")

    return {
        "TAG_NAME": tag_name,
        "TAG": tag
    }

def process_all_files(meshes_dir, benchmark_results_dir, summary_json_path):
    valid_extensions = {'.off', '.obj', '.ply', '.stl', '.STL', '.ts', '.vtp'}
    results = {}
    dataset_file_counts = {}

    for root, _, files in os.walk(meshes_dir):
        for file in files:
            ext = os.path.splitext(file)[1].lower()
            if ext in valid_extensions:
                file_path = os.path.join(root, file)
                file_name = os.path.splitext(file)[0]

                rel_path = os.path.relpath(file_path, meshes_dir)
                sub_dirs = os.path.dirname(rel_path).split(os.sep)
                top_dir = sub_dirs[0] if sub_dirs else ""

                if top_dir not in dataset_file_counts:
                    dataset_file_counts[top_dir] = 0
                dataset_file_counts[top_dir] += 1

                if top_dir not in results:
                    results[top_dir] = {}

                metrics = {
                    "path": file_path,
                    "Performance": get_performance(file_name, benchmark_results_dir),
                    "Quality": get_quality(file_name, benchmark_results_dir),
                    "Robustness": get_robustness(file_name, benchmark_results_dir)
                }

                results[top_dir][file_name] = metrics

    final_json = {
        "Alpha_wrap_3": {
            **results,
            "dataset_file_counts": dataset_file_counts,
            "finished_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
    }

    with open(summary_json_path, 'w', encoding='utf-8') as json_file:
        json.dump(final_json, json_file, indent=4)

def main():
    parser = argparse.ArgumentParser(description="Alpha_wrap_3 benchmark metrics parser")
    parser.add_argument('--meshes-dir', required=True, help='Path to input mesh folder')
    parser.add_argument('--benchmark-results-dir', required=True, help='Path to benchmark output folder')
    parser.add_argument('--summary-json-path', required=True, help='Path to output JSON file')
    args = parser.parse_args()

    current_date = datetime.now().strftime("%Y-%m-%d")
    summary_json_path = os.path.join(
        args.summary_json_path,
        f"Alpha_wrap_3_results_{current_date}.json"
    )

    process_all_files(args.meshes_dir, args.benchmark_results_dir, summary_json_path)
    print(f"Results written to {args.summary_json_path}")

if __name__ == "__main__":
    main()
