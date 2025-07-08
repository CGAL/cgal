#pragma once
/**
 * @file benchmark_tetrahedral_remeshing_common.h
 * @brief Common typedefs, helpers, and utilities for tetrahedral remeshing benchmarks.
 *
 * This header is shared by both uniform and adaptive remeshing benchmark executables.
 */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Surface_mesh/Surface_mesh.h>
//#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Real_timer.h>
//#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
//#include <CGAL/Polygon_mesh_processing/self_intersections.h>
//#include <CGAL/Mesh_triangulation_3.h>
//#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
//#include <CGAL/Mesh_polyhedron_3.h>
//#include <CGAL/Polyhedral_mesh_domain_3.h>
//#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>
#include <CGAL/Memory_sizer.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <optional>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>

namespace benchmarking {

/// CGAL kernel and mesh typedefs
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
#ifdef CGAL_CONCURRENT_TETRAHEDRAL_REMESHING
typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K, CGAL::Parallel_tag> Remeshing_triangulation;
#else
typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;
#endif

/**
 * @brief Appends a metric result to a JSON object, grouped by category and metric label, under 'metrics'.
 * @param results_json The JSON object to append to.
 * @param category The metric category (e.g., "Performance").
 * @param metric_name The metric label (e.g., "Total_Time").
 * @param stat_name The stat (e.g., "Value").
 * @param value The value to store.
 */
inline void append_metric_result(nlohmann::json& results_json, const std::string& category, const std::string& metric_name, const std::string& stat_name, const nlohmann::json& value) {
    results_json["metrics"][category][metric_name][stat_name] = value;
}

/**
 * @brief Appends a run info entry to run_metadata.run_info.
 * @param results_json The JSON object to append to.
 * @param info_name The name of the info (e.g., "Lock_radius").
 * @param value The value to store (can be a dict or value).
 */
inline void append_run_info(nlohmann::json& results_json, const std::string& info_name, const nlohmann::json& value) {
    results_json["run_metadata"]["run_info"][info_name] = value;
}

/**
 * @brief Sets the execution status in run_metadata.status.
 * @param results_json The JSON object to update.
 * @param status The status string (e.g., "success").
 */
inline void append_execution_status(nlohmann::json& results_json, const std::string& status) {
    results_json["run_metadata"]["status"] = status;
}

/**
 * @brief Print an error message and exit the program.
 * @param msg The error message to print.
 * @param code The exit code (default 1).
 */
inline void fatal_error(const std::string& msg, int code = 1) {
    std::cerr << "[ERROR] " << msg << std::endl;
    std::exit(code);
}

/**
 * @brief Display the remeshing technique info.
 * @param num_threads Number of threads used.
 * @param technique String describing the technique.
 */
inline void display_info(int num_threads, const std::string& technique) {
    std::cout << technique << " (threads: " << num_threads << ")" << std::endl;
}

/**
 * @brief Compute the average edge length of a remeshing triangulation.
 * @tparam T Remeshing triangulation type.
 * @param tr The triangulation.
 * @return The average edge length.
 */
template<typename T>
double compute_average_edge_length(const T& tr) {
    double total_edge_length = 0.0;
    std::size_t num_edges = 0;
    for(typename T::Finite_cells_iterator cell_it = tr.finite_cells_begin(); cell_it != tr.finite_cells_end(); ++cell_it) {
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                const auto& p1 = cell_it->vertex(i)->point();
                const auto& p2 = cell_it->vertex(j)->point();
                total_edge_length += std::sqrt(CGAL::squared_distance(p1, p2));
                num_edges++;
            }
        }
    }
    if (num_edges == 0) return 0.0;
    return total_edge_length / num_edges;
}

/**
 * @brief Extract mesh domain info and append to results_json.
 * @param results_json The JSON object to append to.
 * @param tr The remeshing triangulation.
 * @param input_name The name of the input mesh.
 */
inline void write_triangulation_info(nlohmann::json& results_json, const Remeshing_triangulation& tr, const std::string& input_name) {
  append_run_info(results_json, "Mesh",
                  {{"Name", input_name},
                   {"Vertices", tr.number_of_vertices()},
                   {"Facets", tr.number_of_facets()},
                   {"Cells", tr.number_of_cells()}});
}

/**
 * @brief Write selected results to a JSON file.
 * @param results_json The JSON object to write.
 * @param results_json_path Path to the output JSON file.
 */
inline void write_results_json(const nlohmann::json& results_json, const std::string& results_json_path) {

    std::ofstream results_out(results_json_path);
    if (results_out) {
        results_out << results_json.dump(4);
        results_out.close();
        std::cout << "Wrote results to " << results_json_path << std::endl;
    } else {
        fatal_error("Failed to open " + results_json_path + " for writing!");
    }
}

} // namespace benchmarking 