#include "benchmark_tetrahedral_remeshing_common.h"
#include "mesh_quality.h"
#include <CGAL/Tetrahedral_remeshing/internal/elementary_remesh_impl.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <nlohmann/json.hpp>

#define Concurrency_tag CGAL::Sequential_tag
// Define C3t3 type wrapped around Remeshing_triangulation
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3<K>;
using Cb = CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3<K>;
using T3 = CGAL::Triangulation_3<K, CGAL::Triangulation_data_structure_3<Vb,Cb,Concurrency_tag>>;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<T3, int, int>;


int main(int argc, char** argv) {
  using namespace benchmarking;
  using nlohmann::json;
  json results_json;
  std::cout << std::setprecision(17);
  std::cerr << std::setprecision(17);
  if(argc != 7) {
    fatal_error(std::string("Usage: ") + argv[0] +
                " <input_mesh> <num_iterations> <remeshing_target_edge_factor> <smooth_constrained_edges> <num_threads> <results_json_path>");
  }
  std::string input = argv[1];
  int num_iterations = std::stoi(argv[2]);
  double remeshing_target_edge_factor = std::stod(argv[3]);
  bool smooth_constrained_edges = std::stoi(argv[4]) != 0;
  int num_threads = std::stoi(argv[5]);
  std::string results_json_path = argv[6];
  std::filesystem::create_directories(std::filesystem::path(std::filesystem::absolute(results_json_path)).parent_path());

  append_run_info(results_json, "Num_threads", 1);
  append_run_info(results_json, "Lockgrid_size", "N/A");
  append_run_info(results_json, "Lock_radius", "N/A");
  append_run_info(results_json, "Num_work_items_per_batch", "N/A");

  //Remeshing_triangulation tr;
  C3t3 c3t3;
  T3& tr = c3t3.triangulation();
  std::ifstream is(input, std::ios_base::in);
  if(!CGAL::IO::read_MEDIT(is, tr)){
    std::cerr << "Error: Could not read input mesh '" << input << "'" << std::endl;
    fatal_error(std::string("Could not read input mesh '") + input + "'");
  }

  std::cout << "Number of vertices: " << tr.number_of_vertices() << std::endl;
  std::cout << "Number of cells: " << tr.number_of_cells() << std::endl;

  // Extract input_name from input path
  std::string input_name = std::filesystem::path(input).stem().string();
  write_triangulation_info(results_json, c3t3.triangulation(), input_name);

  const double avg_edge_length = compute_average_edge_length(c3t3.triangulation());
  if(avg_edge_length <= 0.0) {
    fatal_error("Could not compute average edge length.");
  }
  double target_edge_length = avg_edge_length * remeshing_target_edge_factor;

  CGAL::Real_timer t;
  append_run_info(results_json, "Technique", "Atomic");
  append_run_info(results_json, "Edge Length", target_edge_length);

  CGAL::Real_timer t_atomic;
  t_atomic.start();

  CGAL::tetrahedral_isotropic_remeshing(c3t3, target_edge_length,
                                        CGAL::parameters::number_of_iterations(num_iterations).smooth_constrained_edges(smooth_constrained_edges));
  t_atomic.stop();

  append_metric_result(results_json, "Performance", "Total_Time", "Value", t_atomic.time());
  generate_quality_metrics(c3t3, results_json);

  append_metric_result(results_json, "Performance", "Memory", "Value", CGAL::Memory_sizer().virtual_size() >> 20);

  // Status
  append_execution_status(results_json, "success");

  // Write JSON
  write_results_json(results_json, results_json_path);
  return 0;
}
