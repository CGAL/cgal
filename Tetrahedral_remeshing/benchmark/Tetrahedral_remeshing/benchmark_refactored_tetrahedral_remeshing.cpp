#define USE_REFACTORED_TETRAHEDRAL_REMESHING
#define CGAL_CONCURRENT_TETRAHEDRAL_REMESHING
//#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE

#include "benchmark_refactored_tetrahedral_remeshing_macros_config.h"
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

 #ifdef CGAL_CONCURRENT_TETRAHEDRAL_REMESHING
   #define USE_THREADSAFE_INCIDENT_CELLS
 #endif

#ifdef CGAL_CONCURRENT_TETRAHEDRAL_REMESHING
  #include <tbb/version.h>
  #if TBB_VERSION_MAJOR >= 2018
    #include <tbb/global_control.h>
  #else
    #include <tbb/task_scheduler_init.h>
  #endif
#endif

#ifdef CGAL_CONCURRENT_TETRAHEDRAL_REMESHING
  #define Concurrency_tag CGAL::Parallel_tag
#else
  #define Concurrency_tag CGAL::Sequential_tag
#endif

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
  //std::filesystem::create_directories(std::filesystem::path(results_json_path).parent_path());
  std::filesystem::create_directories(std::filesystem::path(std::filesystem::absolute(results_json_path)).parent_path());

#ifdef CGAL_CONCURRENT_TETRAHEDRAL_REMESHING
  //Concurrent_tetrahedral_remesher_config::load_config_file(CONCURRENT_TETRAHEDRAL_REMESHING_CONFIG_FILENAME, true);
  if(num_threads <= 0) {
    fatal_error("num_threads must be positive.");
  }

  // Use TBB task scheduler with thread control
  // Note: tbb::global_control is available in TBB 2018+, fallback for older versions
  #if TBB_VERSION_MAJOR >= 2018
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, num_threads);
  #else
    tbb::task_scheduler_init scheduler_init(num_threads);
  #endif

  display_info(num_threads, "Task-scheduler (auto) (atomic)");
  append_run_info(results_json, "Num_threads", num_threads);
  // append_run_info(results_json, "Lockgrid_size", Concurrent_tetrahedral_remesher_config::get().locking_grid_num_cells_per_axis);
  // append_run_info(results_json, "Lock_radius", Concurrent_tetrahedral_remesher_config::get().first_grid_lock_radius);
  // append_run_info(results_json, "Num_work_items_per_batch", Concurrent_tetrahedral_remesher_config::get().num_work_items_per_batch);
#else
  if(num_threads != 1) {
    std::cerr << "Warning: num_threads argument ignored in sequential mode (should be 1)." << std::endl;
  }
  append_run_info(results_json, "Num_threads", 1);
  append_run_info(results_json, "Lockgrid_size", "N/A");
  append_run_info(results_json, "Lock_radius", "N/A");
  append_run_info(results_json, "Num_work_items_per_batch", "N/A");
#endif

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

  //std::cout << "Before remeshing:" << std::endl;
  //for(auto it=c3t3.finite_vertices_begin(); it!=c3t3.finite_vertices_end(); ++it) {
  //  if(it->in_dimension() == -1) {
  //    std::cout << "Vertex " << it->point() << " is not in the complex." << std::endl;
  //  }
  //}
  CGAL::Real_timer t_atomic;
  t_atomic.start();

  // Create C3t3 object wrapped around the remeshed triangulation
  //C3t3 c3t3;
  //c3t3.triangulation() = tr;
  // Create and run atomic remesher
  CGAL::tetrahedral_isotropic_remeshing(c3t3, target_edge_length,
                                        CGAL::parameters::number_of_iterations(num_iterations).smooth_constrained_edges(smooth_constrained_edges));
  t_atomic.stop();
  std::cout << "Remeshing took " << t_atomic.time() << std::endl;

  //std::cout << "After remeshing:" << std::endl;
  //for(auto it=tr.finite_vertices_begin(); it!=tr.finite_vertices_end(); ++it) {
  //  if(it->in_dimension() == -1) {
  //    std::cout << "Vertex " << it->point() << " is not in the complex." << std::endl;
  //  }
  //}
  //std::ofstream os("remeshed_allcells_"+input_name+".mesh");
  //CGAL::IO::write_MEDIT(os, tr, CGAL::parameters::all_cells(true));
  //CGAL::Tetrahedral_remeshing::debug::dump_triangulation_cells(tr, "triangulation_cells.mesh");
  //CGAL::Tetrahedral_remeshing::debug::dump_binary(tr, "binary_remeshed");
  append_metric_result(results_json, "Performance", "Total_Time", "Value", t_atomic.time());
  generate_quality_metrics(c3t3, results_json);

  append_metric_result(results_json, "Performance", "Memory", "Value", CGAL::Memory_sizer().virtual_size() >> 20);

  // Status
  append_execution_status(results_json, "success");

  // Write JSON
  write_results_json(results_json, results_json_path);
  //system("pause");
  return 0;
}
