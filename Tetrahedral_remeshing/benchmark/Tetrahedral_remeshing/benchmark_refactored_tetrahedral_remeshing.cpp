//#define CGAL_REFACTORED_TETRAHEDRAL_REMESHING_DEBUG
#define CGAL_TETRAHEDRAL_REMESHING_USE_ELEMENTARY
#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
//#define USE_THREADSAFE_INCIDENT_CELLS
//#define ENABLE_EDGE_FLIP_DEBUG
#include "benchmark_tetrahedral_remeshing_common.h"
#include "mesh_quality.h"
#include <CGAL/Tetrahedral_remeshing/internal/elementary_remesh_impl.h>
#include <CGAL/Real_timer.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <nlohmann/json.hpp>

#ifdef CGAL_CONCURRENT_TETRAHEDRAL_REMESHING
  #include <tbb/version.h>
  #if TBB_VERSION_MAJOR >= 2018
    #include <tbb/global_control.h>
  #else
    #include <tbb/task_scheduler_init.h>
  #endif
#endif

int main(int argc, char** argv) {
  using namespace benchmarking;
  using nlohmann::json;
  json results_json;
  std::cout << std::setprecision(17);
  std::cerr << std::setprecision(17);
  if(argc != 6) {
    fatal_error(std::string("Usage: ") + argv[0] +
                " <input_mesh> <num_iterations> <remeshing_target_edge_factor> <num_threads> <results_json_path>");
  }
  std::string input = argv[1];
  int num_iterations = std::stoi(argv[2]);
  double remeshing_target_edge_factor = std::stod(argv[3]);
  int num_threads = std::stoi(argv[4]);
  std::string results_json_path = argv[5];
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
  std::cout << "-- Sequential Atomic Tetrahedral Remeshing --" << std::endl;
  if(num_threads != 1) {
    std::cerr << "Warning: num_threads argument ignored in sequential mode (should be 1)." << std::endl;
  }
  append_run_info(results_json, "Num_threads", 1);
  append_run_info(results_json, "Lockgrid_size", "N/A");
  append_run_info(results_json, "Lock_radius", "N/A");
  append_run_info(results_json, "Num_work_items_per_batch", "N/A");
#endif

  Remeshing_triangulation tr;
  std::ifstream is(input, std::ios_base::in);
  if(!CGAL::IO::read_MEDIT(is, tr)) {
    std::cerr << "Error: Could not read input mesh '" << input << "'" << std::endl;
    fatal_error(std::string("Could not read input mesh '") + input + "'");
  }

  // Extract input_name from input path
  std::string input_name = std::filesystem::path(input).stem().string();
  write_triangulation_info(results_json, tr, input_name);

  const double avg_edge_length = compute_average_edge_length(tr);
  if(avg_edge_length <= 0.0) {
    fatal_error("Could not compute average edge length.");
  }
  double target_edge_length = avg_edge_length * remeshing_target_edge_factor;

  CGAL::Real_timer t;
  append_run_info(results_json, "Technique", "Atomic");
  append_run_info(results_json, "Edge Length", target_edge_length);

  CGAL::Real_timer t_atomic;
  t_atomic.start();

  // Create and run atomic remesher
  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length, CGAL::parameters::number_of_iterations(num_iterations));

  t_atomic.stop();
  append_metric_result(results_json, "Performance", "Total_Time", "Value", t_atomic.time());
  generate_quality_metrics(tr, results_json);

  append_metric_result(results_json, "Performance", "Memory", "Value", CGAL::Memory_sizer().virtual_size() >> 20);

  // Status
  append_execution_status(results_json, "success");

  // Write JSON
  write_results_json(results_json, results_json_path);

  return 0;
}