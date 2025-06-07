#include "benchmark_tetrahedral_remeshing_common.h"
#include "mesh_quality.h"
#include <CGAL/Adaptive_remeshing_sizing_field.h>
#include <CGAL/Real_timer.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <nlohmann/json.hpp>

int main(int argc, char** argv) {
    using namespace benchmarking;
    using nlohmann::json;
    json results_json;
    std::cout << std::setprecision(17);
    std::cerr << std::setprecision(17);
    if (argc != 5) {
        fatal_error(std::string("Usage: ") + argv[0] + " <input_mesh> <num_iterations> <num_threads> <results_json_path>");
    }
    std::string input = argv[1];
    int num_iterations = std::stoi(argv[2]);
    int num_threads = std::stoi(argv[3]);
    std::string results_json_path = argv[4];
    std::filesystem::create_directories(std::filesystem::path(results_json_path).parent_path());
#ifdef CGAL_CONCURRENT_TETRAHEDRAL_REMESHING
    Concurrent_tetrahedral_remesher_config::load_config_file(CONCURRENT_TETRAHEDRAL_REMESHING_CONFIG_FILENAME, true);
    if(num_threads <= 0) {
        fatal_error("num_threads must be positive.");
    }
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, num_threads);
    display_info(num_threads, "Task-scheduler (auto) (adaptive)");
    append_run_info(results_json, "Num_threads", num_threads);
    append_run_info(results_json, "Lockgrid_size", Concurrent_tetrahedral_remesher_config::get().locking_grid_num_cells_per_axis);
    append_run_info(results_json, "Lock_radius", Concurrent_tetrahedral_remesher_config::get().first_grid_lock_radius);
    append_run_info(results_json, "Num_work_items_per_batch", Concurrent_tetrahedral_remesher_config::get().num_work_items_per_batch);
#else
    std::cout << "-- Sequential Adaptive Tetrahedral Remeshing --" << std::endl;
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
        fatal_error(std::string("Could not read input mesh '") + input + "'");
    }
    // Extract input_name from input path
    std::string input_name = std::filesystem::path(input).stem().string();
    write_triangulation_info(results_json, tr, input_name);
    CGAL::Real_timer t;
    t.start();
    append_run_info(results_json, "Technique", "Adaptive");
    Remeshing_triangulation temp_tr = tr;
    CGAL::tetrahedral_isotropic_remeshing(temp_tr,
        CGAL::create_adaptive_remeshing_sizing_field(temp_tr),
        CGAL::parameters::number_of_iterations(num_iterations));
    t.stop();
    generate_quality_metrics(temp_tr, results_json);
    append_metric_result(results_json, "Performance", "Memory", "Value", CGAL::Memory_sizer().virtual_size() >> 20);
    append_metric_result(results_json, "Performance", "Total_Time", "Value", t.time());
    // Status
    append_execution_status(results_json, "success");
    // Write JSON
    write_results_json(results_json, results_json_path);
    return 0;
} 