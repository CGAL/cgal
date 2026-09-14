// Standalone experiment: read the already-glued/stitched mesh (exported by the
// domain-decomp extract-remesh-merge benchmark as glued_before_global.mesh)
// and time a SEQUENTIAL (Sequential_tag, no TBB dispatch, no lock grid) single
// global pass over it, for comparison against the parallel global pass
// (Global_Pass_Time in out_bear_8t.json), which VTune showed spends ~31% of
// its CPU time in lock/yield-spin contention.
//
// target_edge_length is passed explicitly (not derived from the input mesh's
// own average edge length) because the input here is already the fine glued
// mesh, not the original coarse bear.mesh.

#define USE_REFACTORED_TETRAHEDRAL_REMESHING
#include "benchmark_refactored_tetrahedral_remeshing_macros_config.h"
#include "benchmark_tetrahedral_remeshing_common.h"
#include "mesh_quality.h"

#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>
#include <CGAL/Real_timer.h>

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

// No CGAL_CONCURRENT_TETRAHEDRAL_REMESHING define -> Sequential_tag.
using K    = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb   = CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3<K>;
using Cb   = CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3<K>;
using T3   = CGAL::Triangulation_3<K, CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Sequential_tag>>;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<T3, int, int>;

int main(int argc, char** argv)
{
  using namespace benchmarking;
  using nlohmann::json;
  json results_json;
  std::cout << std::setprecision(17);

  if (argc < 5 || argc > 6)
  {
    fatal_error(std::string("Usage: ") + argv[0] +
                " <input_mesh> <target_edge_length> <num_iterations>"
                " <results_json_path> [remesh_boundaries: 0|1 (default 1)]");
  }
  const std::string input              = argv[1];
  const double      target_edge_length = std::stod(argv[2]);
  const int         num_iterations     = std::stoi(argv[3]);
  const std::string results_json_path  = argv[4];
  // The domain-decomposition run freezes the whole complex boundary in every
  // bucket, so it never remeshes the model surface. A baseline that DOES
  // remesh it is doing strictly more work, and neither the time nor the
  // quality numbers are then comparable. Pass 0 to freeze it here too.
  const bool        remesh_boundaries  = (argc >= 6) ? (std::stoi(argv[5]) != 0) : true;

  std::filesystem::create_directories(
    std::filesystem::path(std::filesystem::absolute(results_json_path)).parent_path());

  append_run_info(results_json, "Technique", "DomainDecomp-SeqGlobalPass");
  append_run_info(results_json, "Num_threads", 1);
  append_run_info(results_json, "Edge Length", target_edge_length);
  append_run_info(results_json, "Remesh_Boundaries", remesh_boundaries ? 1 : 0);
  std::cout << "[baseline] iterations=" << num_iterations
            << "  remesh_boundaries=" << (remesh_boundaries ? "true" : "false")
            << "  target_edge_length=" << target_edge_length << "\n";

  C3t3 c3t3;
  T3& tr = c3t3.triangulation();
  {
    std::ifstream is(input, std::ios_base::in);
    if (!CGAL::IO::read_MEDIT(is, tr))
      fatal_error(std::string("Could not read input mesh '") + input + "'");
  }
  std::cout << "Number of vertices: " << tr.number_of_vertices() << "\n";
  std::cout << "Number of cells: "    << tr.number_of_cells()    << "\n";
  write_triangulation_info(results_json, tr, std::filesystem::path(input).stem().string());

  CGAL::Real_timer t;
  t.start();
  CGAL::tetrahedral_isotropic_remeshing(
    c3t3, target_edge_length,
    CGAL::parameters::number_of_iterations(num_iterations)
                     .remesh_boundaries(remesh_boundaries)
                     .smooth_constrained_edges(true));
  t.stop();
  std::cout << "[sequential global pass] " << t.time() << "s -> "
            << tr.number_of_finite_cells() << " cells, "
            << tr.number_of_vertices() << " vertices\n";

  append_metric_result(results_json, "Performance", "Total_Time", "Value", t.time());
  generate_quality_metrics(c3t3, results_json);
  append_execution_status(results_json, "success");
  write_results_json(results_json, results_json_path);
  return 0;
}
