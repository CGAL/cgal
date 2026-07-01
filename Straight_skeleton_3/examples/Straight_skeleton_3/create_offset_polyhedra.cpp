#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Straight_skeleton_3/create_offset_polyhedra.h>

#include <CGAL/Real_timer.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/property_map.h>

#include <iostream>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace SS3 = CGAL::Straight_skeletons_3;

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = K::FT;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point_3>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

int main(int argc, char** argv)
{
  if (argc == 1) {
    std::cout << "Usage: " << argv[0] << " <input filename>" << std::endl;
    std::cout << "  <input filename>: Path to input mesh file (required)" << std::endl;
    return EXIT_FAILURE;
  }

  const char* mesh_filename = argv[1];

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(mesh_filename, sm) || CGAL::is_empty(sm) || !is_valid_face_graph(sm)) {
    std::cerr << "Error: failed to read input " << mesh_filename << std::endl;
    return EXIT_FAILURE;
  }

  // let's tolerate untriangulated inputs
  if (!CGAL::is_triangle_mesh(sm)) {
    std::cout << "Triangulating input mesh..." << std::endl;
    PMP::triangulate_faces(sm);
  }

  std::cout << "Input mesh: " << vertices(sm).size() << " NV " << faces(sm).size() << " NF" << std::endl;

  std::vector<FT> save_times = { 1, 2, 3 };

  // every face has weight 3.0, meaning faces move at thrice the default speed
  CGAL::Constant_property_map<face_descriptor, FT> face_weight_map(3);

  std::vector<Mesh> results;
  results.reserve(save_times.size());

  CGAL::Real_timer timer;
  timer.start();

  // Main call
  bool success = CGAL::create_straight_skeleton_and_offset_polyhedra_3(sm, save_times, results, CGAL::parameters::face_weight_map(face_weight_map));

  timer.stop();
  std::cout << "Elapsed: " << timer.time() << std::endl;

  // save the results
  for (std::size_t i=0; i<results.size(); ++i) {
    const FT save_time = save_times[i];
    const Mesh& sm = results[i];

    std::stringstream out_ss;
    out_ss <<  "result_" << save_time << ".obj";

    if (!CGAL::IO::write_polygon_mesh(out_ss.str(), sm, CGAL::parameters::stream_precision(17))) {
      std::cerr << "Error: failed to write result " << out_ss.str() << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "Done" << std::endl;
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
