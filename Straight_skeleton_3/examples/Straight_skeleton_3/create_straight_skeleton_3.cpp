#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Straight_skeleton_3/create_straight_skeleton_3.h>

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
  std::cout.precision(17);
  std::cerr.precision(17);

  if (argc == 1) {
    std::cout << "Usage: " << argv[0] << " <input filename> [time bound]" << std::endl;
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

  const FT time_bound = (argc > 2) ? std::stod(argv[2]) : 1;

  // every face has weight 2.0, meaning faces move at twice the default speed
  CGAL::Constant_property_map<face_descriptor, FT> face_weight_map(2);

  CGAL::Real_timer timer;
  timer.start();

  // Main call
  auto straight_skeleton = CGAL::create_straight_skeleton_3(sm, time_bound,
                                                            CGAL::parameters::face_weight_map(face_weight_map));

  timer.stop();
  std::cout << "Elapsed: " << timer.time() << std::endl;

  if (!straight_skeleton) {
    std::cerr << "Error: failed to create straight skeleton." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Straight skeleton created:\n" << straight_skeleton->to_string() << std::endl;

  std::cout << "Done" << std::endl;
  return EXIT_SUCCESS;
}
