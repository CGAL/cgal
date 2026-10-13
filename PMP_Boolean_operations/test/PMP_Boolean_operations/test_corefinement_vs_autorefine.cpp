#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangle_soup_boolean_operations_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/IO/polygon_mesh_io.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace CFR = PMP::Corefinement;

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
using Polygon = std::vector<std::size_t>;
using Polygons = std::vector<Polygon>;
using Points = std::vector<Point_3>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

using Triangle = std::array<Point_3, 3>;

void test(const std::string& p_filename, const std::string& q_filename)
{
  Mesh P, Q;

  const std::filesystem::path p_path(p_filename);
  const std::filesystem::path q_path(q_filename);
  std::cout << "Testing " << p_path.filename().string() << " and " << q_path.filename().string() << '\n';

  if (!CGAL::IO::read_polygon_mesh(p_filename, P))
  {
    std::cerr << "Failed to read polygon mesh: " << p_filename << '\n';
    std::exit(EXIT_FAILURE);
  }
  if (!CGAL::IO::read_polygon_mesh(q_filename, Q))
  {
    std::cerr << "Failed to read polygon mesh: " << q_filename << '\n';
    std::exit(EXIT_FAILURE);
  }

  // Convert the inputs to triangle soups.
  Points p_points, q_points;
  Polygons p_polygons, q_polygons;
  PMP::polygon_mesh_to_polygon_soup(P, p_points, p_polygons);
  PMP::polygon_mesh_to_polygon_soup(Q, q_points, q_polygons);

  // Results computed directly on triangle soups.
  const std::array<const char*, 4> names = {{
    "union", "intersection", "P - Q", "Q - P"
  }};

  std::array<Points, 4> soup_points;
  std::array<Polygons, 4> soup_polygons;

  PMP::compute_union(p_points, p_polygons, q_points, q_polygons, soup_points[CFR::UNION], soup_polygons[CFR::UNION]);
  PMP::compute_intersection(p_points, p_polygons, q_points, q_polygons, soup_points[CFR::INTERSECTION], soup_polygons[CFR::INTERSECTION]);
  PMP::compute_difference(p_points, p_polygons, q_points, q_polygons, soup_points[CFR::TM1_MINUS_TM2], soup_polygons[CFR::TM1_MINUS_TM2]);
  PMP::compute_difference(q_points, q_polygons, p_points, p_polygons, soup_points[CFR::TM2_MINUS_TM1], soup_polygons[CFR::TM2_MINUS_TM1]);

  // Results computed by the existing mesh API.
  Mesh union_mesh, intersection_mesh, p_minus_q_mesh, q_minus_p_mesh;
  std::array<std::optional<Mesh*>, 4> corefine_outputs{&union_mesh, &intersection_mesh, &p_minus_q_mesh, &q_minus_p_mesh};
  const std::array<bool, 4> corefine_success = PMP::corefine_and_compute_boolean_operations(P, Q, corefine_outputs);

  for (std::size_t i = 0; i < 4; ++i)
  {
    if (!corefine_success[i] ){
      std::cout << names[i] << " is not manifold\n";
      continue;
    }
    std::cout << "Checking " << names[i] << "... ";

    // Convert the mesh result to a soup for comparison.
    Mesh M;
    PMP::polygon_soup_to_polygon_mesh(soup_points[i], soup_polygons[i], M);

    std::vector<face_descriptor> m1_only;
    std::vector<face_descriptor> m2_only;
    PMP::match_faces(*(corefine_outputs[i].value()), M, boost::make_function_output_iterator([](const auto&) {}),
                                                        std::back_inserter(m1_only),
                                                        std::back_inserter(m2_only));
    const bool equal = m1_only.empty() && m2_only.empty();

    if (!equal)
    {
      std::cerr << "FAILED: results are different\n";
      assert(false);
    }
    else
      std::cout << " OK\n";
  }
}

int main(int argc, char** argv)
{
  const std::filesystem::path program_path(argv[0]);

  if (argc < 3)
  {
    std::cerr << "Usage: " << program_path.filename().string() << " P.off Q.off\n";
    return EXIT_FAILURE;
  }

  if(argc == 3)
  {
    test(argv[1], argv[2]);
    return EXIT_SUCCESS;
  }

  if ((argc-1)%2==0)
  {
    // cmd mode
    for (int i=1; i<argc; i+=2)
      test(argv[i], argv[i+1]);
    return EXIT_SUCCESS;
  }

  std::cerr << "An odd number of arguments \n";
  return EXIT_FAILURE;
}