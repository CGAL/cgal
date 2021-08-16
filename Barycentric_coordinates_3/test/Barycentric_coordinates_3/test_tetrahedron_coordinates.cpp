#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/tetrahedron_coordinates.h>

#include "include/utils.h"

// Typedefs.
using SCKER = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

template <typename Kernel>
void test_overloads(){

  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Mesh = typename CGAL::Surface_mesh<Point_3>;

  // tetrahedron
  Mesh tetrahedron;
  std::vector<Point_3> vertices;
  std::vector<Point_3> sample_points;
  std::vector<FT> tetra_coords(4);

  //Interior and surface points
  std::tie(tetrahedron, vertices) = tests::get_irregular_tetrahedron<Kernel, Mesh>();
  tests::random_points_tetrahedron<Kernel>(vertices, std::back_inserter(sample_points), 1000);

  //Add vertices
  for(auto& v : vertices)
    sample_points.push_back(v);

  //Exterior + interior points
  std::vector<Point_3> ext_vertices = {Point_3(0.0, 0.0, 0.0), Point_3(5.0, 0.0, 0.0),
                                       Point_3(0.0, 5.0, 0.0), Point_3(0.0, 0.0, 5.0)};
  tests::random_points_tetrahedron<Kernel>(ext_vertices, std::back_inserter(sample_points), 1000);

  for(auto& point : sample_points){

    const FT x = point.x(), y = point.y(), z = point.z();
    const Point_3 query(x, y, z);
    CGAL::Barycentric_coordinates::tetrahedron_coordinates(vertices[0], vertices[1],
      vertices[2], vertices[3], query, tetra_coords.begin());

    const std::array<FT, 4> tetra_coords_array =
    CGAL::Barycentric_coordinates::tetrahedron_coordinates_in_array(vertices[0], vertices[1],
      vertices[2], vertices[3], query);

    assert(CGAL::abs(1-x-y-z - tetra_coords[0]) == FT(0) &&
           CGAL::abs(x - tetra_coords[1]) == FT(0) &&
           CGAL::abs(y - tetra_coords[2]) == FT(0) &&
           CGAL::abs(z - tetra_coords[3]) == FT(0));

    assert(tetra_coords_array[0] == tetra_coords[0] &&
           tetra_coords_array[1] == tetra_coords[1] &&
           tetra_coords_array[2] == tetra_coords[2] &&
           tetra_coords_array[3] == tetra_coords[3]);
  }
}

int main(){

  // Set cout precision
  std::cout.precision(20);

  std::cout << "SCKER test :" << std::endl;
  test_overloads<SCKER>();
  std::cout << "SCKER PASSED" << std::endl;

  std::cout << "EPICK test :" << std::endl;
  test_overloads<EPICK>();
  std::cout << "EPICK PASSED" << std::endl;

  std::cout << "EPECK test :" << std::endl;
  test_overloads<EPECK>();
  std::cout << "EPECK PASSED" << std::endl;

  return EXIT_SUCCESS;
}
