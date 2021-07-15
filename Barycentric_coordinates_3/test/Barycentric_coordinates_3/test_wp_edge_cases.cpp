#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/Wachspress_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/barycentric_enum_3.h>

#include "include/utils.h"

// Typedefs.
using SCKER = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
using CP3 = CGAL::Barycentric_coordinates::Computation_policy_3;

template<typename Kernel>
void test_overloads() {

  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Mesh = CGAL::Surface_mesh<Point_3>;

  // Cube
  Mesh tetrahedron;
  std::vector<Point_3> tetrahedron_coords;

  std::tie(tetrahedron, tetrahedron_coords) = tests::get_irregular_tetrahedron<Kernel, Mesh>();

  CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Mesh, Kernel> wp_tetrahedron(tetrahedron, CP3::WITH_EDGE_CASES);

  std::vector<FT> wp_coordinates_tetrahedron;
  wp_coordinates_tetrahedron.resize(4);

  wp_tetrahedron(Point_3(0.05, 0.1, 0.85), wp_coordinates_tetrahedron.begin());

  for(auto u : wp_coordinates_tetrahedron){
    std::cout << u << " \n";
  }

}

int main(){

  // Set cout precision
  std::cout.precision(20);

  std::cout << "SCKER test :" << std::endl;
  test_overloads<SCKER>();
  std::cout << "SCKER PASSED" << std::endl;

  /*
  std::cout << "EPICK test :" << std::endl;
  test_overloads<EPICK>();
  std::cout << "EPICK PASSED" << std::endl;

  std::cout << "EPECK test :" << std::endl;
  test_overloads<EPECK>();
  std::cout << "EPECK PASSED" << std::endl;
  */

  return EXIT_SUCCESS;
}