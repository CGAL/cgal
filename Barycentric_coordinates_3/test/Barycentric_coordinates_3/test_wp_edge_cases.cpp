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
  Mesh cube;
  std::vector<Point_3> cube_coords;

  std::tie(cube, cube_coords) = tests::get_hexahedron<Kernel, Mesh>();

  CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Mesh, Kernel> wp_cube(cube, CP3::WITH_EDGE_CASES);

  const FT step  = FT(1) / FT(10);
  const FT scale = FT(10);
  const FT limit = step*scale;

  std::vector<FT> wp_coordinates_cube;
  wp_coordinates_cube.resize(8);

  // Test cube
  //Check for barycenter
  wp_cube(Point_3(FT(1)/FT(2), FT(1)/FT(2), FT(1)/FT(2)), wp_coordinates_cube.begin());
  tests::test_barycenter<Kernel>(wp_coordinates_cube);

  // Sample interior points
  for(FT x = step; x < limit; x += step){
    for(FT y = step; y < limit; y += step){
      for(FT z = step; z < limit; z += step){

        const Point_3 query(x, y, z);
        wp_cube(query, wp_coordinates_cube.begin());

        tests::test_linear_precision<Kernel>(wp_coordinates_cube, cube_coords, query);
        tests::test_partition_of_unity<Kernel>(wp_coordinates_cube);
        tests::test_positivity<Kernel>(wp_coordinates_cube);

        CGAL::Barycentric_coordinates::wachspress_coordinates_3(
          cube, query, wp_coordinates_cube.begin());
      }
    }
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