#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/Wachspress_coordinates_3.h>

#include "include/utils.h"

// Typedefs.
using SCKER = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

template<typename Kernel>
void test_overloads() {

  using FT      = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Mesh =  typename CGAL::Surface_mesh<Point_3>;

  // Meshes
  Mesh cube;
  Mesh prism;

  cube = tests::get_hexahedron<Kernel, Mesh>();
  prism = tests::get_regular_prism<Kernel, Mesh>();

  CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Mesh, Kernel> wp_cube(cube);
  CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Mesh, Kernel> wp_prism(prism);

  const FT step  = FT(1) / FT(100);
  const FT scale = FT(100);
  const FT limit = step*scale;
  const FT tol = tests::get_tolerance<FT>();

  std::vector<FT> wp_coordinates_cube;
  std::vector<FT> wp_coordinates_prism;

  std::size_t count = 0;

  // Test cube
  //Check for barycenter
  wp_cube(Point_3(0.5, 0.5, 0.5), std::back_inserter(wp_coordinates_cube));
  for(std::size_t i = 0; i < 8; i++)
    assert(wp_coordinates_cube[i + count] == FT(0.125));
  count += 8;

  // Sample interior points
  for(FT x = step; x < limit; x += step){
    for(FT y = step; y < limit; y += step){
      for(FT z = step; z < limit; z += step){

        const Point_3 query(x, y, z);

        wp_cube(query, std::back_inserter(wp_coordinates_cube));

        //Test unit partition and coordinates between 0 and 1
        FT sum = FT(0);
        for(std::size_t i = 0; i < 8; i++){

          sum += wp_coordinates_cube[i + count];
          assert(wp_coordinates_cube[i + count] >= FT(0) &&
                 wp_coordinates_cube[i + count] <= FT(1));
        }

        //NOT EQUAL TO 1 WHEN (0.01, 0.01, 0.01)
        assert(CGAL::abs(FT(1) - sum) < tol);

        count += 8;
      }
    }
  }

  count = 0;

}

int main(){

  // Set cout precision
  std::cout.precision(20);

  std::cout << "SCKER test :" << std::endl;
  test_overloads<SCKER>();
  std::cout << "SCKER PASSED :" << std::endl;

  std::cout << "EPICK test :" << std::endl;
  test_overloads<EPICK>();
  std::cout << "EPICK PASSED :" << std::endl;

  //test_overloads<EPECK>();

  return EXIT_SUCCESS;
}