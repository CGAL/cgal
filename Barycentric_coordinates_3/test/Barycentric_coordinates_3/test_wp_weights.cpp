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
  using Vector_3 = typename Kernel::Vector_3;
  using Mesh =  typename CGAL::Surface_mesh<Point_3>;

  // Cube
  Mesh cube;
  std::vector<Point_3> cube_coords;

  std::tie(cube, cube_coords) = tests::get_hexahedron<Kernel, Mesh>();
  CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Mesh, Kernel> wp_cube(cube);

  const FT step  = FT(1) / FT(100);
  const FT scale = FT(100);
  const FT limit = step*scale;
  const FT tol = tests::get_tolerance<FT>();

  std::vector<FT> wp_coordinates_cube;
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

        //Test unit partition, coordinates between 0 and 1 and Linear combination of vertices
        FT sum = FT(0);
        FT x_linear_comb = FT(0);
        FT y_linear_comb = FT(0);
        FT z_linear_comb = FT(0);

        for(std::size_t i = 0; i < 8; i++){

          sum += wp_coordinates_cube[i + count];
          x_linear_comb += cube_coords[i][0] * wp_coordinates_cube[i + count];
          y_linear_comb += cube_coords[i][1] * wp_coordinates_cube[i + count];
          z_linear_comb += cube_coords[i][2] * wp_coordinates_cube[i + count];

          assert(wp_coordinates_cube[i + count] >= FT(0) &&
                 wp_coordinates_cube[i + count] <= FT(1));
        }

        assert(CGAL::abs(FT(1) - sum) < tol);
        assert(CGAL::abs(x_linear_comb - query.x()) < tol &&
               CGAL::abs(y_linear_comb - query.y()) < tol &&
               CGAL::abs(z_linear_comb - query.z()) < tol);


        count += 8;
      }
    }
  }

  count = 0;

  /*
  Mesh tetra;
  std::vector<Point_3> coords;
  std::vector<Point_3> points;
  std::tie(tetra, coords) = tests::get_irregular_tetrahedron<Kernel, Mesh>();

  tests::sample_random_inside_tetrahedron(coords[0], coords[1], coords[2], coords[3],
   std::back_inserter(points), 6);
  */


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

  /*
  std::cout << "EPECK test :" << std::endl;
  test_overloads<EPECK>();
  std::cout << "EPECK PASSED" << std::endl;
  */

  return EXIT_SUCCESS;
}