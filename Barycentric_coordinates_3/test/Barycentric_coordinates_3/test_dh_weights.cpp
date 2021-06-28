#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/Discrete_harmonic_coordinates_3.h>

#include "include/utils.h"

// Typedefs.
using SCKER = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

template<typename Kernel>
void test_overloads() {

  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Mesh = typename CGAL::Surface_mesh<Point_3>;

  // tetrahedron
  Mesh tetrahedron;
  std::vector<Point_3> tetrahedron_coords;

  std::tie(tetrahedron, tetrahedron_coords) = tests::get_irregular_tetrahedron<Kernel, Mesh>();
  CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_3<Mesh, Kernel> dh_tetrahedron(tetrahedron);

  const FT step  = FT(1) / FT(10);
  const FT scale = FT(10);
  const FT limit = step*scale;

  std::vector<FT> dh_coordinates_tetrahedron;
  dh_coordinates_tetrahedron.resize(4);

  // Test cube
  //Check for barycenter
  dh_tetrahedron(Point_3(FT(1)/FT(4), FT(1)/FT(4), FT(1)/FT(4)), dh_coordinates_tetrahedron.begin());
  tests::test_barycenter<Kernel>(dh_coordinates_tetrahedron);

  // Sample interior points
  for(FT x = step; x < limit; x += step){
    for(FT y = step; y < limit; y += step){
      for(FT z = step; z < FT(1) - x - y - step; z+= step){ // Excludes points inside faces

        dh_coordinates_tetrahedron.clear();
        dh_coordinates_tetrahedron.resize(4);

        const Point_3 query(x, y, z);
        dh_tetrahedron(query, dh_coordinates_tetrahedron.begin());

        tests::test_linear_precision<Kernel>(dh_coordinates_tetrahedron, tetrahedron_coords, query);
        tests::test_partition_of_unity<Kernel>(dh_coordinates_tetrahedron);
        tests::test_positivity<Kernel>(dh_coordinates_tetrahedron);
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

  std::cout << "EPICK test :" << std::endl;
  test_overloads<EPICK>();
  std::cout << "EPICK PASSED" << std::endl;

  std::cout << "EPECK test :" << std::endl;
  test_overloads<EPECK>();
  std::cout << "EPECK PASSED" << std::endl;

  return EXIT_SUCCESS;
}
