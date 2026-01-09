#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/Mean_value_coordinates_3.h>

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
  std::vector<Point_3> sample_points;

  std::tie(tetrahedron, tetrahedron_coords) = tests::get_irregular_tetrahedron<Kernel, Mesh>();
  CGAL::Barycentric_coordinates::Mean_value_coordinates_3<Mesh, Kernel> mv_tetrahedron(tetrahedron);
  tests::random_points_tetrahedron<Kernel>(tetrahedron_coords, std::back_inserter(sample_points), 1000);

  std::vector<FT> mv_coordinates_tetrahedron;

  //Check for barycenter
  mv_tetrahedron(Point_3(FT(0.25), FT(0.25), FT(0.25)),
   std::back_inserter(mv_coordinates_tetrahedron));
  tests::test_barycenter<Kernel>(mv_coordinates_tetrahedron);

  for(auto& point : sample_points){

    FT x = point.x();
    FT y = point.y();
    FT z = point.z();

    if(x + y + z == FT(1) || x == FT(0) || y == FT(0) || z == FT(0)) //avoid points inside faces
      continue;

    const Point_3 query(x, y, z);
    mv_coordinates_tetrahedron.clear();
    mv_tetrahedron(query, std::back_inserter(mv_coordinates_tetrahedron));

    tests::test_linear_precision<Kernel>(mv_coordinates_tetrahedron, tetrahedron_coords, query);
    tests::test_partition_of_unity<Kernel>(mv_coordinates_tetrahedron);
    tests::test_positivity<Kernel>(mv_coordinates_tetrahedron);

    mv_coordinates_tetrahedron.clear();
    CGAL::Barycentric_coordinates::mean_value_coordinates_3(
      tetrahedron, query, std::back_inserter(mv_coordinates_tetrahedron));
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
