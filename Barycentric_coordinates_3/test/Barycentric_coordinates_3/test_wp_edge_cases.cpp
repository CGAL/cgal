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
  using Mesh = typename CGAL::Surface_mesh<Point_3>;

  // tetrahedron
  Mesh tetrahedron;
  std::vector<Point_3> tetrahedron_coords;
  std::vector<Point_3> sample_points;

  // Offsets
  const FT tol = tests::get_tolerance<FT>()/FT(10);

  // Sample interior and boundary
  std::tie(tetrahedron, tetrahedron_coords) = tests::get_regular_tetrahedron<Kernel, Mesh>(FT(1.0));
  tests::random_points_tetrahedron<Kernel>(tetrahedron_coords,
   std::back_inserter(sample_points), 100);

  // Face offsets
  std::vector<Point_3> tetrahedron_tol_diff;
  std::tie(std::ignore, tetrahedron_tol_diff) = tests::get_regular_tetrahedron<Kernel, Mesh>(FT(1.0) - tol);
  tests::random_points_tetrahedron<Kernel>(tetrahedron_tol_diff,
   std::back_inserter(sample_points), 100);

  // Vertice offsets
  for(auto& v : tetrahedron_coords){

    sample_points.push_back(Point_3(v.x(), v.y(), v.z()));
  }

  // WP coordinates
  CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Mesh, Kernel> wp_tetrahedron(
    tetrahedron, CP3::WITH_EDGE_CASES);

  std::vector<FT> wp_coordinates_tetrahedron;

  for(auto& point : sample_points){

    const Point_3 query(point.x(), point.y(), point.z());
    wp_coordinates_tetrahedron.clear();
    wp_tetrahedron(query, std::back_inserter(wp_coordinates_tetrahedron));

    tests::test_linear_precision<Kernel>(wp_coordinates_tetrahedron, tetrahedron_coords, query);
    tests::test_partition_of_unity<Kernel>(wp_coordinates_tetrahedron);
    tests::test_positivity<Kernel>(wp_coordinates_tetrahedron);
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
