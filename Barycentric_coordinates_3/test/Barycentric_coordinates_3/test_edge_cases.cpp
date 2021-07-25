#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3.h>

#include "include/utils.h"

// Typedefs.
using SCKER = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
using CP3 = CGAL::Barycentric_coordinates::Computation_policy_3;

template<typename Kernel, typename BC>
void test_coordinate() {

  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Mesh = typename CGAL::Surface_mesh<Point_3>;

  // tetrahedron
  Mesh tetrahedron;
  std::vector<Point_3> tetrahedron_coords;
  std::tie(tetrahedron, tetrahedron_coords) = tests::get_regular_tetrahedron<Kernel, Mesh>(FT(1.0));

  std::vector<Point_3> sample_points;
  std::vector<Point_3> sampled_tetrahedron_coords;
  const FT tol = tests::get_tolerance<FT>();
  int num_interior_samples = 0;

  for(auto& offset : {FT(0), tol, -tol}){

    // Sample interior and boundary
    std::tie(std::ignore, sampled_tetrahedron_coords) =
      tests::get_regular_tetrahedron<Kernel, Mesh>(FT(1.0) + offset);

    tests::random_points_tetrahedron<Kernel>(sampled_tetrahedron_coords,
      std::back_inserter(sample_points), 100);
    num_interior_samples += 100;

    // Add points close to vertices
    if(offset != tol){

      for(auto& v : sampled_tetrahedron_coords){

        sample_points.push_back(v);
        num_interior_samples++;
      }
    }
  }

  // Exterior points
  std::tie(std::ignore, sampled_tetrahedron_coords) =
    tests::get_regular_tetrahedron<Kernel, Mesh>(FT(2.0));

  tests::random_points_tetrahedron<Kernel>(sampled_tetrahedron_coords,
    std::back_inserter(sample_points), 100);

  // BC coordinates
  BC bc_tetrahedron(tetrahedron, CP3::WITH_EDGE_CASES);

  std::vector<FT> bc_coordinates_tetrahedron;
  int num_samples = 0;

  for(auto& point : sample_points){

    const Point_3 query(point.x(), point.y(), point.z());
    bc_coordinates_tetrahedron.clear();
    bc_tetrahedron(query, std::back_inserter(bc_coordinates_tetrahedron));

    tests::test_linear_precision<Kernel>(bc_coordinates_tetrahedron, tetrahedron_coords, query);
    tests::test_partition_of_unity<Kernel>(bc_coordinates_tetrahedron);

    if(num_samples < num_interior_samples)
      tests::test_positivity<Kernel>(bc_coordinates_tetrahedron);

    num_samples++;
  }
}

template<typename Kernel>
void test_overloads(){

  using Point_3 = typename Kernel::Point_3;
  using Mesh = typename CGAL::Surface_mesh<Point_3>;
  using WP = CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Mesh, Kernel>;
  using DH = CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_3<Mesh, Kernel>;
  using MV = CGAL::Barycentric_coordinates::Mean_value_coordinates_3<Mesh, Kernel>;

  std::cout << "WP test: " << std::endl;
  test_coordinate<Kernel, WP>();
  std::cout << "WP passed" << std::endl;

  std::cout << "DH test: " << std::endl;
  test_coordinate<Kernel, DH>();
  std::cout << "DH passed" << std::endl;

  std::cout << "MV test: " << std::endl;
  test_coordinate<Kernel, MV>();
  std::cout << "MV passed" << std::endl;
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
