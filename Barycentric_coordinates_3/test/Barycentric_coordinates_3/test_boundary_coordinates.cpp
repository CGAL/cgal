#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/boundary_coordinates_3.h>

#include "include/utils.h"

// Typedefs.
using SCKER = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

template<typename Kernel>
void test_overloads(){

  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Mesh = typename CGAL::Surface_mesh<Point_3>;

  // tetrahedron
  Mesh tetrahedron;
  std::vector<Point_3> vertices;
  std::vector<Point_3> sample_points;
  std::vector<FT> tetra_coords(4);

  // 500 interior points + 500 surface points
  std::tie(tetrahedron, vertices) = tests::get_irregular_tetrahedron<Kernel, Mesh>();
  tests::random_points_tetrahedron<Kernel>(vertices, std::back_inserter(sample_points), 1000);

  std::size_t num_samples = 0;
  for(auto& point : sample_points){

    const FT x = point.x(), y = point.y(), z = point.z();
    const Point_3 query(x, y, z);
    CGAL::Barycentric_coordinates::boundary_coordinates_3(tetrahedron, query, tetra_coords.begin());

    if(num_samples < 500){
      assert(tetra_coords[0] == FT(0) &&
             tetra_coords[1] == FT(0) &&
             tetra_coords[2] == FT(0) &&
             tetra_coords[3] == FT(0));
    }
    else{
      assert(CGAL::abs(1-x-y-z - tetra_coords[0]) == FT(0) &&
             CGAL::abs(x - tetra_coords[1]) == FT(0) &&
             CGAL::abs(y - tetra_coords[2]) == FT(0) &&
             CGAL::abs(z - tetra_coords[3]) == FT(0));
    }

    num_samples++;
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
