#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/Wachspress_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/Discrete_harmonic_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/Mean_value_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/tetrahedron_coordinates.h>

#include "include/utils.h"

//Typedefs
using Kernel = CGAL::Simple_cartesian<double>;

using FT = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;
using Mesh = typename CGAL::Surface_mesh<Point_3>;

using WPC3 = CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Mesh, Kernel>;
using MVC3 = CGAL::Barycentric_coordinates::Mean_value_coordinates_3<Mesh, Kernel>;
using DHC3 = CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_3<Mesh, Kernel>;

template <typename COORD>
void test_overload(){

  // Regular tetrahedron
  Mesh tetrahedron;
  std::vector<Point_3> vertices;
  std::tie(tetrahedron, vertices) = tests::get_irregular_tetrahedron<Kernel, Mesh>();

  COORD bar(tetrahedron);

  std::vector<FT> tetra_coordinates;
  std::vector<FT> bar_coordinates;

  // Sample points
  const FT step  = FT(1) / FT(20);
  const FT scale = FT(10);
  const FT tol = tests::get_tolerance<FT>();

  std::size_t count = 0;
  const FT limit = scale * step;

  for(FT x = step; x < limit; x += step){
    for(FT y = step; y < limit; y += step){
      for(FT z = step; z < FT(1) - x - y - step; z+= step){ // Excludes points inside faces

        const Point_3 query = Point_3(x, y, z);

        bar_coordinates.clear();
        tetra_coordinates.clear();
        bar(query, std::back_inserter(bar_coordinates));
        CGAL::Barycentric_coordinates::tetrahedron_coordinates(vertices[0], vertices[1],
         vertices[2], vertices[3], query, std::back_inserter(tetra_coordinates));

        tests::test_positivity<Kernel>(bar_coordinates);
        tests::test_positivity<Kernel>(tetra_coordinates);

        tests::test_partition_of_unity<Kernel>(bar_coordinates);
        tests::test_partition_of_unity<Kernel>(tetra_coordinates);

        tests::test_linear_precision<Kernel>(bar_coordinates, vertices, query);
        tests::test_linear_precision<Kernel>(tetra_coordinates, vertices, query);

        assert(
          CGAL::abs(tetra_coordinates[count + 0] - bar_coordinates[count + 0]) < tol &&
          CGAL::abs(tetra_coordinates[count + 1] - bar_coordinates[count + 1]) < tol &&
          CGAL::abs(tetra_coordinates[count + 2] - bar_coordinates[count + 2]) < tol &&
          CGAL::abs(tetra_coordinates[count + 3] - bar_coordinates[count + 3]) < tol);
      }
    }
  }

}

int main(){

  // Set cout precision
  std::cout.precision(20);

  std::cout << "Wachspress: " << std::endl;
  test_overload<WPC3>();
  std::cout << "Wachspress_tetrahedron PASSED" << std::endl;

  std::cout << "Discrete Harmonic: " << std::endl;
  test_overload<DHC3>();
  std::cout << "Discrete_harmonic_tetrahedron PASSED" << std::endl;

  std::cout << "Mean Value: " << std::endl;
  test_overload<MVC3>();
  std::cout << "Mean_value_tetrahedron PASSED" << std::endl;

  std::cout << "test_wp_dh_mv_tetrahedron: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
