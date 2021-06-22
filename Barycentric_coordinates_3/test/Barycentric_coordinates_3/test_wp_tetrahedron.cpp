#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/Wachspress_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/tetrahedron_coordinates.h>

#include "include/utils.h"

//Typedefs
using Kernel = CGAL::Simple_cartesian<double>;

int main(){

  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Mesh = typename CGAL::Surface_mesh<Point_3>;

  // Set cout precision
  std::cout.precision(20);

  // Regular tetrahedron
  Mesh tetrahedron;
  std::vector<Point_3> vertices;
  std::tie(tetrahedron, vertices) = tests::get_irregular_tetrahedron<Kernel, Mesh>();

  CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Mesh, Kernel> ws(tetrahedron);

  std::vector<FT> tetra_coordinates;
  std::vector<FT> wp_coordinates;

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

        wp_coordinates.clear();
        tetra_coordinates.clear();
        tetra_coordinates.resize(4);
        wp_coordinates.resize(4);

        ws(query, wp_coordinates.begin());
        CGAL::Barycentric_coordinates::tetrahedron_coordinates(vertices[0], vertices[1],
         vertices[2], vertices[3], query, tetra_coordinates.begin());

        tests::test_positivity<Kernel>(wp_coordinates);
        tests::test_positivity<Kernel>(tetra_coordinates);

        tests::test_partition_of_unity<Kernel>(wp_coordinates);
        tests::test_partition_of_unity<Kernel>(tetra_coordinates);

        tests::test_linear_precision<Kernel>(wp_coordinates, vertices, query);
        tests::test_linear_precision<Kernel>(tetra_coordinates, vertices, query);

        assert(
          CGAL::abs(tetra_coordinates[count + 0] - wp_coordinates[count + 0]) < tol &&
          CGAL::abs(tetra_coordinates[count + 1] - wp_coordinates[count + 1]) < tol &&
          CGAL::abs(tetra_coordinates[count + 2] - wp_coordinates[count + 2]) < tol &&
          CGAL::abs(tetra_coordinates[count + 3] - wp_coordinates[count + 3]) < tol);
      }
    }
  }

  std::cout << "test_wp_tetrahedron: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
