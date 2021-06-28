#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/Wachspress_coordinates_3.h>
#include <CGAL/Barycentric_coordinates_3/Discrete_harmonic_coordinates_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

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
  Mesh cube;
  std::vector<Point_3> vertices;
  std::tie(cube, vertices) = tests::get_hexahedron<Kernel, Mesh>();

  //CGAL::Polygon_mesh_processing::triangulate_faces(cube);

  CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Mesh, Kernel> ws(cube);
  CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_3<Mesh, Kernel> dh(cube);

  std::vector<FT> cube_coordinates;
  std::vector<FT> wp_coordinates;
  std::vector<FT> dh_coordinates;

  // Sample points
  const FT step  = FT(1) / FT(10);
  const FT scale = FT(10);
  const FT limit = step*scale;
  const FT tol = tests::get_tolerance<FT>();

  std::size_t count = 0;

  for(FT x = step; x < limit - step; x += step){
    for(FT y = step; y < limit - step; y += step){
      for(FT z = step; z < limit - step; z+= step){

        const Point_3 query = Point_3(x, y, z);

        wp_coordinates.clear();
        dh_coordinates.clear();
        dh_coordinates.resize(8);
        wp_coordinates.resize(8);

        ws(query, wp_coordinates.begin());
        dh(query, dh_coordinates.begin());

        tests::test_positivity<Kernel>(wp_coordinates);
        tests::test_positivity<Kernel>(dh_coordinates);

        tests::test_partition_of_unity<Kernel>(wp_coordinates);
        tests::test_partition_of_unity<Kernel>(dh_coordinates);

        tests::test_linear_precision<Kernel>(wp_coordinates, vertices, query);
        tests::test_linear_precision<Kernel>(dh_coordinates, vertices, query);

        for(std::size_t i = 0; i < 8; i++){
          assert(CGAL::abs(dh_coordinates[count + i] - wp_coordinates[count + i]) < tol);
        }

      }
    }
  }

  std::cout << "test_dh_wp: PASSED" << std::endl;
  return EXIT_SUCCESS;
}
