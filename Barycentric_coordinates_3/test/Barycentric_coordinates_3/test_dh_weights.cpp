#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Barycentric_coordinates_3/Discrete_harmonic_coordinates_3.h>

#include "include/utils.h"

// Typedefs.
using Kernel = CGAL::Simple_cartesian<double>;

using FT = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;
using Mesh = typename CGAL::Surface_mesh<Point_3>;

int main(){

  std::cout.precision(20);

  // Cube
  Mesh cube;
  std::vector<Point_3> cube_coords;

  std::tie(cube, cube_coords) = tests::get_hexahedron<Kernel, Mesh>();

  std::vector<FT> dh_coordinates_cube;
  dh_coordinates_cube.resize(8);

  CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_3<Mesh, Kernel> dh_cube(cube);
  dh_cube(Point_3(FT(1)/FT(2), FT(1)/FT(2), FT(1)/FT(2)), dh_coordinates_cube.begin());

  for(auto c : dh_coordinates_cube)
    std::cout << c << "\n";
  tests::test_barycenter<Kernel>(dh_coordinates_cube);
}

