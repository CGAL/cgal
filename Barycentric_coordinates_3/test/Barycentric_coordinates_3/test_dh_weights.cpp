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
  Mesh tetrahedron;
  std::vector<Point_3> tetrahedron_coords;

  std::tie(tetrahedron, tetrahedron_coords) = tests::get_irregular_tetrahedron<Kernel, Mesh>();

  std::vector<FT> dh_coordinates_tetrahedron;
  dh_coordinates_tetrahedron.resize(4);

  CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_3<Mesh, Kernel> dh_tetrahedron(tetrahedron);
  dh_tetrahedron(Point_3(FT(1)/FT(4), FT(1)/FT(4), FT(1)/FT(4)), dh_coordinates_tetrahedron.begin());

  for(auto c : dh_coordinates_tetrahedron)
    std::cout << c << "\n";
  tests::test_barycenter<Kernel>(dh_coordinates_tetrahedron);
}

