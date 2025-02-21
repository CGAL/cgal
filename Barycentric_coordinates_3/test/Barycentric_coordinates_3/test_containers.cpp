#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Barycentric_coordinates_3.h>
#include <CGAL/boost/graph/generators.h>

// Typedefs.
using SCKER = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

template<typename Kernel, typename Mesh, typename OuputContainer>
void test_container(const Mesh& mesh, OuputContainer& coordinates) {

  namespace BC = CGAL::Barycentric_coordinates;
  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;

  const FT h = FT(1) / FT(2);
  const Point_3 centroid(h, h, h);

  coordinates.clear();
  BC::Wachspress_coordinates_3<Mesh, Kernel> wp(mesh);
  wp(centroid, std::back_inserter(coordinates));
  assert(coordinates.front() >= FT(0));
  assert(coordinates.back()  >= FT(0));

  coordinates.clear();
  BC::Discrete_harmonic_coordinates_3<Mesh, Kernel> dh(mesh);
  dh(centroid, std::back_inserter(coordinates));
  assert(coordinates.front() >= FT(0));
  assert(coordinates.back()  >= FT(0));

  coordinates.clear();
  BC::Mean_value_coordinates_3<Mesh, Kernel> mv(mesh);
  mv(centroid, std::back_inserter(coordinates));
  assert(coordinates.front() >= FT(0));
  assert(coordinates.back()  >= FT(0));
}

template<typename Kernel>
void test_containers() {

  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;

  using Polyhedron = CGAL::Polyhedron_3<Kernel>;
  using Surface_mesh = CGAL::Surface_mesh<Point_3>;

  using LCoords = std::list<FT>;
  using VCoords = std::vector<FT>;

  Polyhedron polyhedron;
  Surface_mesh surface_mesh;

  const Point_3 p0(0, 0, 0);
  const Point_3 p1(1, 0, 0);
  const Point_3 p2(1, 1, 0);
  const Point_3 p3(0, 1, 0);

  const Point_3 p4(0, 1, 1);
  const Point_3 p5(0, 0, 1);
  const Point_3 p6(1, 0, 1);
  const Point_3 p7(1, 1, 1);

  CGAL::make_hexahedron(p0, p1, p2, p3, p4, p5, p6, p7, polyhedron, CGAL::parameters::do_not_triangulate_faces(false));
  CGAL::make_hexahedron(p0, p1, p2, p3, p4, p5, p6, p7, surface_mesh, CGAL::parameters::do_not_triangulate_faces(false));

  LCoords lcoords;
  VCoords vcoords;

  test_container<Kernel>(polyhedron, lcoords);
  test_container<Kernel>(polyhedron, vcoords);

  test_container<Kernel>(surface_mesh, lcoords);
  test_container<Kernel>(surface_mesh, vcoords);
}

int main() {

  std::cout << "SCKER test :" << std::endl;
  test_containers<SCKER>();
  std::cout << "SCKER PASSED" << std::endl;

  std::cout << "EPICK test :" << std::endl;
  test_containers<EPICK>();
  std::cout << "EPICK PASSED" << std::endl;

  std::cout << "EPECK test :" << std::endl;
  test_containers<EPECK>();
  std::cout << "EPECK PASSED" << std::endl;
}
