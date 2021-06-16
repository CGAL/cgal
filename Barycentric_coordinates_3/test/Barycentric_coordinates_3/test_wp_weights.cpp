#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

// Typedefs.
using SCKER = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

template<typename Kernel>
void test_overloads() {

  using FT      = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Mesh =  typename CGAL::Surface_mesh<Point_3>;

  // Meshes
  Mesh regular_terahedron;
  Mesh irregular_terahedron;
  Mesh cube;
  Mesh prism;

  // Vertices
  const Point_3 p0(0.0, 0.0, 0.0);
  const Point_3 p1(1.0, 0.0, 0.0);
  const Point_3 p2(0.0, 1.0, 0.0);
  const Point_3 p3(0.0, 0.0, 1.0);
  const Point_3 p4(1.0, 1.0, 0.0);
  const Point_3 p5(1.0, 0.0, 1.0);
  const Point_3 p6(0.0, 1.0, 1.0);
  const Point_3 p7(1.0, 1.0, 1.0);


  //Construct polyhedron
  CGAL::make_tetrahedron(p0, p1, p2, p3, irregular_terahedron);

}

int main(){

  test_overloads<SCKER>();
  test_overloads<EPICK>();
  test_overloads<EPECK>();

  return EXIT_SUCCESS;
}