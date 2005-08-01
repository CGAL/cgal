// file: examples/Polyhedron/polyhedron_prog_simple.C

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_cgm.h>

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::Halfedge_handle        Halfedge_handle;

int main()
{
  Polyhedron p1, p2, p3;
  Halfedge_handle h1 = p1.make_tetrahedron();
  Halfedge_handle h2 = p2.make_tetrahedron();

  Polyhedral_cgm pcgm1(p1);
  Polyhedral_cgm pcgm2(p2);
  Polyhedral_cgm pcgm3.minkowski_sum(cgm0, cgm1);
  pcgm3.get(p3);
  std::cout << p3 << std::endl;
}
