
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Polyhedral_cgm_polyhedron_3.h>
#include <CGAL/Polyhedral_cgm.h>
#include <CGAL/IO/Polyhedral_cgm_iostream.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>

typedef CGAL::Gmpq                              Number_type;
typedef CGAL::Simple_cartesian<Number_type>     Kernel;
typedef CGAL::Polyhedral_cgm<Kernel>            Polyhedral_cgm;
typedef CGAL::Polyhedral_cgm_polyhedron_3<Polyhedral_cgm>
                                                Polyhedron;
typedef Polyhedron::Halfedge_handle             Halfedge_handle;
typedef CGAL::Polyhedral_cgm_initializer<Polyhedral_cgm>
                                                Polyhedral_cgm_initializer;
typedef Kernel::Point_3                         Point_3;

int main()
{
  Point_3 p( 1.0, 0.0, 0.0);
  Point_3 q( 0.0, 1.0, 0.0);
  Point_3 r( 0.0, 0.0, 1.0);
  Point_3 s( 0.0, 0.0, 0.0);

  Polyhedron P1, P2;
  P1.make_tetrahedron(p, q, r, s);
  P2.make_tetrahedron(p, q, r, s);

  Polyhedral_cgm pcgm1;
  Polyhedral_cgm_initializer pcgm_initializer1(pcgm1);
  pcgm_initializer1(P1);
  std::cout << pcgm1.number_of_vertices() << " "
            << pcgm1.number_of_edges() << " "
            << pcgm1.number_of_facets()
            << std::endl;

  Polyhedral_cgm pcgm2;
  Polyhedral_cgm_initializer pcgm_initializer2(pcgm2);
  pcgm_initializer2(P2);
  std::cout << pcgm2.number_of_vertices() << " "
            << pcgm2.number_of_edges() << " "
            << pcgm2.number_of_facets()
            << std::endl;

  Polyhedral_cgm pcgm3;
  pcgm3.minkowski_sum(&pcgm1, &pcgm2);
  std::cout << pcgm3.number_of_vertices() << " "
            << pcgm3.number_of_edges() << " "
            << pcgm3.number_of_facets()
            << std::endl;
  // std::cout << pcgm3 << std::endl;
  // pcgm3.get(p3);
  // std::cout << p3 << std::endl;
}
