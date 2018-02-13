#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/draw_polyhedron.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef CGAL::Polyhedron_3<Kernel>                       Polyhedron;

int main()
{
  std::cerr<<"Loading OFF file ... "<<std::endl;
  Polyhedron P;
  std::cin>>P;

  CGAL::draw(P);

  return EXIT_SUCCESS;
}
