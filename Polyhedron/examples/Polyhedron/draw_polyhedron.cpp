#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/draw_polyhedron.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef CGAL::Polyhedron_3<Kernel>                       Polyhedron;

int main(int argc, char* argv[])
{
  Polyhedron P;
  std::ifstream in1((argc>1)?argv[1]:"data/cross.off");
  in1 >> P;
  CGAL::draw(P);

  return EXIT_SUCCESS;
}
