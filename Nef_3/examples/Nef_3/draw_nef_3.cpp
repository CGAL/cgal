#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/draw_nef_3.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel   Kernel;
typedef CGAL::Polyhedron_3<Kernel>                          Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel>                      Nef_polyhedron;

int main(int argc, char *argv[])
{
  // read OFF file into a polyhedron
  Polyhedron P;
  std::ifstream ifs((argc > 1) ? argv[1] : "data/cross.off");
  ifs >> P;

  // initialize nef from polyhedron
  Nef_polyhedron N(P);

  // draw Nef Polyhedron
  CGAL::draw(N);

  return EXIT_SUCCESS;
}
