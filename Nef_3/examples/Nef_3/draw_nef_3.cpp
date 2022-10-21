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
  Polyhedron P1, P2;
  std::ifstream ifs1((argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cross_quad.off"));
  ifs1 >> P1;
  std::ifstream ifs2((argc > 1) ? argv[1] : CGAL::data_file_path("meshes/beam.off"));
  ifs2 >> P2;

  // initialize nef from polyhedron
  Nef_polyhedron N1(P1);
  Nef_polyhedron N2(P2);

  // draw Nef Polyhedron
  CGAL::draw(N1-N2);

  return EXIT_SUCCESS;
}
