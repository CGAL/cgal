#include <CGAL/leda_integer.h>
typedef leda_integer RT;

#include <CGAL/Homogeneous.h>
typedef CGAL::Homogeneous<RT> Kernel;

#include <CGAL/basic.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include "worst_case_sphere_map_generator.h"

typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Vector_3 Vector_3;
typedef Nef_polyhedron::Aff_transformation_3 Aff_transformation_3;
typedef Nef_polyhedron::Vertex_const_iterator Vertex_const_iterator;

bool cgal_nef3_timer_on = false;

int main(int argc, char* argv[]) {

  CGAL_assertion(argc==2);

  Nef_polyhedron N1,N2;
  std::ostringstream out;
  worst_case_sphere_map_generator<Nef_polyhedron> wgen(std::cout);
  wgen.print(std::atoi(argv[1]));
};
