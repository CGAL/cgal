#include <CGAL/leda_integer.h>
typedef leda_integer RT;

#include <CGAL/Homogeneous.h>
typedef CGAL::Homogeneous<RT> Kernel;

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Timer.h>
#include <sstream>
#include <fstream>
#include "grid_generator.h"

typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Aff_transformation_3 Aff_transformation_3;
typedef Nef_polyhedron::Vector_3 Vector_3;
typedef CGAL::grid_generator<Nef_polyhedron> ggen;

bool cgal_nef3_timer_on = false;

int main(int argc, char* argv[]) {

  assert(argc>1 && argc < 5);
  
  int nx = argc>2 ? std::atoi(argv[2]) : 2;
  int ny = argc>3 ? std::atoi(argv[3]) : nx;

  std::ifstream in(argv[1]);
  Nef_polyhedron Nin;
  in >> Nin;

  std::ostringstream out1;
  ggen g(out1, Nin,false);
  g.print(nx,1,1);
  std::istringstream in1(out1.str());
  Nef_polyhedron N1;
  in1 >> N1;
  CGAL_assertion(N1.is_valid());

  Nin.transform(Aff_transformation_3(0,-1,0,
				     1,0,0,
				     0,0,1,1));

  std::ostringstream out2;
  ggen g2(out2, Nin,false);
  g2.print(1,ny,1);
  std::istringstream in2(out2.str());
  Nef_polyhedron N2;
  in2 >> N2;
  CGAL_assertion(N2.is_valid());

  N1.transform(Aff_transformation_3(328354,0,0,
				    0,328304,-5730,
				    0,5730,328304,328354));
  N2.transform(Aff_transformation_3(328354,0,0,
				    0,328304,-5730,
				    0,5730,328304,328354));

  cgal_nef3_timer_on = true;

  N1 = N1.join(N2);
}
