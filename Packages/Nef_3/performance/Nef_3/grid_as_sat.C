#include <CGAL/leda_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include "grid_generator.h"
#include "sat_writer.h"

typedef leda_integer RT;
typedef CGAL::Homogeneous<RT> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3;
typedef CGAL::sat_writer<Nef_polyhedron_3> Sat_writer;
typedef CGAL::grid_generator<Nef_polyhedron_3> ggen;

int main(int argc, char* argv[]) {

  assert(argc>1 && argc<6);
  
  int nx = argc>2 ? std::atoi(argv[2]) : 2;
  int ny = argc>3 ? std::atoi(argv[3]) : 2;
  int nz = argc>4 ? std::atoi(argv[4]) : 2;

  std::ifstream fin(argv[1]);
  Nef_polyhedron_3 Nin;
  fin >> Nin;
  
  std::ostringstream out;
  ggen g(out, Nin);
  g.print(nx,ny,nz);

  std::istringstream in(out.str());
  Nef_polyhedron_3 N;
  in >> N;
  
  Sat_writer SW(std::cout,N);
  SW.print();
}
