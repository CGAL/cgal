#include <CGAL/leda_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include "tetrahedron_generator.h"
#include "sat_writer.h"

typedef leda_integer RT;
typedef CGAL::Homogeneous<RT> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3;
typedef CGAL::sat_writer<Nef_polyhedron_3> Sat_writer;
typedef tetrahedron_generator<Kernel> tgen;

int main(int argc, char* argv[]) {

  assert(argc < 7);
  
  int nx = argc>1 ? std::atoi(argv[1]) : 2;
  int ny = argc>2 ? std::atoi(argv[2]) : 2;
  int nz = argc>3 ? std::atoi(argv[3]) : 2;
  int s  = argc>4 ? std::atoi(argv[4]) : 5;

  std::ostringstream out;
  tgen t(out);
  //  if(argc>5)
  //    t.create_tetrahedra(nx,ny,nz,s,std::atoi(argv[5]));
  //  else 
    t.create_tetrahedra(nx,ny,nz,s);

  std::istringstream in(out.str());
  Nef_polyhedron_3 N;
  in >> N;
  
  Sat_writer SW(std::cout,N);
  SW.print();
}
