#include <CGAL/leda_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include "sat_writer.h"

typedef leda_integer RT;
typedef CGAL::Homogeneous<RT> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3;
typedef CGAL::sat_writer<Nef_polyhedron_3> Sat_writer;
int main(int argc, char argv[]) {

  Nef_polyhedron_3 N;
  std::cin >> N;

  Sat_writer SW(std::cout,N);
  SW.print();
}
