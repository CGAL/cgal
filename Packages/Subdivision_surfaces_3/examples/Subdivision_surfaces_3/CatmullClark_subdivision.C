// file: examples/Subdivision_surfaces_3/CatmullClark_subdivision.C

#include <CGAL/Subdivision_surfaces_3.h>

#include <iostream>
#include <fstream>

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Cartesian<double>            Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;

using namespace std;
using namespace CGAL;

int main(int argc, char **argv) {
  if (argc != 3) { 
    cout << "Usage: CatmullClark_subdivision filename d" << endl; 
    cout << "       filename: the input mash (.off)" << endl; 
    cout << "       d: the depth of the subdivision (0 < d < 10)" << endl; 
    exit(1);
  }

  ifstream in(argv[1]);
  int d = argv[2][0] - '0';

  Polyhedron P;
  in >> P; // read the .off

  Subdivision_surfaces_3<Polyhedron>::CatmullClark_subdivision(P,d);

  cout << P; // write the .off
  
  return 0;
}
