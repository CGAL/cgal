#include <CGAL/Simple_cartesian.h>
#include <CGAL/Subdivision_method_3.h>

#include <iostream>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;

using namespace std;
using namespace CGAL;

int main(int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: Loop_subdivision d < filename" << endl;
    cout << "       d: the depth of the subdivision (0 < d < 10)" << endl;
    cout << "       filename: the input mesh (.off)" << endl;
    return 0;
  }

  int d = argv[1][0] - '0';

  Polyhedron P;
  cin >> P; // read the .off

  Subdivision_method_3::Loop_subdivision(P,d);

  cout << P; // write the .off

  return 0;
}
