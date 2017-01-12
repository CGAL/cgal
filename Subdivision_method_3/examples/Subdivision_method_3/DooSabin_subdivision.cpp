#include <CGAL/Simple_cartesian.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Subdivision_method_3.h>

#include <iostream>
#include <boost/lexical_cast.hpp>

typedef CGAL::Simple_cartesian<double>     Kernel;
//typedef CGAL::Polyhedron_3<Kernel>         PolygonMesh;
typedef CGAL::Surface_mesh<Kernel::Point_3> PolygonMesh;

using namespace std;
using namespace CGAL;

int main(int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: DooSabin_subdivision d < filename" << endl;
    cout << "       d: the depth of the subdivision" << endl;
    cout << "       filename: the input mesh (.off)" << endl;
    return 0;
  }

  int d = boost::lexical_cast<int>(argv[1]);

  PolygonMesh pmesh;
  cin >> pmesh;

  Subdivision_method_3::DooSabin_subdivision(pmesh,d);

  cout << pmesh;

  return 0;
}
