#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Subdivision_method_3.h>

#include <iostream>
#include <boost/lexical_cast.hpp>

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> PolygonMesh;

using namespace std;
using namespace CGAL;

int main(int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: Loop_subdivision d < filename" << endl;
    cout << "       d: the depth of the subdivision" << endl;
    cout << "       filename: the input mesh (.off)" << endl;
    return 0;
  }

  unsigned int d = boost::lexical_cast<unsigned int>(argv[1]);

  PolygonMesh pmesh;
  cin >> pmesh;

  Subdivision_method_3::Loop_subdivision(pmesh,d);

  cout << pmesh;

  return 0;
}
