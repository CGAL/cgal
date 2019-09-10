#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/subdivision_method_3.h>
#include <CGAL/Timer.h>

#include <boost/lexical_cast.hpp>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>         Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3>    PolygonMesh;

using namespace std;
using namespace CGAL;
namespace params = CGAL::parameters;

int main(int argc, char** argv) {
  if (argc > 4) {
    cerr << "Usage: CatmullClark_subdivision [d] [filename_in] [filename_out] \n";
    cerr << "         d -- the depth of the subdivision (default: 1) \n";
    cerr << "         filename_in -- the input mesh (.off) (default: data/quint_tris.off) \n";
    cerr << "         filename_out -- the output mesh (.off) (default: result.off)" << endl;
    return 1;
  }

  int d = (argc > 1) ? boost::lexical_cast<int>(argv[1]) : 1;
  const char* in_file = (argc > 2) ? argv[2] : "data/quint_tris.off";
  const char* out_file = (argc > 3) ? argv[3] : "result.off";

  PolygonMesh pmesh;
  std::ifstream in(in_file);
  if(in.fail()) {
    std::cerr << "Could not open input file " << in_file << std::endl;
    return 1;
  }
  in >> pmesh;

  Timer t;
  t.start();
  Subdivision_method_3::CatmullClark_subdivision(pmesh, params::number_of_iterations(d));
  std::cerr << "Done (" << t.time() << " s)" << std::endl;

  std::ofstream out(out_file);
  out << pmesh;

  return 0;
}
