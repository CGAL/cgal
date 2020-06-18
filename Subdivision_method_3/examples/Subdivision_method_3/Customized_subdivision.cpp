#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/subdivision_method_3.h>
#include <CGAL/Timer.h>

#include <boost/lexical_cast.hpp>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> PolygonMesh;
using namespace std;
using namespace CGAL;
namespace params = CGAL::parameters;

// ======================================================================
template <class Poly>
class WLoop_mask_3 {
  typedef Poly                                         PolygonMesh;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<PolygonMesh, vertex_point_t>::type Vertex_pmap;
  typedef typename boost::property_traits<Vertex_pmap>::value_type Point;
  typedef typename boost::property_traits<Vertex_pmap>::reference Point_ref;

  PolygonMesh& pmesh;
  Vertex_pmap vpm;

public:
  WLoop_mask_3(PolygonMesh& pmesh)
    : pmesh(pmesh), vpm(get(CGAL::vertex_point, pmesh))
  {}

  void edge_node(halfedge_descriptor hd, Point& pt) {
    Point_ref p1 = get(vpm, target(hd,pmesh));
    Point_ref p2 = get(vpm, target(opposite(hd,pmesh),pmesh));
    Point_ref f1 = get(vpm, target(next(hd,pmesh),pmesh));
    Point_ref f2 = get(vpm, target(next(opposite(hd,pmesh),pmesh),pmesh));

    pt = Point((3*(p1[0]+p2[0])+f1[0]+f2[0])/8,
               (3*(p1[1]+p2[1])+f1[1]+f2[1])/8,
               (3*(p1[2]+p2[2])+f1[2]+f2[2])/8 );
  }
  void vertex_node(vertex_descriptor vd, Point& pt) {
    double R[] = {0.0, 0.0, 0.0};
    Point_ref S = get(vpm,vd);

    std::size_t n = 0;
    for(halfedge_descriptor hd : halfedges_around_target(vd, pmesh)){
      ++n;
      Point_ref p = get(vpm, target(opposite(hd,pmesh),pmesh));
      R[0] += p[0];         R[1] += p[1];         R[2] += p[2];
    }

    if (n == 6) {
      pt = Point((10*S[0]+R[0])/16, (10*S[1]+R[1])/16, (10*S[2]+R[2])/16);
    } else if (n == 3) {
      double B = (5.0/8.0 - std::sqrt(3+2*std::cos(6.283/n))/64.0)/n;
      double A = 1-n*B;
      pt = Point((A*S[0]+B*R[0]), (A*S[1]+B*R[1]), (A*S[2]+B*R[2]));
    } else {
      double B = 3.0/8.0/n;
      double A = 1-n*B;
      pt = Point((A*S[0]+B*R[0]), (A*S[1]+B*R[1]), (A*S[2]+B*R[2]));
    }
  }

  void border_node(halfedge_descriptor hd, Point& ept, Point& vpt) {
    Point_ref ep1 = get(vpm, target(hd,pmesh));
    Point_ref ep2 = get(vpm, target(opposite(hd,pmesh),pmesh));
    ept = Point((ep1[0]+ep2[0])/2, (ep1[1]+ep2[1])/2, (ep1[2]+ep2[2])/2);

    Halfedge_around_target_circulator<Poly> vcir(hd,pmesh);
    Point_ref vp1  = get(vpm, target(opposite(*vcir,pmesh),pmesh));
    Point_ref vp0  = get(vpm, target(*vcir,pmesh));
    --vcir;
    Point_ref vp_1 = get(vpm,target(opposite(*vcir,pmesh),pmesh));
    vpt = Point((vp_1[0] + 6*vp0[0] + vp1[0])/8,
                (vp_1[1] + 6*vp0[1] + vp1[1])/8,
                (vp_1[2] + 6*vp0[2] + vp1[2])/8 );
  }
};

int main(int argc, char **argv) {
  if (argc > 4) {
    cerr << "Usage: Customized_subdivision [d] [filename_in] [filename_out] \n";
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
  Subdivision_method_3::PTQ(pmesh, WLoop_mask_3<PolygonMesh>(pmesh), params::number_of_iterations(d));
  std::cerr << "Done (" << t.time() << " s)" << std::endl;

  std::ofstream out(out_file);
  out << pmesh;

  return 0;
}
