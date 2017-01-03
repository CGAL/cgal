#include <CGAL/Simple_cartesian.h>
#include <CGAL/Subdivision_method_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <boost/foreach.hpp>

#include <iostream>

typedef CGAL::Simple_cartesian<double>     Kernel;
//typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef CGAL::Surface_mesh<Kernel::Point_3> Polyhedron;
using namespace std;
using namespace CGAL;

// ======================================================================
template <class Poly>
class WLoop_mask_3 {
  typedef Poly                                         Polyhedron;

  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;

  typedef Halfedge_around_target_circulator<Poly> Halfedge_around_target_circulator;

  typedef typename boost::property_map<Polyhedron, vertex_point_t>::type Vertex_pmap;
  typedef typename boost::property_traits<Vertex_pmap>::value_type Point;
  
  Polyhedron& polyhedron;
  Vertex_pmap vpm;

public:
  WLoop_mask_3(Polyhedron& polyhedron)
    : polyhedron(polyhedron), vpm(get(CGAL::vertex_point, polyhedron))
  {}

  void edge_node(halfedge_descriptor hd, Point& pt) {
    Point& p1 = get(vpm, target(hd,polyhedron));
    Point& p2 = get(vpm, target(opposite(hd,polyhedron),polyhedron));
    Point& f1 = get(vpm, target(next(hd,polyhedron),polyhedron));
    Point& f2 = get(vpm, target(next(opposite(hd,polyhedron),polyhedron),polyhedron));

    pt = Point((3*(p1[0]+p2[0])+f1[0]+f2[0])/8,
	       (3*(p1[1]+p2[1])+f1[1]+f2[1])/8,
	       (3*(p1[2]+p2[2])+f1[2]+f2[2])/8 );
  }
  void vertex_node(vertex_descriptor vd, Point& pt) {
    double R[] = {0.0, 0.0, 0.0};
    Point& S = get(vpm,vd);

    std::size_t n = 0;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(vd, polyhedron)){
      ++n;
      Point& p = get(vpm, target(opposite(hd,polyhedron),polyhedron));
      R[0] += p[0]; 	R[1] += p[1]; 	R[2] += p[2];
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
    Point& ep1 = get(vpm, target(hd,polyhedron));
    Point& ep2 = get(vpm, target(opposite(hd,polyhedron),polyhedron));
    ept = Point((ep1[0]+ep2[0])/2, (ep1[1]+ep2[1])/2, (ep1[2]+ep2[2])/2);

    Halfedge_around_target_circulator vcir(hd,polyhedron);
    Point& vp1  = get(vpm, target(opposite(*vcir,polyhedron),polyhedron));
    Point& vp0  = get(vpm, target(*vcir,polyhedron));
    --vcir;
    Point& vp_1 = get(vpm,target(opposite(*vcir,polyhedron),polyhedron));
    vpt = Point((vp_1[0] + 6*vp0[0] + vp1[0])/8,
                (vp_1[1] + 6*vp0[1] + vp1[1])/8,
                (vp_1[2] + 6*vp0[2] + vp1[2])/8 );
  }
};

int main(int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: Customized_subdivision d < filename" << endl;
    cout << "       d: the depth of the subdivision (0 < d < 10)" << endl;
    cout << "       filename: the input mesh (.off)" << endl;
    return 0;
  }

  int d = argv[1][0] - '0';

  Polyhedron P;
  cin >> P; // read the .off

  Subdivision_method_3::PTQ(P, WLoop_mask_3<Polyhedron>(P), d);

  cout << P; // write the .off

  return 0;
}
