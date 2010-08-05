#include <CGAL/Simple_cartesian.h>
#include <CGAL/Subdivision_method_3.h>

#include <cstdio>
#include <iostream>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;

using namespace std;
using namespace CGAL;

// ======================================================================
template <class Poly>
class WLoop_mask_3 {
  typedef Poly                                         Polyhedron;

  typedef typename Polyhedron::Vertex_iterator         Vertex_iterator;
  typedef typename Polyhedron::Halfedge_iterator       Halfedge_iterator;
  typedef typename Polyhedron::Facet_iterator          Facet_iterator;

  typedef typename Polyhedron::Halfedge_around_facet_circulator
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator
                                            Halfedge_around_vertex_circulator;

  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;

  typedef typename Kernel::FT                          FT;
  typedef typename Kernel::Point_3                     Point;
  typedef typename Kernel::Vector_3                    Vector;

public:
  void edge_node(Halfedge_iterator eitr, Point& pt) {
    Point& p1 = eitr->vertex()->point();
    Point& p2 = eitr->opposite()->vertex()->point();
    Point& f1 = eitr->next()->vertex()->point();
    Point& f2 = eitr->opposite()->next()->vertex()->point();

    pt = Point((3*(p1[0]+p2[0])+f1[0]+f2[0])/8,
	       (3*(p1[1]+p2[1])+f1[1]+f2[1])/8,
	       (3*(p1[2]+p2[2])+f1[2]+f2[2])/8 );
  }
  void vertex_node(Vertex_iterator vitr, Point& pt) {
    double R[] = {0.0, 0.0, 0.0};
    Point& S = vitr->point();

    Halfedge_around_vertex_circulator vcir = vitr->vertex_begin();
    std::size_t n = circulator_size(vcir);
    for (std::size_t i = 0; i < n; i++, ++vcir) {
      Point& p = vcir->opposite()->vertex()->point();
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
  void border_node(Halfedge_iterator eitr, Point& ept, Point& vpt) {
    Point& ep1 = eitr->vertex()->point();
    Point& ep2 = eitr->opposite()->vertex()->point();
    ept = Point((ep1[0]+ep2[0])/2, (ep1[1]+ep2[1])/2, (ep1[2]+ep2[2])/2);

    Halfedge_around_vertex_circulator vcir = eitr->vertex_begin();
    Point& vp1  = vcir->opposite()->vertex()->point();
    Point& vp0  = vcir->vertex()->point();
    Point& vp_1 = (--vcir)->opposite()->vertex()->point();
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

  Subdivision_method_3::PTQ(P, WLoop_mask_3<Polyhedron>(), d);

  cout << P; // write the .off

  return 0;
}
