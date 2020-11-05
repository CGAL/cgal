#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

typedef CGAL::Exact_integer NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Point_3* point_iterator;
typedef std::pair<point_iterator,point_iterator> point_range;
typedef std::list<point_range> polyline;

int main() {

  Point_3 pl[4] = {Point_3(0,0,0), Point_3(1,0,0),
                   Point_3(1,1,0), Point_3(0,1,0)};
  polyline poly;
  poly.push_back(point_range(pl,pl+4));
  Nef_polyhedron N(poly.begin(), poly.end(), Nef_polyhedron::Polylines_tag());
  std::cout << N;
}
