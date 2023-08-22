#include <CGAL/Exact_rational.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Bounded_kernel.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Lazy_exact_nt<CGAL::Exact_rational> FT;
typedef CGAL::Simple_cartesian<FT> Kernel;
typedef CGAL::Bounded_kernel<Kernel> Bounded_kernel;
typedef CGAL::Nef_polyhedron_2<Bounded_kernel> Nef_polyhedron;
typedef Nef_polyhedron::Point Point;


int main(int argc, char *argv[])
{
  Point p1[] = {
    Point(0,0),
    Point(100,100),
    Point(0,100)
  };

  Point p2[] = {
    Point(100, 100),
    Point(110, 100),
    Point(100, 110)
  };


  std::list<std::pair<Point*, Point*> > polygons1;
  polygons1.push_back(std::make_pair(p1, p1 + sizeof(p1) / sizeof(Point)));

  Nef_polyhedron poly1(polygons1.begin(), polygons1.end(), Nef_polyhedron::POLYGONS);
  poly1.explorer().check_integrity_and_topological_planarity();

  Nef_polyhedron poly2(p2, p2 + sizeof(p2) / sizeof(Point));
  poly2.explorer().check_integrity_and_topological_planarity();

  Nef_polyhedron intersect = poly1.intersection(poly2);
  intersect.explorer().check_integrity_and_topological_planarity();

  Nef_polyhedron comp = intersect.complement();  //leads to crash/exception since topological plane map of intersect is not correct!
  comp.explorer().check_integrity_and_topological_planarity();
  return 0;
}
