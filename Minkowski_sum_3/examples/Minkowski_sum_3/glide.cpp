#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/minkowski_sum_3.h>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>     Nef_polyhedron;
typedef Kernel::Point_3 Point_3;
typedef Point_3* point_iterator;
typedef std::pair<point_iterator,point_iterator> 
  point_range;
typedef std::list<point_range> polyline;

int main() 
{
  Nef_polyhedron N0;
  std::cin >> N0;
  Point_3 pl[6] = 
    {Point_3(-100,0,0), 
     Point_3(40,-70,0),
     Point_3(40,50,40),
     Point_3(-90,-60,60),
     Point_3(0, 0, -100),
     Point_3(30,0,150)
  };

  polyline poly;
  poly.push_back(point_range(pl,pl+6));
  Nef_polyhedron N1(poly.begin(), poly.end(), Nef_polyhedron::Polylines_tag());
  Nef_polyhedron result = CGAL::minkowski_sum_3(N0, N1);
}
