//! \file examples/Boolean_set_operations_2/example2.C
// Computing the intersection of two triangles.

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/General_polygon_with_holes_2.h>

#include <list>

typedef CGAL::Gmpq                                      Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef Kernel::Point_2                                 Point;
typedef CGAL::Polygon_2<Kernel>                         Polygon;
typedef CGAL::General_polygon_with_holes_2<Polygon>     Polygon_with_holes;

int main(int argc, char * argv[])
{
  Polygon p1, p2;
  p1.push_back(Point(0,0));     p1.push_back(Point(1.5,1.5));
  p1.push_back(Point(2.5,0.5)); p1.push_back(Point(3.5,1.5));
  p1.push_back(Point(5,0));
  p2.push_back(Point(0,2));     p2.push_back(Point(5,2));
  p2.push_back(Point(3.5,0.5)); p2.push_back(Point(2.5,1.5));
  p2.push_back(Point(1.5,0.5));
  std::list<Polygon_with_holes> result;
  CGAL::intersection(p1, p2, std::back_inserter(result));

  // Export to standard output:
  std::copy(result.begin(), result.end(),
            std::ostream_iterator<Polygon_with_holes>(std::cout, "\n"));
  std::cout << std::endl;
  
  return 0;
}
