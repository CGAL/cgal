//! \file examples/Boolean_set_operations_2/ex_intersection.C
// Computing the intersection of two hard coded polygons.

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/General_polygon_with_holes_2.h>

#include <list>

typedef CGAL::Gmpq                                      Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef Kernel::Point_2                                 Point_2;
typedef CGAL::Polygon_2<Kernel>                         Polygon_2;
typedef CGAL::General_polygon_with_holes_2<Polygon_2>   Polygon_with_holes_2;

int main(int argc, char * argv[])
{
  Polygon_2 p1, p2;
  p1.push_back(Point_2(0,0));     p1.push_back(Point_2(1.5,1.5));
  p1.push_back(Point_2(2.5,0.5)); p1.push_back(Point_2(3.5,1.5));
  p1.push_back(Point_2(5,0));
  p2.push_back(Point_2(0,2));     p2.push_back(Point_2(5,2));
  p2.push_back(Point_2(3.5,0.5)); p2.push_back(Point_2(2.5,1.5));
  p2.push_back(Point_2(1.5,0.5));
  std::list<Polygon_with_holes_2> result;
  CGAL::intersection(p1, p2, std::back_inserter(result));
  Polygon_with_holes_2 pwh;
  (void) CGAL::join(p1, p2, pwh);
  result.push_back(pwh);

  std::copy(result.begin(), result.end(),       // export to standard output
            std::ostream_iterator<Polygon_with_holes_2>(std::cout, "\n"));
  std::cout << std::endl;
  
  return 0;
}
