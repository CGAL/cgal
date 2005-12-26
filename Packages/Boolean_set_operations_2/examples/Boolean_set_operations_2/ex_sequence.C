//! \file examples/Boolean_set_operations_2/ex_sequence.C
// Performing a sequence of Boolean Set-Operations.

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/General_polygon_with_holes_2.h>

#include <list>

typedef CGAL::Gmpq                                      Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef Kernel::Point_2                                 Point_2;
typedef CGAL::Polygon_2<Kernel>                         Polygon_2;
typedef CGAL::General_polygon_with_holes_2<Polygon_2>   Polygon_with_holes_2;
typedef CGAL::Gps_segment_traits_2<Kernel>              Traits;
typedef CGAL::General_polygon_set_2<Traits>             Polygon_set_2;

int main(int argc, char * argv[])
{
  Polygon_2 p1, p2, rect;
  p1.push_back(Point_2(0,1));
  p1.push_back(Point_2(2,0));
  p1.push_back(Point_2(1,1));
  p1.push_back(Point_2(2,2));

  p2.push_back(Point_2(3,1));
  p2.push_back(Point_2(1,2));
  p2.push_back(Point_2(2,1));
  p2.push_back(Point_2(1,0));

  rect.push_back(Point_2(0,0));
  rect.push_back(Point_2(3,0));
  rect.push_back(Point_2(3,2));
  rect.push_back(Point_2(0,2));
  
  Polygon_set_2 ps;
  ps.insert(p1);
  ps.join(p2);                  // Compute the union
  ps.complement();              // Compute the complement
  ps.intersection(rect);        // Compute the intersection
  
  std::list<Polygon_with_holes_2> result;
  ps.polygons_with_holes(std::back_inserter(result));

  std::copy(result.begin(), result.end(),       // export to standard output
            std::ostream_iterator<Polygon_with_holes_2>(std::cout, "\n"));
  std::cout << std::endl;
  
  return 0;
}
