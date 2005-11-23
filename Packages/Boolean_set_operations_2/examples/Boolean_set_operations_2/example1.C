
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>

typedef CGAL::Gmpq                                      Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef Kernel::Point_2                                 Point;
typedef CGAL::Polygon_2<Kernel>                         Polygon;

int main(int argc, char * argv[])
{
  Polygon p1, p2;
  p1.push_back(Point(-1,1));
  p1.push_back(Point(0,-1));
  p1.push_back(Point(1,1));
  p2.push_back(Point(-1,-1));
  p2.push_back(Point(1,-1));
  p2.push_back(Point(0,1));
  std::cout << (CGAL::do_intersect(p1, p2) ? "TRUE" : "FALSE") << std::endl;
  
  return 0;
}
