//-----------------------------------------------------------------------//
// This is the polygon example from the reference manual.
//-----------------------------------------------------------------------//

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <list>

typedef CGAL::Cartesian<double> R;
typedef CGAL::Polygon_traits_2<R> Traits;
typedef Traits::Point_2 Point;
typedef std::list<Point> Container;
typedef CGAL::Polygon_2<Traits,Container> Polygon;

#include <iostream>

int main()
{
  Polygon p;

  p.push_back(Point(0,0));
  p.push_back(Point(1,0));
  p.push_back(Point(1,1));
  p.push_back(Point(0,1));

  std::cout << "The polygon is " << 
    (p.is_convex() ? "" : "not ") << "convex." << std::endl;

  return 0;
}

