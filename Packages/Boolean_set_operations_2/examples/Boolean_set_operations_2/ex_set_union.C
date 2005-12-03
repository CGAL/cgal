//! \file examples/Boolean_set_operations_2/ex_set_union.C
// Computing the union of a set of circles.

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <list>
#include <stdlib.h>
#include <math.h>

typedef CGAL::Gmpq                                      NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Gps_circle_segment_traits_2<Kernel>       Traits;
typedef Traits::Polygon_2                               Polygon;
typedef Traits::Polygon_with_holes_2                    Polygon_with_holes;
typedef Traits::Curve_2                                 Curve;
typedef Traits::X_monotone_curve_2                      X_monotone_curve;
typedef Kernel::Point_2                                 Point;
typedef Kernel::Circle_2                                Circle;

void circle_2_polygon(Circle circle, Polygon & polygon)
{
  Traits traits;
  Curve curve(circle);
  std::list<CGAL::Object> objects;
  traits.make_x_monotone_2_object()(curve, std::back_inserter(objects));
  std::list<CGAL::Object>::iterator i = objects.begin();
  X_monotone_curve xcurve;
  CGAL::assign(xcurve, *i++);   polygon.push_back(xcurve);
  CGAL::assign(xcurve, *i);     polygon.push_back(xcurve);
}

int main(int argc, char * argv[])
{
  const double pi = atan(1) * 4;

  unsigned int circles_num = 8;
  if (argc > 1) circles_num = atoi(argv[1]);

  double circles_num_reciep = 1 / circles_num;
  double radius = 1;
  double frac = 2 * pi * circles_num_reciep;
  std::list<Polygon> polygons;
  for (unsigned int i = 0; i < circles_num; ++i) {
    double angle = frac * i;
    double x = radius * sin(angle);
    double y = radius * cos(angle);
    Point center = Point(x, y);
    Circle circle(center, radius);
    Polygon polygon;
    circle_2_polygon(circle, polygon);
    polygons.push_back(polygon);
  }
  
  std::list<Polygon_with_holes> result;
  join(polygons.begin(), polygons.end(), std::back_inserter(result));

  std::copy(result.begin(), result.end(),       // export to standard output
            std::ostream_iterator<Polygon_with_holes>(std::cout, "\n"));
  std::cout << std::endl;
  
  return 0;
}
