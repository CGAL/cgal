//! \file examples/Boolean_set_operations_2/ex_set_union.C
// Computing the union of a set of circles.

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <list>
#include <stdlib.h>
#include <cmath>

typedef CGAL::Gmpq                                      NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Gps_circle_segment_traits_2<Kernel>       Traits;
typedef Traits::Polygon_2                               Polygon_2;
typedef Traits::Polygon_with_holes_2                    Polygon_with_holes_2;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef Kernel::Point_2                                 Point_2;
typedef Kernel::Circle_2                                Circle_2;

void circle_2_polygon(Circle_2 circle, Polygon_2 & polygon)
{
  Traits traits;
  Curve_2 curve(circle);
  std::list<CGAL::Object> objects;
  traits.make_x_monotone_2_object()(curve, std::back_inserter(objects));
  std::list<CGAL::Object>::iterator i = objects.begin();
  X_monotone_curve_2 xcurve;
  CGAL::assign(xcurve, *i++);   polygon.push_back(xcurve);
  CGAL::assign(xcurve, *i);     polygon.push_back(xcurve);
}

int main(int argc, char * argv[])
{
  const double pi = std::atan(1.0) * 4;

  unsigned int circles_num = 8;
  if (argc > 1) circles_num = atoi(argv[1]);

  double circles_num_reciep = 1.0 / circles_num;
  double radius = 1;
  double frac = 2 * pi * circles_num_reciep;
  std::list<Polygon_2> polygons;
  for (unsigned int i = 0; i < circles_num; ++i) {
    double angle = frac * i;
    double x = radius * std::sin(angle);
    double y = radius * std::cos(angle);
    Point_2 center = Point_2(x, y);
    Circle_2 circle(center, radius);
    Polygon_2 polygon;
    circle_2_polygon(circle, polygon);
    polygons.push_back(polygon);
  }
  
  std::list<Polygon_with_holes_2> result;
  join(polygons.begin(), polygons.end(), std::back_inserter(result));

  std::copy(result.begin(), result.end(),       // export to standard output
            std::ostream_iterator<Polygon_with_holes_2>(std::cout, "\n"));
  std::cout << std::endl;
  
  return 0;
}
