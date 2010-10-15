/*! \file set_union.cpp
 * Computing the union of a set of circles.
 */

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Lazy_exact_nt.h>

#include <list>
#include <cstdlib>
#include <cmath>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point_2;
typedef Kernel::Circle_2                                  Circle_2;
typedef CGAL::Gps_circle_segment_traits_2<Kernel>         Traits_2;

typedef CGAL::General_polygon_set_2<Traits_2>             Polygon_set_2;
typedef Traits_2::Polygon_2                               Polygon_2;
typedef Traits_2::Polygon_with_holes_2                    Polygon_with_holes_2;
typedef Traits_2::Curve_2                                 Curve_2;
typedef Traits_2::X_monotone_curve_2                      X_monotone_curve_2;

// Construct a polygon from a circle.
Polygon_2 construct_polygon (const Circle_2& circle)
{
  // Subdivide the circle into two x-monotone arcs.
  Traits_2 traits;
  Curve_2 curve (circle);
  std::list<CGAL::Object> objects;
  traits.make_x_monotone_2_object() (curve, std::back_inserter(objects));
  CGAL_assertion (objects.size() == 2);

  // Construct the polygon.
  Polygon_2 pgn;
  X_monotone_curve_2 arc;
  std::list<CGAL::Object>::iterator iter;

  for (iter = objects.begin(); iter != objects.end(); ++iter) {
    CGAL::assign (arc, *iter);
    pgn.push_back (arc);
  }

  return pgn;
}

// The main program:
int main (int argc, char* argv[])
{
  // Read the number of circles from the command line.
  unsigned int n_circles = 8;
  if (argc > 1) n_circles = std::atoi(argv[1]);

  // Create the circles, equally spaced of the circle x^2 + y^2 = 1.
  const double pi = std::atan(1.0) * 4;
  const double n_circles_reciep = 1.0 / n_circles;
  const double radius = 1;
  const double frac = 2 * pi * n_circles_reciep;
  std::list<Polygon_2> circles;
  unsigned int k;
  for (k = 0; k < n_circles; k++) {
    const double angle = frac * k;
    const double x = radius * std::sin(angle);
    const double y = radius * std::cos(angle);
    Point_2 center = Point_2(x, y);
    Circle_2 circle(center, radius);

    circles.push_back (construct_polygon (circle));
  }

  // Compute the union aggragately.
  std::list<Polygon_with_holes_2> res;
  CGAL::join (circles.begin(), circles.end(), std::back_inserter (res));

  // Print the result.
  std::copy (res.begin(), res.end(),
             std::ostream_iterator<Polygon_with_holes_2>(std::cout, "\n"));
  std::cout << std::endl;

  return 0;
}
