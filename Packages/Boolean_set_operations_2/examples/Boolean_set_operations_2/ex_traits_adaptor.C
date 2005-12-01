//! \file examples/Boolean_set_operations_2/ex_traits_adaptor.C
// Using the traits adaptor to generate a traits of conics.

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/Gps_traits_adaptor_2.h>

#include <list>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef Rat_kernel::Point_2                             Rat_point_2;
typedef Rat_kernel::Segment_2                           Rat_segment_2;
typedef Rat_kernel::Circle_2                            Rat_circle_2;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits>
                                                        Arr_traits;
typedef CGAL::General_polygon_2<Arr_traits>             Polygon;
typedef CGAL::Gps_traits_adaptor_2<Arr_traits,Polygon>  Traits;
typedef CGAL::General_polygon_set_2<Traits>             Polygon_set;
typedef Traits::Polygon_with_holes_2                    Polygon_with_holes;
typedef Traits::Curve_2                                 Curve;
typedef Traits::X_monotone_curve_2                      X_monotone_curve;
typedef Traits::Point_2                                 Point;

void conic_2_polygon(Curve & conic, Polygon & polygon)
{
  Traits traits;
  std::list<CGAL::Object> objects;
  traits.make_x_monotone_2_object()(conic, std::back_inserter(objects));
  std::list<CGAL::Object>::iterator i;
  for (i = objects.begin(); i != objects.begin(); ++i) {
    X_monotone_curve xcurve;
    CGAL::assign(xcurve, *i);
    polygon.push_back(xcurve);
  }
}

int main(int argc, char * argv[])
{
  // Construct a parabolic arc supported by a parabola: x^2 + y - 2 = 0),
  // and whose endpoints are on the line y = 0
  Curve parabola1 =
    Curve(1, 0, 0, 0, 1, -2, CGAL::COUNTERCLOCKWISE,
          Point(1.41, 0),       // Approximation of the source.
          0, 0, 0, 0, 1, 0,     // The line: y = 0.
          Point(-1.41, 0),      // Approximation of the target.
          0, 0, 0, 0, 1, 0);
  // Construct a parabolic arc supported by a parabola: x^2 - y - 2 = 0),
  // and whose endpoints are on the line y = 0
  Curve parabola2 =
    Curve(1, 0, 0, 0, -1, -2, CGAL::COUNTERCLOCKWISE,
          Point(-1.41, 0),      // Approximation of the source.
          0, 0, 0, 0, 1, 0,     // The line: y = 0.
          Point(1.41, 0),       // Approximation of the target.
          0, 0, 0, 0, 1, 0);
  
  Polygon p1, p2;
  conic_2_polygon(parabola1, p1);
  conic_2_polygon(parabola2, p1);

  Curve ellipse = Curve(1, 4, 0, 0, 0, -4);
  conic_2_polygon(ellipse, p2);
  std::list<Polygon_with_holes> result;
  CGAL::intersection(p1, p2, std::back_inserter(result));

  std::copy(result.begin(), result.end(),       // export to standard output
            std::ostream_iterator<Polygon_with_holes>(std::cout, "\n"));
  std::cout << std::endl;
  
  return 0;
}
