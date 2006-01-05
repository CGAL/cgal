//! \file examples/Boolean_set_operations_2/ex_traits_adaptor.C
// Using the traits adaptor to generate a traits of conics.

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <list>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits>
                                                        Arr_traits;
typedef CGAL::General_polygon_2<Arr_traits>             Polygon_2;
typedef CGAL::Gps_traits_2<Arr_traits,Polygon_2>        Traits;
typedef Traits::Polygon_with_holes_2                    Polygon_with_holes_2;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef Traits::Point_2                                 Point_2;

void conic_2_polygon(Curve_2 & conic, Polygon_2 & polygon)
{
  Traits traits;
  std::list<CGAL::Object> objects;
  traits.make_x_monotone_2_object()(conic, std::back_inserter(objects));
  std::list<CGAL::Object>::iterator i;
  for (i = objects.begin(); i != objects.end(); ++i) {
    X_monotone_curve_2 xcurve;
    CGAL::assign(xcurve, *i);
    polygon.push_back(xcurve);
  }
}

int main(int argc, char * argv[])
{
  // Construct a parabolic arc supported by a parabola: x^2 + 2y - 4 = 0),
  // and whose endpoints are on the line y = 0
  Curve_2 parabola1 =
    Curve_2(1, 0, 0, 0, 2, -4, CGAL::COUNTERCLOCKWISE,
            Point_2(2, 0), Point_2(-2, 0));
  // Construct a parabolic arc supported by a parabola: x^2 - 2y - 4 = 0),
  // and whose endpoints are on the line y = 0
  Curve_2 parabola2 =
    Curve_2(1, 0, 0, 0, -2, -4, CGAL::COUNTERCLOCKWISE,
            Point_2(-2, 0), Point_2(2, 0));
  
  Polygon_2 p1, p2;
  conic_2_polygon(parabola1, p1);
  conic_2_polygon(parabola2, p1);

  Curve_2 ellipse = Curve_2(1, 9, 0, 0, 0, -9);
  conic_2_polygon(ellipse, p2);
  std::list<Polygon_with_holes_2> result;
  CGAL::intersection(p1, p2, std::back_inserter(result));

  std::copy(result.begin(), result.end(),       // export to standard output
            std::ostream_iterator<Polygon_with_holes_2>(std::cout, "\n"));
  std::cout << std::endl;
  
  return 0;
}
