/*! \file traits_adapter.cpp
 * Using the traits adaptor to generate a traits of conics.
 */

#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl;
  return (0);
}
#else

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <list>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;

typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel,Nt_traits>
                                                        Conic_traits_2;

typedef CGAL::General_polygon_2<Conic_traits_2>         Polygon_2;
typedef CGAL::Gps_traits_2<Conic_traits_2, Polygon_2>   Traits_2;
typedef Traits_2::General_polygon_with_holes_2          Polygon_with_holes_2;
typedef Traits_2::Curve_2                               Curve_2;
typedef Traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef Traits_2::Point_2                               Point_2;

// Insert a conic arc as a polygon edge: Subdivide the arc into x-monotone
// sub-arcs and append these sub-arcs as polygon edges.
void append_conic_arc (Polygon_2& polygon, const Curve_2& arc)
{
  Conic_traits_2                    traits;
  std::list<CGAL::Object>           objects;
  std::list<CGAL::Object>::iterator it;
  X_monotone_curve_2                xarc;

  traits.make_x_monotone_2_object() (arc, std::back_inserter(objects));
  for (it = objects.begin(); it != objects.end(); ++it)
  {
    if (CGAL::assign (xarc, *it))
      polygon.push_back (xarc);
  }
}

int main ()
{
  // Construct a parabolic arc supported by a parabola: x^2 + 2y - 4 = 0,
  // and whose endpoints lie on the line y = 0:
  Curve_2 parabola1 = Curve_2 (1, 0, 0, 0, 2, -4, CGAL::COUNTERCLOCKWISE,
                               Point_2(2, 0), Point_2(-2, 0));

  // Construct a parabolic arc supported by a parabola: x^2 - 2y - 4 = 0,
  // and whose endpoints lie on the line y = 0:
  Curve_2 parabola2 = Curve_2 (1, 0, 0, 0, -2, -4, CGAL::COUNTERCLOCKWISE,
                               Point_2(-2, 0), Point_2(2, 0));

  // Construct a polygon from these two parabolic arcs.
  Polygon_2 P;
  append_conic_arc (P, parabola1);
  append_conic_arc (P, parabola2);

  // Construct a polygon that corresponds to the ellipse: x^2 + 9y^2 - 9 = 0:
  Polygon_2 Q;
  append_conic_arc (Q, Curve_2 (-1, -9, 0, 0, 0, 9));

  // Compute the intersection of the two polygons.
  std::list<Polygon_with_holes_2> res;
  CGAL::intersection (P, Q, std::back_inserter(res));

  std::copy (res.begin(), res.end(),       // export to standard output
             std::ostream_iterator<Polygon_with_holes_2>(std::cout, "\n"));
  std::cout << std::endl;

  return (0);
}

#endif
