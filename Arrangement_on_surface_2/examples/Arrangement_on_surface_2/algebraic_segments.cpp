//! \file examples/Arrangement_on_surface_2/algebraic_segments.cpp
// Constructing an arrangement of algebraic segments.

#include <iostream>
#include <CGAL/config.h>

#if (!CGAL_USE_CORE) && (!CGAL_USE_LEDA) && (!(CGAL_USE_GMP && CGAL_USE_MPFI))
int main ()
{
  std::cout << "Sorry, this example needs CORE, LEDA, or GMP+MPFI ..."
            << std::endl;
  return 0;
}
#else

#include <CGAL/basic.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>

#include "integer_type.h"
#include "arr_print.h"

typedef CGAL::Arr_algebraic_segment_traits_2<Integer> Traits;
typedef CGAL::Arrangement_2<Traits>                   Arrangement;
typedef Traits::Curve_2                               Curve;
typedef Traits::Polynomial_2                          Polynomial;
typedef Traits::Algebraic_real_1                      Algebraic_real;
typedef Traits::X_monotone_curve_2                    X_monotone_curve;
typedef Traits::Point_2                               Point;

typedef boost::variant<Point, X_monotone_curve>       Make_x_monotone_result;

int main() {
  Traits traits;

  auto make_xmon = traits.make_x_monotone_2_object();
  auto ctr_cv = traits.construct_curve_2_object();
  auto ctr_pt = traits.construct_point_2_object();
  auto construct_xseg = traits.construct_x_monotone_segment_2_object();

  Polynomial x = CGAL::shift(Polynomial(1), 1, 0);
  Polynomial y = CGAL::shift(Polynomial(1), 1, 1);

  // Construct a curve (C1) with the equation x^4+y^3-1=0.
  Curve cv1 = ctr_cv(CGAL::ipower(x, 4) + CGAL::ipower(y, 3) - 1);
  // Construct all x-monotone segments using the Make_x_mononotone functor.
  std::vector<Make_x_monotone_result> pre_segs;
  make_xmon(cv1, std::back_inserter(pre_segs));
  // Cast all CGAL::Objects into X_monotone_segment
  // (the vector might also contain Point objects for isolated points,
  // but not in this case).
  std::vector<X_monotone_curve> segs;
  for(size_t i = 0; i < pre_segs.size(); ++i) {
    auto* curr_p = boost::get<X_monotone_curve>(&pre_segs[i]);
    CGAL_assertion(curr_p);
    segs.push_back(*curr_p);
  }
  // Construct an ellipse (C2) with the equation 2*x^2+5*y^2-7=0.
  Curve cv2 = ctr_cv(2*CGAL::ipower(x,2)+5*CGAL::ipower(y,2)-7);

  // Construct point on the upper arc (counting of arc numbers starts with 0).
  Point p11 = ctr_pt(Algebraic_real(0), cv2, 1);

  construct_xseg(cv2, p11, Traits::POINT_IN_INTERIOR,
                 std::back_inserter(segs));

  // Construct a vertical cusp (C3) with the equation x^2-y^3=0.
  Curve cv3 = ctr_cv(CGAL::ipower(x, 2)-CGAL::ipower(y, 3));

  // Construct a segment containing the cusp point.
  // This adds two X_monotone_curve objects to the vector,
  // because the cusp is a critical point.
  Point p21 = ctr_pt(Algebraic_real(-2), cv3, 0);
  Point p22 = ctr_pt(Algebraic_real(2), cv3, 0);
  construct_xseg(cv3 ,p21, p22, std::back_inserter(segs));

  // Construct an unbounded curve, starting at x=3.
  Point p23 = ctr_pt(Algebraic_real(3), cv3, 0);
  construct_xseg(cv3, p23, Traits::MIN_ENDPOINT, std::back_inserter(segs));

  // Construct another conic (C4) with the equation y^2-x^2+1=0.
  Curve cv4 = ctr_cv(CGAL::ipower(y,2)-CGAL::ipower(x,2)+1);

  Point p31 = ctr_pt(Algebraic_real(2), cv4, 1);
  construct_xseg(cv4, p31, Traits::MAX_ENDPOINT, std::back_inserter(segs));

  // Construct a vertical segment (C5).
  Curve cv5 = ctr_cv(x);
  Point v1 = ctr_pt(Algebraic_real(0), cv3, 0);
  Point v2 = ctr_pt(Algebraic_real(0), cv2, 1);
  construct_xseg(cv5, v1, v2, std::back_inserter(segs));

  Arrangement arr(&traits);
  CGAL::insert(arr, segs.begin(), segs.end());

  // Add some isolated points (must be wrapped into CGAL::Object).
  std::vector<CGAL::Object> isolated_pts;
  isolated_pts.push_back(CGAL::make_object(ctr_pt(Algebraic_real(2), cv4, 0)));
  isolated_pts.push_back(CGAL::make_object(ctr_pt(Integer(1), Integer(2))));
  isolated_pts.push_back(CGAL::make_object(ctr_pt(Algebraic_real(-1),
                                                  Algebraic_real(2))));
  CGAL::insert(arr, isolated_pts.begin(), isolated_pts.end());

  print_arrangement_size(arr);                  // print the arrangement size

  return 0;
}

#endif
