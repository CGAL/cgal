// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file : include/CGAL/Circular_arc_traits_checker.h

#ifndef CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_CHECKER_H
#define CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_CHECKER_H

// TODO :
// - Similar to Pm_checker_traits, but it doesn't require conversion
//   between the 2 types (especially endpoints).
// - Maybe this should move to TAU's CGAL packages (?)
// - It's not specific to circular arcs...

#include <CGAL/basic.h>
#include <cassert>
#include <utility>

namespace CGAL {

/// Traits checker : compares the output of all predicates by running two
/// traits in parallel.

template < typename Traits_1, typename Traits_2 >
class Circular_arc_traits_checker {

  Traits_1 t1;
  Traits_2 t2;

public:

  typedef Traits_1                                     First_traits;
  typedef Traits_2                                     Second_traits;

  typedef std::pair<typename Traits_1::Curve,
                    typename Traits_2::Curve>          Curve;
  typedef std::pair<typename Traits_1::Point,
                    typename Traits_2::Point>          Point;

  typedef Curve                                        Curve_2;
  typedef Curve                                        X_curve;
  typedef Curve                                        X_monotone_curve_2;
  typedef Point                                        Point_2;

  typedef typename Traits_1::Has_left_category         Has_left_category;

  Circular_arc_traits_checker()
  {
    // "Check" that the categories are the same...
    typename Traits_1::Has_left_category cat1;
    typename Traits_2::Has_left_category cat2 = cat1;
  }

  CGAL::Comparison_result
  compare_x(const Point_2 &p, const Point_2 &q) const
  {
    CGAL::Comparison_result ret1 = t1.compare_x(p.first, q.first);
    CGAL::Comparison_result ret2 = t2.compare_x(p.second, q.second);
    assert( ret1 == ret2);
    return ret1;
  }

  CGAL::Comparison_result
  compare_xy(const Point_2 &p, const Point_2 &q) const
  {
    CGAL::Comparison_result ret1 = t1.compare_xy(p.first, q.first);
    CGAL::Comparison_result ret2 = t2.compare_xy(p.second, q.second);
    assert( ret1 == ret2);
    return ret1;
  }

  bool
  point_equal(const Point_2 &p, const Point_2 &q) const
  {
    bool ret1 = t1.point_equal(p.first, q.first);
    bool ret2 = t2.point_equal(p.second, q.second);
    assert( ret1 == ret2);
    return ret1;
  }

  bool
  curve_equal(const X_monotone_curve_2 & cv1,
              const X_monotone_curve_2 & cv2) const
  {
    bool ret1 = t1.curve_equal(cv1.first, cv2.first);
    bool ret2 = t2.curve_equal(cv1.second, cv2.second);
    assert( ret1 == ret2);
    return ret1;
   }

  Point_2
  curve_source(const X_monotone_curve_2 &cv) const
  {
    // [Approximate] comparison for constructions ?
    return Point_2(t1.curve_source(cv.first),
                   t2.curve_source(cv.second));
  }

  Point_2
  curve_target(const X_monotone_curve_2 &cv) const
  {
    return Point_2(t1.curve_target(cv.first),
                   t2.curve_target(cv.second));
  }

  bool
  curve_is_vertical(const X_monotone_curve_2 &cv) const
  {
    bool ret1 = t1.curve_is_vertical(cv.first);
    bool ret2 = t2.curve_is_vertical(cv.second);
    assert( ret1 == ret2);
    return ret1;
  }

  bool
  point_in_x_range(const X_monotone_curve_2 & cv,
                   const Point_2 & p) const
  {
    bool ret1 = t1.point_in_x_range(cv.first, p.first);
    bool ret2 = t2.point_in_x_range(cv.second, p.second);
    assert( ret1 == ret2);
    return ret1;
  }

  CGAL::Comparison_result
  curve_compare_y_at_x(const Point_2 & p,
                       const X_monotone_curve_2 & cv) const
  {
    CGAL::Comparison_result ret1 = t1.curve_compare_y_at_x(p.first, cv.first);
    CGAL::Comparison_result ret2 = t2.curve_compare_y_at_x(p.second, cv.second);
    assert( ret1 == ret2);
    return ret1;
  }

  CGAL::Comparison_result
  curves_compare_y_at_x(const X_monotone_curve_2 & cv1,
                        const X_monotone_curve_2 & cv2,
                        const Point_2 & p) const
  {
    CGAL::Comparison_result ret1, ret2;
    ret1 = t1.curves_compare_y_at_x(cv1.first, cv2.first, p.first);
    ret2 = t2.curves_compare_y_at_x(cv1.second, cv2.second, p.second);
    assert( ret1 == ret2);
    return ret1;
  }

  CGAL::Comparison_result
  curves_compare_y_at_x_right(const X_monotone_curve_2 & cv1,
                              const X_monotone_curve_2 & cv2,
                              const Point_2 & p) const
  {
    CGAL::Comparison_result ret1, ret2;
    ret1 = t1.curves_compare_y_at_x_right(cv1.first, cv2.first, p.first);
    ret2 = t2.curves_compare_y_at_x_right(cv1.second, cv2.second, p.second);
    assert( ret1 == ret2);
    return ret1;
  }

  Point_2
  point_reflect_in_x_and_y(const Point_2 & p) const
  {
    return Point_2(t1.point_reflect_in_x_and_y(p.first),
                   t2.point_reflect_in_x_and_y(p.second));
  }

  X_monotone_curve_2
  curve_reflect_in_x_and_y(const X_monotone_curve_2 & cv) const
  {
    return X_monotone_curve_2(t1.curve_reflect_in_x_and_y(cv.first),
                              t2.curve_reflect_in_x_and_y(cv.second));
  }

  template < typename OutputIterator >
  OutputIterator
  curve_make_x_monotone(const X_monotone_curve_2 &cv, OutputIterator res) const
  {
    std::vector<typename X_monotone_curve_2::first_type>  V1;
    std::vector<typename X_monotone_curve_2::second_type> V2;

    t1.curve_make_x_monotone(cv.first,  std::back_inserter(V1));
    t2.curve_make_x_monotone(cv.second, std::back_inserter(V2));

    assert(V1.size() == V2.size());

    for (unsigned i = 0; i != V1.size(); ++i)
      *res++ = X_monotone_curve_2(V1[i], V2[i]);

    return res;
  }

  void
  curve_split(const X_monotone_curve_2 & cv,
              X_monotone_curve_2 & cv1,
              X_monotone_curve_2 & cv2,
              const Point_2 & p) const
  {
    t1.curve_split(cv.first, cv1.first, cv2.first, p.first);
    t2.curve_split(cv.second, cv1.second, cv2.second, p.second);
  }

  bool
  curves_overlap(const X_monotone_curve_2 &cv1,
                 const X_monotone_curve_2 &cv2) const
  {
    bool ret1 = t1.curves_overlap(cv1.first, cv2.first);
    bool ret2 = t2.curves_overlap(cv1.second, cv2.second);
    assert( ret1 == ret2);
    return ret1;
  }

  bool
  nearest_intersection_to_right(const X_monotone_curve_2 & cv1,
                                const X_monotone_curve_2 & cv2,
                                const Point_2 & pt,
                                Point_2 & p1,
                                Point_2 & p2) const
  {
    bool ret1, ret2;
    ret1 = t1.nearest_intersection_to_right(cv1.first, cv2.first,
                                            pt.first, p1.first, p2.first);
    ret2 = t2.nearest_intersection_to_right(cv1.second, cv2.second,
                                            pt.second, p1.second, p2.second);
    assert( ret1 == ret2);
    return ret1;
  }

  bool
  is_x_monotone(const Curve_2 &cv) const
  {
    bool ret1 = t1.is_x_monotone(cv.first);
    bool ret2 = t2.is_x_monotone(cv.second);
    assert( ret1 == ret2);
    return ret1;
  }

  X_monotone_curve_2
  curve_opposite(const X_monotone_curve_2 & cv) const
  {
    return X_monotone_curve_2(t1.curve_opposite(cv.first),
                              t2.curve_opposite(cv.second));
  }

};

} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_CHECKER_H
