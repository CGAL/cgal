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

// file : include/CGAL/Circular_arc_traits_tracer.h

#ifndef CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_TRACER_H
#define CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_TRACER_H

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Profile_counter.h>

namespace CGAL {

/// Wrapper around a traits class for Arrangement which just forwards the
/// calls to its nested traits, and prints debug info for each argument and
/// return value.  Compiling with -DCGAL_PROFILE dumps the number of calls
/// to each predicate at the end of the program.

// to get better profiles for the traits interface :
#ifndef CGAL_NO_INLINE
#  define CGAL_NO_INLINE __attribute__((__noinline__))
//#define CGAL_NO_INLINE
#endif

// FIXME : should it be the responsibility of the tracer
//         to add the Debug_id stuff ?
template < typename Traits, int debug_level = 2 >
class Circular_arc_traits_tracer {

// debug_level :
// 0 = no debug messages
// 1 = important messages
// 2 = maybe important messages
// 3 = all calls

  Traits traits;

  mutable int nested_level; /// for detection of nested calls
  mutable int call_number;  /// identifies the call number

// FIXME : transform these macros into member functions ?
#define CGAL_DEBUG(level, args) if (level <= debug_level) \
   { std::cout << "Call #" << call_number++ << " " << __FUNCTION__ \
               << "  ( " << args << " )"; \
     if (nested_level != 0) \
       std::cout << "   [ nested level " << nested_level << " ]"; \
     ++nested_level; \
     std::cout << std::endl; }

#define CGAL_DEBUG_RET(level, args) if (level <= debug_level) \
   { std::cout << "         returns : " << args; \
     --nested_level; \
     if (nested_level != 0) \
       std::cout << "   [ nested level " << nested_level << " ]"; \
     std::cout << std::endl; }

public:

  Circular_arc_traits_tracer(const Traits &t = Traits())
    : traits(t), nested_level(0), call_number(0) {}

  typedef typename Traits::Curve                 Curve;
  typedef typename Traits::Curve_2               Curve_2;
  typedef typename Traits::X_curve               X_curve;
  typedef typename Traits::X_monotone_curve_2    X_monotone_curve_2;

  typedef typename Traits::Point                 Point;
  typedef typename Traits::Point_2               Point_2;

  typedef typename Traits::Has_left_category     Has_left_category;

  CGAL::Comparison_result
  CGAL_NO_INLINE
  compare_x(const Point_2 &p, const Point_2 &q) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(2, p.id() << ", " << q.id());
    CGAL::Comparison_result ret = traits.compare_x(p, q);
    CGAL_DEBUG_RET(2, ret);
    return ret;
  }

  CGAL::Comparison_result
  CGAL_NO_INLINE
  compare_xy(const Point_2 &p, const Point_2 &q) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(2, p.id() << ", " << q.id());
    CGAL::Comparison_result ret = traits.compare_xy(p, q);
    CGAL_DEBUG_RET(2, ret);
    return ret;
  }

  bool
  CGAL_NO_INLINE
  point_equal(const Point_2 &p, const Point_2 &q) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(3, p.id() << ", " << q.id());
    bool ret = traits.point_equal(p, q);
    CGAL_DEBUG_RET(3, ret);
    return ret;
  }

  bool
  CGAL_NO_INLINE
  curve_equal(const X_monotone_curve_2 & cv1,
              const X_monotone_curve_2 & cv2) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(1, cv1.id() << ", " << cv2.id());
    bool ret = traits.curve_equal(cv1, cv2);
    CGAL_DEBUG_RET(1, ret);
    return ret;
  }

  const Point_2 &
  CGAL_NO_INLINE
  curve_source(const X_monotone_curve_2 &cv) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(3, cv.id());
    const Point_2 & ret = traits.curve_source(cv);
    CGAL_DEBUG_RET(3, ret.id());
    return ret;
  }

  const Point_2 &
  CGAL_NO_INLINE
  curve_target(const X_monotone_curve_2 &cv) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(3, cv.id());
    const Point_2 & ret = cv.target();
    CGAL_DEBUG_RET(3, ret.id());
    return ret;
  }

  bool
  CGAL_NO_INLINE
  curve_is_vertical(const X_monotone_curve_2 &cv) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(3, cv.id());
    bool ret = traits.curve_is_vertical(cv);
    CGAL_DEBUG_RET(3, ret);
    return ret;
  }

  bool
  CGAL_NO_INLINE
  point_in_x_range(const X_monotone_curve_2 & cv,
                   const Point_2 & p) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(2, cv.id() << ", " << p.id());
    bool ret = traits.point_in_x_range(cv, p);
    CGAL_DEBUG_RET(2, ret);
    return ret;
  }

  CGAL::Comparison_result
  CGAL_NO_INLINE
  curve_compare_y_at_x(const Point_2 & p,
                       const X_monotone_curve_2 & cv) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(1, p.id() << ", " << cv.id());
    CGAL::Comparison_result ret = traits.curve_compare_y_at_x(p, cv);
    CGAL_DEBUG_RET(1, ret);
    return ret;
  }

  CGAL::Comparison_result
  CGAL_NO_INLINE
  curves_compare_y_at_x(const X_monotone_curve_2 & cv1,
                        const X_monotone_curve_2 & cv2,
                        const Point_2 & p) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(1, cv1.id() << ", " << cv2.id() << ", " << p.id());
    CGAL::Comparison_result ret = traits.curves_compare_y_at_x(cv1, cv2, p);
    CGAL_DEBUG_RET(1, ret);
    return ret;
  }

  CGAL::Comparison_result
  CGAL_NO_INLINE
  curves_compare_y_at_x_right(const X_monotone_curve_2 & cv1,
                              const X_monotone_curve_2 & cv2,
                              const Point_2 & p) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(1, cv1.id() << ", " << cv2.id() << ", " << p.id());
    CGAL::Comparison_result ret;
    ret = traits.curves_compare_y_at_x_right(cv1, cv2, p);
    CGAL_DEBUG_RET(1, ret);
    return ret;
  }

  Point_2
  CGAL_NO_INLINE
  point_reflect_in_x_and_y(const Point_2 & p) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(1, p.id());
    Point_2 ret = traits.point_reflect_in_x_and_y(p);
    CGAL_DEBUG_RET(1, ret.id());
    return ret;
  }

  X_monotone_curve_2
  CGAL_NO_INLINE
  curve_reflect_in_x_and_y(const X_monotone_curve_2 & cv) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(1, cv.id());
    X_monotone_curve_2 ret = traits.curve_reflect_in_x_and_y(cv);
    CGAL_DEBUG_RET(1, ret.id());
    return ret;
  }

  template < typename OutputIterator >
  OutputIterator
  CGAL_NO_INLINE
  curve_make_x_monotone(const Curve_2 &cv, OutputIterator res) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(1, cv.id());
    CGAL_DEBUG_RET(1, "output iterator...");
    return traits.curve_make_x_monotone(cv, res);
  }

  void
  CGAL_NO_INLINE
  curve_split(const X_monotone_curve_2 & cv,
              X_monotone_curve_2 & cv1,
              X_monotone_curve_2 & cv2,
              const Point_2 & p) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(1, cv.id() << ", " << p.id());
    traits.curve_split(cv, cv1, cv2, p);
    CGAL_DEBUG_RET(1, cv1.id() << ", " << cv2.id());
  }

  bool
  CGAL_NO_INLINE
  curves_overlap(const X_monotone_curve_2 &cv1,
                 const X_monotone_curve_2 &cv2) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(1, cv1.id() << ", " << cv2.id());
    bool ret = traits.curves_overlap(cv1, cv2);
    CGAL_DEBUG_RET(1, ret);
    return ret;
  }

  bool
  CGAL_NO_INLINE
  nearest_intersection_to_right(const X_monotone_curve_2 & cv1,
                                const X_monotone_curve_2 & cv2,
                                const Point_2 & pt,
                                Point_2 & p1,
                                Point_2 & p2) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(1, cv1.id() << ", " << cv2.id() << ", " << pt.id());
    bool ret= traits.nearest_intersection_to_right(cv1, cv2, pt, p1, p2);
    CGAL_DEBUG_RET(1, ret << ", " << p1.id() << ", " << p2.id());
    return ret;
  }

  bool
  CGAL_NO_INLINE
  is_x_monotone(const Curve_2 & cv) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(2, cv.id());
    bool ret = traits.is_x_monotone(cv);
    CGAL_DEBUG_RET(2, ret);
    return ret;
  }

  X_monotone_curve_2
  CGAL_NO_INLINE
  curve_opposite(const X_monotone_curve_2 & cv) const
  {
    CGAL_PROFILER(__FUNCTION__);
    CGAL_DEBUG(2, cv.id());
    X_monotone_curve_2 ret = traits.curve_opposite(cv);
    CGAL_DEBUG_RET(2, ret.id());
    return cv2;
  }

};

} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_TRACER_H
