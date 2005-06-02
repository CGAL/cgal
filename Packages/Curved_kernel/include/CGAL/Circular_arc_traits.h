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

// file : include/CGAL/Circular_arc_traits.h

#ifndef CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_H
#define CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_H

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/global_functions_on_circular_arcs_2.h>

namespace CGAL {

/// Traits class for CGAL::Arrangement_2 (and similar) based on a CurvedKernel.

template < typename CurvedKernel >
class Circular_arc_traits {

  CurvedKernel ck;

public:

  typedef CurvedKernel Kernel;
  typedef typename CurvedKernel::Circular_arc_2  Curve_2;
  typedef typename CurvedKernel::Circular_arc_2  X_monotone_curve_2;

  typedef typename CurvedKernel::Circular_arc_endpoint_2      Point;
  typedef typename CurvedKernel::Circular_arc_endpoint_2      Point_2;

  typedef CGAL::Tag_false                        Has_left_category;
  typedef CGAL::Tag_false 			 Has_merge_category;

  Circular_arc_traits(const CurvedKernel &k = CurvedKernel())
    : ck(k) {}

  typedef typename CurvedKernel::Compare_x_2           Compare_x_2;
  typedef typename CurvedKernel::Compare_xy_2          Compare_xy_2;
  typedef typename CurvedKernel::Compare_y_at_x_2      Compare_y_at_x_2;
  typedef typename CurvedKernel::Compare_y_to_right_2  Compare_y_at_x_right_2; 
  typedef typename CurvedKernel::Equal_2               Equal_2;
  typedef typename CurvedKernel::Make_x_monotone_2     Make_x_monotone_2;
  typedef typename CurvedKernel::Split_2               Split_2;
  typedef typename CurvedKernel::Construct_intersections_2 Intersect_2;

  class Construct_min_vertex_2
  {
  public:
    const Point_2& operator() (const X_monotone_curve_2 & cv) const
    {
      assert( CurvedKernel().compare_xy_2_object()(cv.left(),cv.right())==CGAL::SMALLER);
      return (cv.left());
    }
  };

  class Construct_max_vertex_2
  {
  public:
    /*!
     * Get the right endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The right endpoint.
     */
    const Point_2& operator() (const X_monotone_curve_2 & cv) const
    {
      assert( CurvedKernel().compare_xy_2_object()(cv.left(),cv.right())==CGAL::SMALLER);
      return (cv.right());
    }
  };

  class Is_vertical_2
  {
  public:
    typedef bool result_type;

    // TO BE IMPLEMENTED !!!!!!!
    bool operator() (const X_monotone_curve_2& cv) const
    {
      return false;
    }
  };

  
  Compare_x_2 compare_x_2_object() const
  { return ck.compare_x_2_object(); }

  Compare_xy_2 compare_xy_2_object() const
  { return ck.compare_xy_2_object(); }

  Compare_y_at_x_2 compare_y_at_x_2_object() const 
  { return ck.compare_y_at_x_2_object(); }

  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const 
  { return ck.compare_y_to_right_2_object(); }

  Equal_2 equal_2_object() const
  { return ck.equal_2_object(); }

  Make_x_monotone_2 make_x_monotone_2_object() const
  { return ck.make_x_monotone_2_object(); }

  Split_2 split_2_object() const
  { return ck.split_2_object(); }

  Intersect_2 intersect_2_object() const
    { return ck.construct_intersections_2_object(); }

  Construct_min_vertex_2 construct_min_vertex_2_object() const
    { return Construct_min_vertex_2(); }

  Construct_max_vertex_2 construct_max_vertex_2_object() const
    { return Construct_max_vertex_2(); }

  Is_vertical_2 is_vertical_2_object() const
    { return Is_vertical_2();}

  //////////////////////////////////////////////
  CGAL::Comparison_result
  compare_x(const Point_2 &p, const Point_2 &q) const
  {
    return ck.compare_x_2_object()(p, q);
  }

  CGAL::Comparison_result
  compare_xy(const Point_2 &p, const Point_2 &q) const
  {
    return ck.compare_xy_2_object()(p, q);
  }

  bool
  point_equal(const Point_2 &p, const Point_2 &q) const
  {
    return ck.equal_2_object()(p, q);
  }

  bool
  curve_equal(const X_monotone_curve_2 & cv1,
              const X_monotone_curve_2 & cv2) const
  {
    return ck.equal_2_object()(cv1, cv2);
  }

  const Point_2 &
  curve_source(const X_monotone_curve_2 &cv) const
  {
    return cv.source();
  }

  const Point_2 &
  curve_target(const X_monotone_curve_2 &cv) const
  {
    return cv.target();
  }

  bool
  curve_is_vertical(const X_monotone_curve_2 &/*cv*/) const
  {
    return false;
  }

  bool
  point_in_x_range(const X_monotone_curve_2 & cv,
                   const Point_2 & p) const
  {
    return ck.in_range_2_object()(cv, p);
  }

  CGAL::Comparison_result
  curve_compare_y_at_x(const Point_2 & p,
                       const X_monotone_curve_2 & cv) const
  {
    return ck.compare_y_at_x_2_object()(p, cv);
  }

  // That's one of those painful requirements...
  CGAL::Comparison_result
  curves_compare_y_at_x(const X_monotone_curve_2 & cv1,
                        const X_monotone_curve_2 & cv2,
                        const Point_2 & p) const
  {
    assert(point_in_range(cv1, p));
    assert(point_in_range(cv2, p));

    // The point is probably on one of the two curves for the sweep.
    // Otherwise we don't know how to implement it exactly...
    // At least for the sweep that should be true.

    // We know how to compare a point and a curve, so let's try to see
    // if this is enough to determine the result.
    CGAL::Comparison_result cmp1 = curve_compare_y_at_x(p, cv1);
    CGAL::Comparison_result cmp2 = curve_compare_y_at_x(p, cv2);

    if (cmp1 != cmp2 || cmp1 == CGAL::EQUAL)
      // We can determine the result exactly in this case.
      return CGAL::compare(cmp2, cmp1);

    // What we cannot determine is the case where cv1 and cv2 are on
    // the same side of p (either strictly above, or strictly below).
    std::cout << __FUNCTION__ << " : case of approximate evaluation"
              << std::endl;
    return CGAL::compare(cv1.approximate_y_at(p),
                            cv2.approximate_y_at(p));
  }

  CGAL::Comparison_result
  curves_compare_y_at_x_right(const X_monotone_curve_2 & cv1,
                              const X_monotone_curve_2 & cv2,
                              const Point_2 & p) const
  {
    return ck.compare_y_to_right_2_object()(cv1, cv2, p);
  }

private:
  // For internal use in the 2 predicates that follow.
  typename CurvedKernel::Circle_2
  opposite_circle(const typename CurvedKernel::Circle_2 &c) const
  {
    typedef typename CurvedKernel::Circle_2                 Circle_2;
    typedef typename CurvedKernel::Linear_kernel::Point_2   Pt_2;
    Pt_2 new_center (-c.center().x(), -c.center().y());
    return Circle_2(new_center, c.squared_radius());
  }

public:

  // Used by the incremental building (for point location ?)
  Point_2
  point_reflect_in_x_and_y(const Point_2 & p) const
  {
    assert(CGAL::is_valid(p.x()) && CGAL::is_valid(p.y()));
    return Point_2 (opposite_circle(p.circle(0)),
                    opposite_circle(p.circle(1)),
                    ! p.is_left());
  }

  // Used by the incremental building (for point location ?)
  X_monotone_curve_2
  curve_reflect_in_x_and_y(const X_monotone_curve_2 & cv) const
  {
    X_monotone_curve_2 ret (opposite_circle(cv.supporting_circle()),
                            opposite_circle(cv.begin().circle(1)),
                                          ! cv.begin().is_left(),
                            opposite_circle(cv.end().circle(1)),
                                          ! cv.end().is_left());
    if (cv.swap_source())
      ret.swap_source_and_target();

    return ret;
  }

  template < typename OutputIterator >
  OutputIterator
  curve_make_x_monotone(const Curve_2 &cv, OutputIterator res) const
  {
    return ck.make_x_monotone_2_object()(cv, res);
  }

  void
  curve_split(const X_monotone_curve_2 & cv,
              X_monotone_curve_2 & cv1,
              X_monotone_curve_2 & cv2,
              const Point_2 & p) const
  {
    typename CurvedKernel::Split_2::circle_result_type
             pair = ck.split_2_object()(cv, p);
    cv1 = pair.first;
    cv2 = pair.second;
  }

  bool
  curves_overlap(const X_monotone_curve_2 &cv1,
                 const X_monotone_curve_2 &cv2) const
  {
    return ck.do_overlap_2_object()(cv1, cv2);
  }

  bool
  nearest_intersection_to_right(const X_monotone_curve_2 & cv1,
                                const X_monotone_curve_2 & cv2,
                                const Point_2 & pt,
                                Point_2 & p1,
                                Point_2 & p2) const
  {
    return ck.nearest_intersection_to_right_2_object()(cv1, cv2, pt, p1, p2);
  }

  // Undocumented requirement.
  // This one and curve_opposite() will hopefully go to the trash soon.
  bool
  is_x_monotone(const Curve_2 & cv) const
  {
    return cv.is_x_monotone();
  }

  X_monotone_curve_2
  curve_opposite(const X_monotone_curve_2 & cv) const
  {
    X_monotone_curve_2 ret = cv;
    ret.swap_source_and_target();
    return ret;
  }

};

} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_CIRCULAR_ARC_TRAITS_H
