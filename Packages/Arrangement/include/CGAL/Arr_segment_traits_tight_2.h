// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/Arr_segment_traits_tight_2.h
// package       : arr (1.03)
// source        :
// revision      :
// revision_date :
// author(s)     : Iddo Hanniel <hanniel@math.tau.ac.il>
//                 Efi Fogel    <efif@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin)
// chapter       : Arrangement_2
//
// ======================================================================
#ifndef CGAL_ARR_SEGMENT_TRAITS_TIGHT_H
#define CGAL_ARR_SEGMENT_TRAITS_TIGHT_H

#include <CGAL/tags.h>
#include <list>

CGAL_BEGIN_NAMESPACE

template <class Kernel_>
class Arr_segment_traits_2
{
public:
  typedef Kernel_                               Kernel;
  typedef int                                   Info_face;
  typedef int                                   Info_edge;
  typedef int                                   Info_vertex;
  
  // Categories:
  // #define HAS_LEFT_NOT
#if !defined(HAS_LEFT_NOT)
  typedef Tag_true                              Has_left_category;
#else
  typedef Tag_false                             Has_left_category;
#endif
    
  // Traits objects
  typedef typename Kernel::Point_2              Point_2;
  typedef typename Kernel::Segment_2            X_curve_2;

  typedef X_curve_2                             Curve_2;

  // Obsolete, for backward compatibility
  typedef Point_2                               Point;
  typedef X_curve_2                             X_curve;
  typedef Curve_2                               Curve;

protected:
  typedef typename Kernel::Is_vertical_2        Is_vertical_2;
  typedef typename Kernel::Construct_vertex_2   Construct_vertex_2;
  typedef typename Kernel::Less_x_2             Less_x_2;
  typedef typename Kernel::Equal_2              Equal_2;
  typedef typename Kernel::Construct_vertex_2   Construct_vertex_2;
  typedef typename Kernel::Construct_segment_2  Construct_segment_2;
  typedef typename Kernel::Compare_x_2          Compare_x_2;
  typedef typename Kernel::Compare_y_2          Compare_y_2;
  typedef typename Kernel::Compare_xy_2         Compare_xy_2;
  typedef typename Kernel::Compare_y_at_x_2     Compare_y_at_x_2;
  typedef typename Kernel::Is_vertical_2        Is_vertical_2;
  typedef typename Kernel::Orientation_2        Orientation_2;
  
public:
  Arr_segment_traits_2() {}

   * \param p1 the first point
   * \param p2 the second point
   * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2); 
   *         SMALLER if x(p1) < x(p2), or if x(p1) = x(p2) and y(p1) < y(p2);
   *         or else EQUAL
   */
  Comparison_result compare_xy(const Point_2 & p1, const Point_2 & p2) const
  { return compare_xy_2_object()(p1, p2); }

  /*! curve_is_vertical()
   * \param cv the curve
   * \return true iff the curve is vertical
   *
   * \todo replace indirect use curve_is_vertical() with is_vertical_2()
   */
  bool curve_is_vertical(const X_curve_2 & cv) const 
  { return is_vertical_2_object()(cv); }

  /*! curve_is_in_x_range()
   * \param cv the curve
   * \param q the point
   * \return true if q is in the x range of cv
   *
   * \todo Intorduce Is_in_x_range_2() or perhaps Is_in_x_closed_range_2()
   * in kernel. Currently, this is implemented using existing traits (kernel)
   * functions (curve_source(), curve_target()) that return the source and
   * target points by value, which is not as efficient as possible.
   */
  bool curve_is_in_x_range(const X_curve_2 & cv, const Point_2 & q) const
  {
    Construct_vertex_2 construct_vertex = construct_vertex_2_object();
    const Point_2 & source = construct_vertex(cv, 0);
    const Point_2 & target = construct_vertex(cv, 1);
    Less_x_2 less_x = less_x_2_object();
    return !((less_x(source, q) && less_x(target, q)) ||
             (less_x(q, source) && less_x(q, target)));
  }

  /*! curve_compare_at_x() compares the y-coordinate of two given curves at
   * the x-coordinate of a given point.
   * \param cv1 the first curve
   * \param cv2 the second curve
   * \param q the point
   * \return LARGER if cv1(x(q)) > cv2(x(q)); SMALLER if cv1(x(q)) < cv2(x(q));
   * or else EQUAL.
   * \pre The point q is in the x range of the two curves.
   *
   * \todo replace indirect use curve_compare_at_x() with compare_y_at_x_2()
   */
  Comparison_result curve_compare_at_x(const X_curve_2 & cv1, 
				       const X_curve_2 & cv2, 
				       const Point_2 & q) const
  {
    CGAL_precondition(curve_is_in_x_range(cv1, q));
    CGAL_precondition(curve_is_in_x_range(cv2, q));

    return compare_y_at_x_2_object()(q, cv1, cv2);
  }

#if !defined(HAS_LEFT_NOT)
  /*! curve_compare_at_x_left() compares the y value of two curves in an
   * epsilon environment to the left of the x value of the input point
   * Preconditions: The point q is in the x range of the two curves, and both
   * of them must be also be defined to its left. The two curves must also
   * intersect at x(q).
   */
  Comparison_result curve_compare_at_x_left(const X_curve_2 & cv1,
                                            const X_curve_2 & cv2, 
                                            const Point_2 & q) const 
  {
    // The two curves must not be vertical.
    CGAL_precondition(! curve_is_vertical(cv1));
    CGAL_precondition(! curve_is_vertical(cv2));

    // The two curve must be defined at q and also to its left.
    CGAL_precondition_code(
        Construct_vertex_2 construct_vertex = construct_vertex_2_object();
	Less_x_2 less_x = less_x_2_object();
	const Point_2 & source1 = construct_vertex(cv1, 0);
	const Point_2 & target1 = construct_vertex(cv1, 1);
	const Point_2 & source2 = construct_vertex(cv2, 0);
	const Point_2 & target2 = construct_vertex(cv2, 1);
	);

    CGAL_precondition (less_x(source1, q) || less_x(target1, q));
    CGAL_precondition (!(less_x(source1, q) && less_x(target1, q)));
    
    CGAL_precondition (less_x(source2, q) || less_x(target2, q));
    CGAL_precondition (!(less_x(source2, q) && less_x(target2, q)));
    
    // Since the curves are continuous, if they are not equal at q, the same
    // result also applies to q's left.
    CGAL_precondition (compare_y_at_x_2_object()(q, cv1, cv2) == EQUAL);
    
    // <cv2> and <cv1> meet at a point with the same x-coordinate as q
    // compare their derivatives.
    return compare_slope_2_object()(cv2, cv1);
  }
#else
  /*! point_reflect_in_x_and_y() reflects the given point about the origin
   */
  Point_2 point_reflect_in_x_and_y(const Point_2 & pt) const
  {
    Point_2 org = construct_point_2_object()(ORIGIN);      
    typename Kernel::Vector_2 v = construct_vector_2_object()(pt, org);
    Point_2 reflected_pt(v);
    return reflected_pt;
  }

  /*! curve_reflect_in_x_and_y reflects the given curve about the origin
   */
  X_curve_2 curve_reflect_in_x_and_y(const X_curve_2 & cv) const
  {
    X_curve_2 reflected_cv(point_reflect_in_x_and_y ( cv.source()),
                           point_reflect_in_x_and_y ( cv.target()));
    return reflected_cv;
  }
#endif
    
  /*! curve_compare_at_x_right() compares the y value of two curves in an
   * epsilon environment to the right of the x value of the input point
   * Preconditions: The point q is in the x range of the two curves, and both
   * of them must be also be defined to its right. The two curves must also
   * intersect at x(q).
   */
  Comparison_result curve_compare_at_x_right(const X_curve_2 & cv1,
                                             const X_curve_2 & cv2, 
                                             const Point_2 & q) const 
  {
    // The two curves must not be vertical.
    CGAL_precondition(! curve_is_vertical(cv1));
    CGAL_precondition(! curve_is_vertical(cv2));

    // The two curve must be defined at q and also to its right.
    CGAL_precondition_code(
        Construct_vertex_2 construct_vertex = construct_vertex_2_object();
	Less_x_2 less_x = less_x_2_object();
	const Point_2 & source1 = construct_vertex(cv1, 0);
	const Point_2 & target1 = construct_vertex(cv1, 1);
	const Point_2 & source2 = construct_vertex(cv2, 0);
	const Point_2 & target2 = construct_vertex(cv2, 1);
	);

    CGAL_precondition (less_x(q, source1) || less_x(q, target1));
    CGAL_precondition (!(less_x(q, source1) && less_x(q, target1)));
    
    CGAL_precondition (less_x(q, source2) || less_x(q, target2));
    CGAL_precondition (!(less_x(q, source2) && less_x(q, target2)));
    
    // Since the curves are continuous, if they are not equal at q, the same
    // result also applies to q's left.
    CGAL_precondition (curve_compare_at_x(cv1, cv2, q) == EQUAL);     
    
    // <cv1> and <cv2> meet at a point with the same x-coordinate as q
    // compare their derivatives
    return compare_slope_2_object()(cv1, cv2);
  }
    
  /*! Return the curve-point status of the input objects.
   * \pre p must be in the x-range of cv.
   */
  Comparison_result curve_get_point_status (const X_curve_2 & cv, 
					    const Point_2 & p) const
  {
    CGAL_precondition(curve_is_in_x_range(cv, p));

    Comparison_result res = compare_y_at_x_2_object()(p, cv);
    if (res == SMALLER)
      return (LARGER);
    else if (res == LARGER)
      return (SMALLER);
    return (EQUAL);
  }

  /*! \todo replace indirect use curve_is_same() with equal_2()
   */
  bool curve_is_same(const X_curve_2 & cv1,const X_curve_2 & cv2) const
  {
    Equal_2 equal = equal_2_object();
    const X_curve_2 & ocv1 = construct_opposite_segment_2_object()(cv1);
    return equal(cv1, cv2) || equal(ocv1, cv2);
  }

  /*! \todo replace indirect use point_is_same() with equal_2()
   */
  bool point_is_same(const Point_2 & p1, const Point_2 & p2) const
  {
    return equal_2_object()(p1, p2);
  }
  
  /*! \todo replace indirect use curve_source() with construct_vertex_2()
   */
  Point_2 curve_source(const X_curve_2 & cv) const 
  { return construct_vertex_2_object()(cv, 0); }

  /*! \todo replace indirect use curve_target() with construct_vertex_2()
   */
  Point_2 curve_target(const X_curve_2 & cv) const 
  { return construct_vertex_2_object()(cv, 1); }

  ///////////////////////////////////////////////////////////////////////////
  // Planar map with intersection requirement:
  ///////////////////////////////////////////////////////////////////////////

  /*! is_x_monotone()
   * \return true if the given curve is an x-monotone curve. False, otherwise.
   * For segments, this is always true
   */
  bool is_x_monotone(const Curve_2 &) {return true;}
  
  /*! make_x_monotone() cuts the given curve into x-monotone subcurves and
   * stores them in the given list. The order in which they are inserted into
   * the list defines their order in the hierarchy tree.
   * Segments are x_monotone!
   */
  void make_x_monotone(const Curve_2 & /*cv*/, std::list<Curve_2>& /*l*/) {} 

  /*! curve_flip() flips a given curve
   * \param cv the curve
   * \return a segment with source and target point interchanged
   * \todo replace indirect use curve_flip() with opposite()
   */
  X_curve_2 curve_flip(const X_curve_2 & cv) const
  {
    return construct_opposite_segment_2_object()(cv);
  }
 
  /*! curve_split() splits a given curve at a given split point into two
   * sub-curves.
   * \param cv the curve to split
   * \param c1 the output first part of the split curve.Its source is the
   * source of the original curve.
   * \param c2 the output second part of the split curve. Its target is the
   * target of the original curve
   * \param split_pt
   * \pre split_pt is on cv but is not an endpoint.
   */
  void curve_split(const X_curve_2 & cv, X_curve_2 & c1, X_curve_2 & c2, 
                   const Point_2 & split_pt)
  {
    //split curve at split point (x coordinate) into c1 and c2
    CGAL_precondition(curve_get_point_status(cv, split_pt) == EQUAL);
    CGAL_precondition_code(Compare_xy_2 compare_xy = compare_xy_2_object());
    CGAL_precondition(compare_xy(curve_source(cv), split_pt) != EQUAL);
    CGAL_precondition(compare_xy(curve_target(cv), split_pt) != EQUAL);
    
    Construct_vertex_2 construct_vertex = construct_vertex_2_object();
    const Point_2 & source = construct_vertex(cv, 0);
    const Point_2 & target = construct_vertex(cv, 1);
    Construct_segment_2 construct_segment = construct_segment_2_object();
    c1 = construct_segment(source, split_pt);
    c2 = construct_segment(split_pt, target);
  }

  /*! nearest_intersection_to_right() finds the nearest intersection point of
   * two given curves to the right of a given point. Nearest is defined as the
   * lexicographically nearest not including the point itself with one
   * exception explained bellow..
   * If the intersection of the two curves is an X_curve_2, that is,
   * there is an overlapping subcurve, then if the the source and target of the
   * subcurve are strickly to the right, they are returned through two
   * other point references p1 and p2. If pt is between the source and target
   * of the overlapping subcurve, or pt is its left endpoint, pt and the target
   * of the right endpoint of the subcurve are returned through p1 and p2 
   * respectively.
   * If the intersection of the two curves is a point to the right of pt, pt
   * is returned through the p1 and p2.
   * \param c1 the first curve
   * \param c2 the second curve
   * \param pt the point to compare against
   * \param p1 the first point reference
   * \param p2 the second point reference
   * \return true if c1 and c2 do intersect to the right of pt. Otherwise,
   * false
   */
  bool nearest_intersection_to_right(const X_curve_2 & c1,
                                     const X_curve_2 & c2,
                                     const Point_2 & pt,
                                     Point_2 & p1, Point_2 & p2) const
  {
    Object res = intersect_2_object()(c1, c2);

    // Empty object is returned - no intersection.
    if (res.is_empty()) return (false);

    // Intersection is a point
    Point_2 ip;
    if (assign(p1,res)) {
      // the intersection is a point:
      if (compare_xy_2_object()(p1, pt) == LARGER) {
        p2 = p1;
        return true;
      }
      return false;
    }
    
    // Intersection is a segment
    X_curve_2 seg;
    if (assign(seg, res)) {
      // the intersection is a curve:
      Construct_vertex_2 construct_vertex = construct_vertex_2_object();
      const Point_2 & src = construct_vertex(seg, 0);
      const Point_2 & trg = construct_vertex(seg, 1);
      Compare_xy_2 compare_xy = compare_xy_2_object();
      Comparison_result src_pt = compare_xy(src, pt);
      Comparison_result trg_pt = compare_xy(trg, pt);

      if (src_pt == LARGER && trg_pt == LARGER) {
        // the subcurve is completely to the right:
        p1 = src;
        p2 = trg;
        return true;
      }
      
      if (trg_pt != LARGER && src_pt == LARGER) {
        // target is to the left, source is to the right:
        p1 = pt;
        p2 = src;
        return true;
      }

      if (src_pt != LARGER && trg_pt == LARGER) {
        // source is to the left, target is to the right:
        p1 = pt;
        p2 = trg;
        return true;
      }

      // the subcurve is completely to the left:
      return false;
    }

    // the curves do not intersect:
    return false;
  }

  /*! curves_overlap() test overlapping between two given curves
   * \patam c1 the first curve
   * \patam c2 the second curve
   * \return true if c1 and c2 overlap in a one-dimensional subcurve
   * (i.e., not in a finite number of points). Otherwise, false.
   * \todo end point coincidence instead of intersection!
   */
  bool curves_overlap(const X_curve_2 & cv1, const X_curve_2 & cv2) const
  {
    Construct_vertex_2 construct_vertex = construct_vertex_2_object();
    const Point_2 & src2 = construct_vertex(cv2, 0);
    const Point_2 & trg2 = construct_vertex(cv2, 1);
    const Point_2 & src1 = construct_vertex(cv1, 0);
    const Point_2 & trg1 = construct_vertex(cv1, 1);

    Orientation_2 orient = orientation_2_object();
    if ((!orient(src1, trg1, src2) == COLLINEAR) ||
        (!orient(src1, trg1, trg2) == COLLINEAR))
      return false;

    Is_vertical_2 is_vertical = is_vertical_2_object();
    if (is_vertical(cv1)) {
      if (is_vertical(cv2)) {
        Compare_y_2 compare_y = compare_y_2_object();
        Comparison_result res_ss = compare_y (src1, src2);
        Comparison_result res_st = compare_y (src1, trg2);
        if (res_ss == SMALLER) {
          if (res_st == LARGER) return true;
          if (compare_y (trg1, src2) == LARGER) return true;
          return (compare_y (trg1, trg2) == LARGER);
        }

        if (res_ss == LARGER) {
          if (res_st == SMALLER) return true;
          if (compare_y (trg1, src2) == SMALLER) return true;
          return (compare_y (trg1, trg2) == SMALLER);
        }

        // res_ss == EQUAL
        if (res_st == SMALLER)
          return (compare_y (trg1, src2) == LARGER);
        return (compare_y (trg1, src2) == SMALLER);
      }
      return false;
    }
    if (is_vertical(cv2)) return false;

    Compare_x_2 compare_x = compare_x_2_object();
    Comparison_result res_ss = compare_x (src1, src2);
    Comparison_result res_st = compare_x (src1, trg2);
    if (res_ss == SMALLER) {
      if (res_st == LARGER) return true;
      if (compare_x (trg1, src2) == LARGER) return true;
      return (compare_x (trg1, trg2) == LARGER);
    }

    if (res_ss == LARGER) {
      if (res_st == SMALLER) return true;
      if (compare_x (trg1, src2) == SMALLER) return true;
      return (compare_x (trg1, trg2) == SMALLER);
    }

    // res_ss == EQUAL
    if (res_st == SMALLER)
      return (compare_x (trg1, src2) == LARGER);
    return (compare_x (trg1, src2) == SMALLER);
  }

};

CGAL_END_NAMESPACE

#endif
