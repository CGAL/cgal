// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Eran Leiserowitz <leiserow@post.tau.ac.il>
//                 Efi Fogel    <efif@post.tau.ac.il>

#ifndef CGAL_ARR_SEGMENT_DEBUG_EXACT_TRAITS_H
#define CGAL_ARR_SEGMENT_DEBUG_EXACT_TRAITS_H


#include <CGAL/Pm_segment_debug_traits_2.h>
#include <CGAL/intersections.h>

#include <list>



CGAL_BEGIN_NAMESPACE

template <class Kernel_>
class My_arr_segment_traits_2 : public Pm_segment_traits_2<Kernel_>
{
public:

  class My_curve_2{
  public:
    typedef typename Kernel_::Point_2     Point_2;
    typedef typename Kernel::Segment_2   Segment;

    Segment seg;

    My_curve_2()
    {}

    My_curve_2(Segment s): seg(s)
    {}

    My_curve_2(Point_2 p1, Point_2 p2): seg(p1,p2)
    {}
  };



public:
  typedef Kernel_                               Kernel;
  typedef int                                   Info_face;
  typedef int                                   Info_edge;
  typedef int                                   Info_vertex;
  
  typedef Pm_segment_traits_2<Kernel>           Base;

  typedef typename Base::Has_left_category      Has_left_category;
  typedef typename Base::Has_reflect_category   Has_reflect_category;
    
  typedef typename Base::Point_2                Point_2;
  typedef typename Base::X_monotone_curve_2     X_monotone_curve_2;
  typedef My_curve_2                            Curve_2;

  // Obsolete, for backward compatibility
  typedef Point_2                               Point;
  typedef X_monotone_curve_2                    X_curve;
  typedef Curve_2                               Curve;

protected:
  typedef typename Kernel::Construct_vertex_2   Construct_vertex_2;
  typedef typename Kernel::Construct_segment_2  Construct_segment_2;
  typedef typename Kernel::Compare_x_2          Compare_x_2;
  typedef typename Kernel::Compare_y_2          Compare_y_2;
  typedef typename Kernel::Compare_xy_2         Compare_xy_2;
  typedef typename Kernel::Compare_y_at_x_2     Compare_y_at_x_2;
  typedef typename Kernel::Is_vertical_2        Is_vertical_2;
  typedef typename Kernel::Equal_2              Equal_2;
  typedef typename Kernel::Orientation_2        Orientation_2;
  
public:
  My_arr_segment_traits_2() : Base() { }

  /*! is_x_monotone()
   * \return true if the given curve is an x-monotone curve. False, otherwise.
   * For segments, this is always true
   */
//  bool is_x_monotone(const Curve_2 &) {return true;}
  
  /*! curve_make_x_monotone() cuts the given curve into x-monotone subcurves
   * and inserts them to the given output iterator. The order in which they
   * are inserted defines their order in the hierarchy tree.
   * While segments are x_monotone, still need to pass them out.
   * \param cv the input curve
   * \param o the output iterator
   * \return the past-the-end iterator
   */
  template<class OutputIterator>
  OutputIterator curve_make_x_monotone(const Curve_2 & cv,
                                       OutputIterator o) const
  {
    *o++ = cv.seg;
    return o;
  } 

  /*! curve_opposite() flips a given curve
   * \param cv the curve
   * \return a segment with source and target point interchanged
   */
  X_monotone_curve_2 curve_opposite(const X_monotone_curve_2 & cv) const
  {
    return construct_opposite_segment_2_object()(cv.seg);
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
  void curve_split(const X_monotone_curve_2 & cv,
                   X_monotone_curve_2 & c1, X_monotone_curve_2 & c2, 
                   const Point_2 & split_pt) const
  {
    //split curve at split point (x coordinate) into c1 and c2
    CGAL_precondition(curve_compare_y_at_x(split_pt, cv) == EQUAL);
    CGAL_precondition_code(Equal_2 is_equal = equal_2_object());
    CGAL_precondition(!is_equal(curve_source(cv), split_pt));
    CGAL_precondition(!is_equal(curve_target(cv), split_pt));
    
    Construct_vertex_2 construct_vertex = construct_vertex_2_object();
    const Point_2 & source = construct_vertex(cv.seg, 0);
    const Point_2 & target = construct_vertex(cv.seg, 1);
    Construct_segment_2 construct_segment = construct_segment_2_object();
    c1 = construct_segment(source, split_pt);
    c2 = construct_segment(split_pt, target);
  }

  /*! nearest_intersection_to_right() finds the nearest intersection point of
   * two given curves to the right of a given point. Nearest is defined as the
   * lexicographically nearest not including the point itself with one
   * exception explained bellow..
   * If the intersection of the two curves is an X_monotone_curve_2, that is,
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
  bool nearest_intersection_to_right(const X_monotone_curve_2 & c1,
                                     const X_monotone_curve_2 & c2,
                                     const Point_2 & pt,
                                     Point_2 & p1, Point_2 & p2) const
  {
    Object res = intersect_2_object()(c1.seg, c2.seg);

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
    X_monotone_curve_2 seg;
    if (assign(seg.seg, res)) {
      // the intersection is a curve:
      Construct_vertex_2 construct_vertex = construct_vertex_2_object();
      const Point_2 & src = construct_vertex(seg.seg, 0);
      const Point_2 & trg = construct_vertex(seg.seg, 1);
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

  /*! nearest_intersection_to_left() finds the nearest intersection point of
   * two given curves to the left of a given point. Nearest is defined as the
   * lexicographically nearest not including the point itself with one
   * exception explained bellow..
   * If the intersection of the two curves is an X_monotone_curve_2, that is,
   * there is an overlapping subcurve, then if the the source and target of the
   * subcurve are strickly to the left, they are returned through two
   * other point references p1 and p2. If pt is between the source and target
   * of the overlapping subcurve, or pt is its left endpoint, pt and the target
   * of the left endpoint of the subcurve are returned through p1 and p2 
   * respectively.
   * If the intersection of the two curves is a point to the left of pt, pt
   * is returned through the p1 and p2.
   * \param c1 the first curve
   * \param c2 the second curve
   * \param pt the point to compare against
   * \param p1 the first point reference
   * \param p2 the second point reference
   * \return true if c1 and c2 do intersect to the left of pt. Otherwise,
   * false
   */
  bool nearest_intersection_to_left(const X_monotone_curve_2 & c1,
                                    const X_monotone_curve_2 & c2,
                                    const Point_2 & pt,
                                    Point_2 & p1, Point_2 & p2) const
  {


    Object res = intersect_2_object()(c1.seg, c2.seg);

    // Empty object is returned - no intersection.
    if (res.is_empty()) return (false);

    // Intersection is a point
    if (assign(p1,res)) {
      // the intersection is a point:
      if (compare_xy_2_object()(p1, pt) == SMALLER) {
        p2 = p1;
        return true;
      }
      return false;
    }

    // Intersection is a segment
    bool same_curve=false;
    X_monotone_curve_2 seg;
    if (assign(seg.seg, res)) {
      if(same_curve)
	seg=c1.seg;
      // the intersection is a curve:
      Construct_vertex_2 construct_vertex = construct_vertex_2_object();
      const Point_2 & src = construct_vertex(seg.seg, 0);
      const Point_2 & trg = construct_vertex(seg.seg, 1);
      Compare_xy_2 compare_xy = compare_xy_2_object();
      Comparison_result src_pt = compare_xy(src, pt);
      Comparison_result trg_pt = compare_xy(trg, pt);

      if (src_pt == SMALLER && trg_pt == SMALLER) {
        // the subcurve is completely to the left:
        p1 = src;
        p2 = trg;
        return true;
      }
      
      if (trg_pt != SMALLER && src_pt == SMALLER) {
        // target is to the right, source is to the left:
        p1 = pt;
        p2 = src;
        return true;
      }

      if (src_pt != SMALLER && trg_pt == SMALLER) {
        // source is to the right, target is to the left:
        p1 = pt;
        p2 = trg;
        return true;
      }

      // the subcurve is completely to the right:
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
   */
  bool curves_overlap(const X_monotone_curve_2 & cv1,
                      const X_monotone_curve_2 & cv2) const
  {
    Construct_vertex_2 construct_vertex = construct_vertex_2_object();
    const Point_2 & src2 = construct_vertex(cv2.seg, 0);
    const Point_2 & trg2 = construct_vertex(cv2.seg, 1);
    const Point_2 & src1 = construct_vertex(cv1.seg, 0);
    const Point_2 & trg1 = construct_vertex(cv1.seg, 1);

    Orientation_2 orient = orientation_2_object();
    if ((!orient(src1, trg1, src2) == COLLINEAR) ||
        (!orient(src1, trg1, trg2) == COLLINEAR))
      return false;

    Is_vertical_2 is_vertical = is_vertical_2_object();
    if (is_vertical(cv1.seg)) {
      if (is_vertical(cv2.seg)) {
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
    if (is_vertical(cv2.seg)) return false;

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

  Point_2 curve_source(const Curve_2 & cv) const 
  { return construct_vertex_2_object()(cv.seg, 0); }
  
  Point_2 curve_target(const Curve_2 & cv) const 
  { return construct_vertex_2_object()(cv.seg, 1); }

  Point_2 curve_source(const X_monotone_curve_2 & cv) const 
  { return construct_vertex_2_object()(cv.seg, 0); }
  
  Point_2 curve_target(const X_monotone_curve_2 & cv) const 
  { return construct_vertex_2_object()(cv.seg, 1); }

};

CGAL_END_NAMESPACE

#endif
