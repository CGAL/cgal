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
// file          : include/CGAL/Arr_segment_traits_2.h
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
#ifndef CGAL_ARR_SEGMENT_EXACT_TRAITS_H
#define CGAL_ARR_SEGMENT_EXACT_TRAITS_H

#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/intersections.h>

#include <list>

CGAL_BEGIN_NAMESPACE

template <class Kernel_>
class Arr_segment_traits_2 : public Pm_segment_traits_2<Kernel_>
{
public:
  typedef Kernel_                               Kernel;
  typedef int                                   Info_face;
  typedef int                                   Info_edge;
  typedef int                                   Info_vertex;
  
  typedef Pm_segment_traits_2<Kernel>           Base;

  typedef typename Base::Has_left_category      Has_left_category;
    
  typedef typename Base::Point_2                Point_2;
  typedef typename Base::X_curve_2              X_curve_2;
  typedef X_curve_2                             Curve_2;

  // Obsolete, for backward compatibility
  typedef typename Base::Curve_point_status     Curve_point_status;
  typedef Point_2                               Point;
  typedef X_curve_2                             X_curve;
  typedef Curve_2                               Curve;

protected:
  typedef typename Kernel::Construct_vertex_2   Construct_vertex_2;
  typedef typename Kernel::Construct_segment_2  Construct_segment_2;
  typedef typename Kernel::Compare_xy_2         Compare_xy_2;

public:
  Arr_segment_traits_2() : Base() { }

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
    CGAL_precondition(curve_get_point_status(cv, split_pt) == ON_CURVE);
    CGAL_precondition_code(Compare_xy_2 compare_xy = compare_xy_2_object());
    CGAL_precondition(compare_xy(curve_source(cv), split_pt) != EQUAL);
    CGAL_precondition(compare_xy(curve_target(cv), split_pt) != EQUAL);
    
    Construct_vertex_2 construct_vertex = construct_vertex_2_object();
    Construct_segment_2 construct_segment = construct_segment_2_object();
    const Point_2 & source = construct_vertex(cv, 0);
    const Point_2 & target = construct_vertex(cv, 1);
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

  /*! nearest_intersection_to_left() finds the nearest intersection point of
   * two given curves to the left of a given point. Nearest is defined as the
   * lexicographically nearest not including the point itself with one
   * exception explained bellow..
   * If the intersection of the two curves is an X_curve_2, that is,
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
  bool nearest_intersection_to_left(const X_curve_2 & c1,
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
      if (compare_xy_2_object()(p1, pt) == SMALLER) {
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
   * \todo end point coincidence instead of intersection!
   */
  bool curves_overlap(const X_curve_2 & c1, const X_curve_2 & c2) const
  {
    Object res = intersect_2_object()(c1, c2);
    X_curve_2 seg;
    return (assign(seg, res) != 0);
  }

};

CGAL_END_NAMESPACE

#endif
