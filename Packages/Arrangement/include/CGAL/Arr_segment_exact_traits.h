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
// file          : include/CGAL/Arr_segment_exact_traits.h
// package       : arr (1.03)
// source        :
// revision      :
// revision_date :
// author(s)     : Iddo Hanniel <hanniel@math.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin)
// chapter       : Arrangement_2
//
// ======================================================================
#ifndef CGAL_ARR_SEGMENT_EXACT_TRAITS_H
#define CGAL_ARR_SEGMENT_EXACT_TRAITS_H

#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>

#include <list>

CGAL_BEGIN_NAMESPACE

template <class Kernel_>
class Arr_segment_exact_traits : public Pm_segment_traits_2<Kernel_>
{
public:
  typedef Kernel_                       Kernel;

  typedef int                           Info_face;
  typedef int                           Info_edge;
  typedef int                           Info_vertex;
  
  typedef Pm_segment_traits_2<Kernel>   Base;
  
  typedef typename Base::Point_2        Point_2;
  typedef typename Base::X_curve_2      X_curve_2;
  typedef X_curve_2                     Curve_2;

  // Obsolete, for backward compatibility
  typedef typename Base::Curve_point_status     Curve_point_status;
  typedef Point_2                       Point;
  typedef X_curve_2                     X_curve;
  typedef Curve_2                       Curve;

public:
  Arr_segment_exact_traits() : Base() { }

  bool is_x_monotone(const Curve_2 &) {return true;}
  
  // segments are x_monotone:
  void make_x_monotone(const Curve_2 & /*cv*/, std::list<Curve_2>& /*l*/) {} 

  /*!
   */
  X_curve_2 curve_flip(const X_curve_2 & cv) const {
    Construct_vertex_2 construct_vertex = construct_vertex_2_object();
    return X_curve_2(construct_vertex(cv, 1), construct_vertex(cv, 0));
  }
 
  /*!
   */
  void curve_split(const X_curve_2 & cv, X_curve_2 & c1, X_curve_2 & c2, 
                   const Point_2 & split_pt)
  {
    Construct_vertex_2 construct_vertex = construct_vertex_2_object();

    //split curve at split point (x coordinate) into c1 and c2
    CGAL_precondition(curve_get_point_status(cv, split_pt) == ON_CURVE);
    CGAL_precondition(compare_lexicographically_xy(curve_source(cv), split_pt)
		      != EQUAL);
    CGAL_precondition(compare_lexicographically_xy(curve_target(cv), split_pt)
		      != EQUAL);
    
    c1 = X_curve_2(curve_source(cv), split_pt);
    c2 = X_curve_2(split_pt, curve_target(cv));
  }

  /*!
   */
  //returns true iff the intersection is strictly right of pt
  bool do_intersect_to_right(const X_curve_2 & c1, const X_curve_2 & c2,
                             const Point_2 & pt) const 
  {
    Point_2 ip;
    X_curve_2 seg;
    Object res = intersection(c1, c2);

    if (assign(ip, res)) {
      return (compare_lexicographically_xy(ip, pt) == LARGER);
    }
    if (assign(seg, res)) {
      return ((compare_lexicographically_xy(seg.source(), pt) == LARGER) ||
              (compare_lexicographically_xy(seg.target(), pt) == LARGER) );
    }
     
    return false; //don't intersect at all
  }

  /*!
   */
  bool nearest_intersection_to_right(const X_curve_2 & c1,
                                     const X_curve_2 & c2,
                                     const Point_2 & pt,
                                     Point_2 & p1, Point_2 & p2) const
  {
    if (!do_intersect_to_right(c1, c2, pt)) return false;
    X_curve_2 seg;
    Object res = intersection(c1, c2);
    if (assign(seg, res)) {
      //p1, p2 will always be ordered left,right (make seg left to right)
      if (compare_lexicographically_xy(curve_source(seg), curve_target(seg))
	  == LARGER)
          seg = curve_flip(seg);

      if (compare_lexicographically_xy(curve_target(seg), pt) == LARGER) {
        p2 = curve_target(seg);
        if (compare_lexicographically_xy(curve_source(seg), pt) == LARGER)
          p1 = curve_source(seg);
        else
          p1 = pt;
        return true;
      }
      else {
        return false;
      }
    }
    
    if (assign(p1,res)) {
      if (compare_lexicographically_xy(p1, pt) == LARGER) {
        p2 = p1;
        return true;
      }

      return false;
    }

    return false;
  }

  /*!
   */
  X_curve_2 curve_reflect_in_x_and_y(const X_curve_2 & cv) const
  {
    // use hx(), hy(), hw() in order to support both Homogeneous and Cartesian
    X_curve_2 reflected_cv(point_reflect_in_x_and_y ( cv.source()),
                           point_reflect_in_x_and_y ( cv.target()));
    return reflected_cv;
  }


  /*!
   */
  Point_2 point_reflect_in_x_and_y(const Point_2 & pt) const
  {
    // use hx(), hy(), hw() in order to support both Homogeneous and Cartesian
    Point_2 reflected_pt(-pt.hx(), -pt.hy(), pt.hw());
    return reflected_pt;
  }

  /*!
   */
  //the following function is needed to deal with overlaps
  bool curves_overlap(const X_curve_2 & cv1, const X_curve_2 & cv2) const 
  {
    X_curve_2 seg;
    Object res = intersection(cv1, cv2);
    return (assign(seg, res) != 0);
  }

};

CGAL_END_NAMESPACE

#endif
