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
#include <CGAL/Pm_segment_exact_traits.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>

#include <list>

CGAL_BEGIN_NAMESPACE

template <class R>
class Arr_segment_exact_traits : public Pm_segment_exact_traits<R>
{
public:
  typedef int			Info_face;
  typedef int	      		Info_edge;
  typedef int	       	       	Info_vertex;
  
  typedef Pm_segment_exact_traits<R> Base;
  
  typedef typename Base::Curve_status Curve_status;
  typedef typename Base::Curve_point_status Curve_point_status;
  
  typedef typename Base::Point Point;
  typedef typename Base::X_curve X_curve;

typedef Segment_2<R>      	Curve;

  
public:
  Arr_segment_exact_traits() : Base() 
  { }

  bool is_x_monotone(const Curve&) {return true;}
  
  //segments are x_monotone :
  void make_x_monotone(const Curve& /*cv*/, std::list<Curve>& /*l*/) {} 

  X_curve curve_flip(const X_curve& cv) const {
    return X_curve(cv.target(), cv.source());
  }
 
  void curve_split(const X_curve& cv, X_curve& c1, X_curve& c2, 
                   const Point& split_pt)
  {
    //split curve at split point (x coordinate) into c1 and c2
    CGAL_precondition(curve_get_point_status(cv,split_pt)==ON_CURVE);
    CGAL_precondition(compare_lexicographically_xy(curve_source(cv),split_pt)
		      != EQUAL);
    CGAL_precondition(compare_lexicographically_xy(curve_target(cv),split_pt)
		      != EQUAL);
    
    c1=X_curve(curve_source(cv),split_pt);
    c2=X_curve(split_pt,curve_target(cv));
  }

  //returns true iff the intersection is strictly right of pt
  bool do_intersect_to_right(const X_curve& c1, const X_curve& c2,
                             const Point& pt) const 
  {
    Object res;
    Point ip;
    X_curve seg;
    res=intersection(c1,c2);

    if (assign(ip,res)) {
      return (compare_lexicographically_xy(ip,pt)==LARGER);
    }
    if (assign(seg,res)) {
      return ( (compare_lexicographically_xy(seg.source(),pt)==LARGER) ||
               (compare_lexicographically_xy(seg.target(),pt)==LARGER) );
    }
     
    return false; //don't intersect at all
  }


  bool nearest_intersection_to_right(const X_curve& c1,
                                      const X_curve& c2,
                                      const Point& pt,
                                     Point& p1,
                                     Point& p2) const
  {
    if (!do_intersect_to_right(c1,c2,pt)) return false;
    Object res;
    X_curve seg;
    res=intersection(c1,c2);
    if (assign(seg,res)) {
      //p1, p2 will always be ordered left,right (make seg left to right)
      if (compare_lexicographically_xy(curve_source(seg), curve_target(seg))
	  == LARGER)
          seg=curve_flip(seg);

      if (compare_lexicographically_xy(curve_target(seg),pt)==LARGER) {
        p2=curve_target(seg);
        if (compare_lexicographically_xy(curve_source(seg),pt)==LARGER)
          p1=curve_source(seg);
        else
          p1=pt;
        return true;
      }
      else {
        return false;
      }
    }
    
    if (assign(p1,res)) {
      if (compare_lexicographically_xy(p1,pt)==LARGER) {
        p2=p1;
        return true;
      }

      return false;
    }

    return false;
  }

  X_curve curve_reflect_in_x_and_y (const X_curve& cv) const
  {
    // use hx(), hy(), hw() in order to support both Homogeneous and Cartesian
    X_curve reflected_cv( point_reflect_in_x_and_y ( cv.source()),
			  point_reflect_in_x_and_y ( cv.target()));
    return reflected_cv;
  }


  Point point_reflect_in_x_and_y (const Point& pt) const
  {
    // use hx(), hy(), hw() in order to support both Homogeneous and Cartesian
    Point reflected_pt( -pt.hx(), -pt.hy(), pt.hw());
    return reflected_pt;
  }

  //the following function is needed to deal with overlaps
  bool curves_overlap(const X_curve& cv1, const X_curve& cv2) const 
  {
    Object res;
    X_curve seg;
    res=intersection(cv1,cv2);
    return (assign(seg,res)!=0);
  }

};

CGAL_END_NAMESPACE

#endif // CGAL_PM_SEGMENT_EPSILON_TRAITS_H



















