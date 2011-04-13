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
// release		 : 
// release_date  : 1999, October 13
//
// file 		 : include/CGAL/Td_traits.h
// package		 : Trapezoidal decomposition 2
// source		 : 
// revision 	 : 
// revision_date : 
// author(s)	 : Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// maintainer(s) : Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// coordinator	 : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter		 : 
// ======================================================================
#ifndef CGAL_TD_TRAITS_H
#define CGAL_TD_TRAITS_H

CGAL_BEGIN_NAMESPACE

template <class Pm_traits_,class X_curve_> class Td_traits : public Pm_traits_
{
public:
  typedef Pm_traits_ Traits_base;  
  typedef X_curve_ X_curve;
  typedef Td_traits<Traits_base,X_curve> Self;
  typedef typename Traits_base::Point Point;
  typedef X_curve* X_curve_ptr;
  typedef X_curve& X_curve_ref;
  typedef const X_curve& X_curve_const_ref;
  typedef Td_X_trapezoid<Self> X_trapezoid;
  typedef X_trapezoid* X_trapezoid_ptr;
  typedef X_trapezoid& X_trapezoid_ref;
  typedef const X_trapezoid& X_trapezoid_const_ref;
  
  Td_traits(const Traits_base& t) : Traits_base(t){}
  Td_traits() {}
  
protected:
  typedef X_trapezoid_const_ref const_ref;
  
public:
  /*
    note:
    The traits assume that the trapezoid is active,non empty,
    and planar, that is no two curves intersect in non degenerate curve.
  */

  inline bool curve_is_unbounded(const X_curve& cv) const {
/* compare curve with static unbounded curve */
    return cv.identical(CURVE_AT_INFINITY);
  }
inline bool trapezoid_bottom_curve_is_same(X_trapezoid_const_ref left,
					   X_trapezoid_const_ref right) const
/* returns true if bottom curves of input are the same */
{
  if (left.is_bottom_unbounded()) return right.is_bottom_unbounded();
  if (right.is_bottom_unbounded()) return false;
  return curve_is_same(left.bottom(),right.bottom());
}

inline bool trapezoid_top_curve_is_same(X_trapezoid_const_ref left,
					X_trapezoid_const_ref right) const
/* returns true if top curves of input are the same */
{
  if (left.is_top_unbounded()) return right.is_top_unbounded();
  if (right.is_top_unbounded()) return false;
  return curve_is_same(left.top(),right.top());
}

  bool is_degenerate(const_ref tr) const {
    return is_degenerate_point(tr) || 
      !tr.is_top_unbounded() && 
      !tr.is_bottom_unbounded() && 
      curve_is_same(tr.bottom(),tr.top());
  }		
  
  bool is_degenerate_point(const_ref tr) const
  {
    return !tr.is_left_unbounded() && !tr.is_right_unbounded() && 
      point_is_same(tr.left(),tr.right());
  }
  bool is_degenerate_curve(const_ref tr) const
  {
    return !tr.is_top_unbounded() && !tr.is_bottom_unbounded() && 
      curve_is_same(tr.bottom(),tr.top())&&!is_degenerate_point(tr);
  }
  bool is_vertical(const_ref tr) const
  {
    return !tr.is_left_unbounded() && 
      !tr.is_right_unbounded() && 
      point_is_same_x(tr.left(),tr.right());
  }
  
  /* Description:
     returns whether point is inside trapezoid using lexicographic order */
  
  bool is_inside(const_ref tr,const Point& p) const
  {
    return	
      (tr.is_left_unbounded()||
       point_is_left_low(tr.left(),p))&&
      (tr.is_right_unbounded()||
       point_is_right_top(tr.right(),p))&&
      (tr.is_bottom_unbounded()||
       curve_get_point_status(tr.bottom(),p)==Traits_base::ABOVE_CURVE)&&
      (tr.is_top_unbounded()||
       curve_get_point_status(tr.top(),p)==Traits_base::UNDER_CURVE);
  }
  bool is_in_closure(const_ref tr,const Point& p) const
  {
    // test left and right sides
    if ((tr.is_left_unbounded()||
         !point_is_left_low(p,tr.left()))&&
        (tr.is_right_unbounded()||
         !point_is_right_top(p,tr.right())))
      {
        // test bottom side
        if (!tr.is_bottom_unbounded()) 
          {
            typename Traits_base::Curve_point_status 
              s=curve_get_point_status(tr.bottom(),p);
            if (s!=Traits_base::ABOVE_CURVE&&s!=Traits_base::ON_CURVE)
              return false;
          }
        // test top side
        if (!tr.is_top_unbounded())
          {
            typename Traits_base::Curve_point_status 
              s=curve_get_point_status(tr.top(),p);
            if (s!=Traits_base::UNDER_CURVE&&s!=Traits_base::ON_CURVE)
              return false;
          }
        return true;
      }
    return false;
  }
  
  
public:
  static const Point& get_point_at_left_top_infinity();
  static const Point& get_point_at_right_bottom_infinity();
  static const X_curve& get_curve_at_infinity();
private:
  static Point POINT_AT_LEFT_TOP_INFINITY;
  static Point POINT_AT_RIGHT_BOTTOM_INFINITY;
  static X_curve CURVE_AT_INFINITY;
};

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Td_traits.C>
#endif

#endif //CGAL_TD_TRAITS_H








