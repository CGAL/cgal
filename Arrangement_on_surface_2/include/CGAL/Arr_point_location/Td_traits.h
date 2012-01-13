// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)	 : Oren Nechushtan <theoren@math.tau.ac.il>
#ifndef CGAL_TD_TRAITS_H
#define CGAL_TD_TRAITS_H

#include <CGAL/Arr_point_location/Td_X_trapezoid.h>

namespace CGAL {

template <class Pm_traits_,class X_curve_> class Td_traits : public Pm_traits_
{
public:
  typedef Pm_traits_                      Traits_base;  
  typedef X_curve_                        X_curve;
  typedef Td_traits<Traits_base,X_curve>  Self;
  typedef typename Traits_base::Point_2   Point;
  typedef X_curve*                        X_curve_ptr;
  typedef X_curve&                        X_curve_ref;
  typedef const X_curve&                  X_curve_const_ref;
  typedef Td_X_trapezoid<Self>            X_trapezoid;
  typedef X_trapezoid*                    X_trapezoid_ptr;
  typedef X_trapezoid&                    X_trapezoid_ref;
  typedef const X_trapezoid&              X_trapezoid_const_ref;
  
  Td_traits(const Traits_base& t) 
    : Traits_base(t)
  {
  }

  Td_traits() 
  {
  }

  ~Td_traits(void)
  {
    if (POINT_AT_LEFT_TOP_INFINITY) {
      delete POINT_AT_LEFT_TOP_INFINITY;
      POINT_AT_LEFT_TOP_INFINITY = 0;
    }
    if (POINT_AT_RIGHT_BOTTOM_INFINITY) {
      delete POINT_AT_RIGHT_BOTTOM_INFINITY;
      POINT_AT_RIGHT_BOTTOM_INFINITY = 0;
    }
    if (CURVE_AT_INFINITY) {
      delete CURVE_AT_INFINITY;
      CURVE_AT_INFINITY = 0;
    }
  }

protected:
  typedef X_trapezoid_const_ref           const_ref;
  
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

  inline bool trapezoid_bottom_curve_equal(X_trapezoid_const_ref left,
					   X_trapezoid_const_ref right) const
  /* returns true if bottom curves of input are the same */
  {
    if (left.is_bottom_unbounded())
      return (right.is_bottom_unbounded());
  
    if (right.is_bottom_unbounded()) 
      return (false);
    
    return (this->equal_2_object()(left.bottom(),right.bottom()));
  }

  inline bool trapezoid_top_curve_equal(X_trapezoid_const_ref left,
					X_trapezoid_const_ref right) const
  /* returns true if top curves of input are the same */
  {
    if (left.is_top_unbounded()) 
      return (right.is_top_unbounded());
    
    if (right.is_top_unbounded()) 
      return (false);
    
    return (this->equal_2_object()(left.top(),right.top()));
  }

  //returns true if the trapezoid is a point or a curve
  bool is_degenerate(const_ref tr) const
  {
    return (is_degenerate_point(tr) || 
	    (!tr.is_top_unbounded() && 
	     !tr.is_bottom_unbounded() && 
	     this->equal_2_object()(tr.bottom(),tr.top())));
  }		
  
  //returns true if the trapezoid is a point 
  bool is_degenerate_point(const_ref tr) const
  {
    return (!tr.is_left_unbounded() && 
	    !tr.is_right_unbounded() && 
	    this->equal_2_object()(tr.left(),tr.right()));
  }

  //returns true if the trapezoid is a curve 
  bool is_degenerate_curve(const_ref tr) const
  {
    return (!tr.is_top_unbounded() && 
	    !tr.is_bottom_unbounded() && 
	    this->equal_2_object()(tr.bottom(), tr.top()) && 
	    !is_degenerate_point(tr));
  }

  //returns true if the trapezoid is vertical 
  bool is_vertical(const_ref tr) const
  {
    return (!tr.is_left_unbounded() && 
	    !tr.is_right_unbounded() && 
	    (this->compare_x_2_object()(tr.left(),tr.right())== EQUAL));
  }
  
  /* Description:
     returns whether point is inside trapezoid using lexicographic order */
    bool is_inside(const_ref tr,const Point& p) const
  {
    return	
      (tr.is_left_unbounded() ||
       (this->compare_xy_2_object()(tr.left(),p)==SMALLER)) &&
      (tr.is_right_unbounded() ||
       (this->compare_xy_2_object()(tr.right(),p)==LARGER)) &&
      (tr.is_bottom_unbounded() ||
       this->compare_y_at_x_2_object()(p, tr.bottom()) == LARGER) &&
      (tr.is_top_unbounded()||
       this->compare_y_at_x_2_object()(p, tr.top()) == SMALLER);
  }

  // returns true if the point is inside the closure of the trapezoid 
  // (inlcude all boundaries)
  bool is_in_closure(const_ref tr,const Point& p) const
  {
    // test left and right sides
    if ((tr.is_left_unbounded()||
         !(this->compare_xy_2_object()(p,tr.left())==SMALLER)) &&
        (tr.is_right_unbounded()||
         !(this->compare_xy_2_object()(p,tr.right())==LARGER)))
      {
        // test bottom side
        if (!tr.is_bottom_unbounded()) 
	{
	  if (this->compare_y_at_x_2_object()(p, tr.bottom()) == SMALLER)
	    return false;
	}
        // test top side
        if (!tr.is_top_unbounded())
	{
	  if (this->compare_y_at_x_2_object()(p, tr.top()) == LARGER)
	    return false;
	}
        return true;
      }
    return false;
  }
  
  
public:
  static const Point& point_at_left_top_infinity();
  static const Point& point_at_right_bottom_infinity();
  static const X_curve& curve_at_infinity();
private:
  static Point * POINT_AT_LEFT_TOP_INFINITY;
  static Point * POINT_AT_RIGHT_BOTTOM_INFINITY;
  static X_curve * CURVE_AT_INFINITY;
};

} //namespace CGAL

#include <CGAL/Arr_point_location/Td_traits_functions.h>

#endif
