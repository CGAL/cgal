// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-26 $
// release_date  : $CGAL_Date: 2001/01/05 $
//
// file          : include/CGAL/Pm_leda_segment_traits_2.h
// package       : pm (5.43)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Eyal flato  	  <flato@math.tau.ac.il>
//		   Iddo hanniel	<hanniel@math.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_PM_LEDA_SEGMENT_TRAITS_2_H
#define CGAL_PM_LEDA_SEGMENT_TRAITS_2_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#include <CGAL/rat_leda_in_CGAL_2.h>

#include <cmath>

#ifdef CGAL_PM_LEDA_SEGMENT_EXACT_TRAITS_DEBUG
#include <CGAL/Pm_traits_operations_count.h>
#define OP_CONSUME(x, y) OP.Consume(x, y)
#else
#define OP_CONSUME(x, y) 
#endif

CGAL_BEGIN_NAMESPACE

class Pm_leda_segment_traits_2
{
public:
	
  typedef leda_rat_segment X_curve_2;
  typedef leda_rat_point   Point_2;
  typedef leda_rat_vector  Vector_2;

  // Obsolete, for backward compatibility
  typedef leda_rational	   NumberType;
  typedef X_curve_2        X_curve;
  typedef Point_2          Point;
  typedef Vector_2         Vector;
	
protected:
  typedef enum 
    {
      CURVE_VERTICAL_UP = 0,
      CURVE_VERTICAL_DOWN = 2,
      CURVE_LEFT = 3,
      CURVE_RIGHT = 1
    } Curve_status;
public:	

  typedef enum
    {
      UNDER_CURVE = -1,
      ABOVE_CURVE = 1,
      ON_CURVE = 2,
      CURVE_NOT_IN_RANGE = 0
      //CURVE_VERTICAL = 3
    } Curve_point_status;	
	
#ifdef CGAL_PM_LEDA_SEGMENT_EXACT_TRAITS_DEBUG
  COperationsConsumptionCount OP;
#endif
	
public:
	
  Pm_leda_segment_traits_2()
#ifdef CGAL_PM_LEDA_SEGMENT_EXACT_TRAITS_DEBUG
    : OP(opNumOperations)
#endif
  {
  }
	
  ~Pm_leda_segment_traits_2() 
  {
#ifdef CGAL_PM_LEDA_SEGMENT_EXACT_TRAITS_DEBUG
    std::ofstream os("OpCount.txt");
    OP.Write(os);
#endif
  }
	
  Point_2 curve_source(const X_curve_2 & cv) const 
  { 
    return cv.source(); 
  }
	
  Point_2 curve_target(const X_curve_2 & cv) const 
  {
    return cv.target(); 
  }
	
  bool curve_is_vertical(const X_curve_2 & cv) const 
  {
    OP_CONSUME(opCompare, 1);
    if (cv.is_trivial())
      return true;
    return cv.is_vertical();
  }
	
  bool curve_is_in_x_range(const X_curve_2 & cv, const Point_2 & q) const
  { 
    OP_CONSUME(opCompare, 3);
    return !( is_right(q, rightmost(cv.source(), cv.target())) ||
              is_left(q, leftmost(cv.source(), cv.target()))	 );
  }
	
  bool curve_is_in_y_range(const X_curve_2 &cv, const Point_2 & q) const
  { 
    OP_CONSUME(opCompare, 3);
    bool r = !( is_lower(q, lowest(cv.source(), cv.target())) ||
                is_higher(q, highest(cv.source(), cv.target())) );
    return r;
  }
	
  Orientation orientation(const Point_2 &p, const Point_2 &q, const Point_2 &r)
    const
  {
    OP_CONSUME(opOrientation, 1);
    return CGAL::orientation(p, q, r);
  }
	
  Curve_point_status 
  curve_get_point_status(const X_curve_2 &cv, const Point_2 & p) 
    const
  {
    if (!curve_is_in_x_range(cv, p))
      return CURVE_NOT_IN_RANGE;
    if (!curve_is_vertical(cv))
    {
      Orientation o;
      if (compare_x(cv.source(), cv.target()) < 0)
        o = orientation(cv.source(), cv.target(), p);
      else
        o = orientation(cv.target(), cv.source(), p);
			
      if (o < 0)
        return UNDER_CURVE;
      if (o > 0)
        return ABOVE_CURVE;
      if (o == 0)
        return ON_CURVE;
      return ON_CURVE;
    }
    else
    {
      //return CURVE_VERTICAL; 
      //&&& - copied from exact-traits ask Oren & Iddo 1/4/99
			
      if (is_lower(p,lowest(cv.source(),cv.target())))
        return UNDER_CURVE;
      if (is_higher(p,highest(cv.source(),cv.target())))
        return ABOVE_CURVE;
      // if (curve_is_in_y_range(cv,p))
      return ON_CURVE;
    }                    
  }
	
  Comparison_result 
  curve_compare_at_x(const X_curve_2 &cv1, const X_curve_2 &cv2, 
                     const Point_2 &q) 
    const 
  {
    if ((!curve_is_in_x_range(cv1, q)) || 
        (!curve_is_in_x_range(cv2, q)))
      return EQUAL;
		
    int res;
    // bug ??? in LEDA - 
    // cmp_segments_at_xcoord returns wrong answer if
    // cv1 (or cv2) are from right to left
    // cv1_ and cv2_ are the same as cv1 and cv2 - 
    //   oriented from left to right
    X_curve_2 cv1_ = cv1;
    X_curve_2 cv2_ = cv2;
    if (lexicographically_xy_larger(cv1.source(), cv1.target()))
      cv1_ = cv1.reversal();
    if (lexicographically_xy_larger(cv2.source(), cv2.target()))
      cv2_ = cv2.reversal();
		
                
    // checking verical curves.
    if (curve_is_vertical(cv1_)) {
                  
      if (curve_is_vertical(cv2_)) {
        // both cv1 and cv2 are vertical
        if ( is_lower(cv1_.target(), cv2_.source()) )
          return SMALLER;
        if ( is_higher(cv1_.source(), cv2_.target()) )
          return LARGER;
        return EQUAL; // overlapping. 
      } // end  both cv1 and cv2 are vertical.
      else { // only cv1 is vertical.
        if (orientation(cv2_.source(), 
                        cv2_.target(), 
                        cv1_.source()) > 0 )
          return LARGER;
                    
        if (orientation(cv2_.source(), 
                        cv2_.target(), 
                        cv1_.target()) < 0)
          return SMALLER;

        return EQUAL;
      }
    }
                
    if (curve_is_vertical(cv2_)) { // only cv2 is vertical.
      if (orientation(cv1_.source(), 
                      cv1_.target(), 
                      cv2_.source()) > 0 )
        return SMALLER;
                    
      if (orientation(cv1_.source(), 
                      cv1_.target(), 
                      cv2_.target()) < 0)
        return LARGER;

      return EQUAL;  
    }
                  
    // end checking verical curves.

    OP_CONSUME(opSegmentX, 1);
    res = cmp_segments_at_xcoord(cv1_, cv2_, q);
		
    if (res < 0) 
      return SMALLER;
    if (res > 0) 
      return LARGER;
    return EQUAL;
  }
	
	
  Comparison_result 
  curve_compare_at_x_left(const X_curve_2 &cv1, 
                          const X_curve_2 &cv2, 
                          const Point_2 &q) const 
  {
    if (curve_is_vertical(cv1) || (curve_is_vertical(cv2))) 
      return EQUAL;
    if (!is_left(leftmost(cv1.source(), cv1.target()), q)) 
      return EQUAL;
    if (!is_left(leftmost(cv2.source(), cv2.target()), q)) 
      return EQUAL;
		
    Comparison_result r = curve_compare_at_x(cv1, cv2, q);
		
    if ( r != EQUAL)
      return r;     // since the curve is continous 
		
    // <cv2> and <cv1> meet at a point with the same x-coordinate 
    // as q
    return (Comparison_result)
      orientation(leftmost(cv2.source(), cv2.target()), q,
                  leftmost(cv1.source(), cv1.target()) );
  }
	
  Comparison_result 
  curve_compare_at_x_right(const X_curve_2 &cv1, 
                           const X_curve_2 &cv2, 
                           const Point_2 & q) const 
  {
    if (curve_is_vertical(cv1) || (curve_is_vertical(cv2))) 
      return EQUAL;
    if (!is_right(rightmost(cv1.source(), cv1.target()), q)) 
      return EQUAL;
    if (!is_right(rightmost(cv2.source(), cv2.target()), q)) 
      return EQUAL;
		
    Comparison_result r = curve_compare_at_x(cv1, cv2, q);
		
    if ( r != EQUAL)
      return r;     // since the curve is continous (?)
		
    // <cv1> and <cv2> meet at a point with the same x-coordinate 
    // as q
    return (Comparison_result)
      orientation(q, rightmost(cv2.source(), cv2.target()),  
                  rightmost(cv1.source(), cv1.target()) );
  }
	
	
  X_curve_2 curve_flip(const X_curve_2 &cv) const
  {
    return cv.reversal();
  }
	

  bool curve_is_between_cw(const X_curve_2 &cv, 
                           const X_curve_2 &first, 
                           const X_curve_2 &second, 
                           const Point_2 &point)	const
    // TRUE if cv is between first and second in cw direction
    // precondition: this, first and second have a common endpoint
    // precondition: first, second, this are pairwise interior 
    //disjoint
  {
    Point_2 p = point;
		
    X_curve_2 cv0 = first;
    X_curve_2 cv1 = second;
    X_curve_2 cvx = cv;
		
    if ( !is_same(cv0.source(),p) ) cv0 = curve_flip(cv0);
    if ( !is_same(cv1.source(),p) ) cv1 = curve_flip(cv1);
    if ( !is_same(cvx.source(),p) ) cvx = curve_flip(cvx);
		
    Point_2 p0(cv0.target()), p1(cv1.target()), px(cvx.target());
		
    int or0=orientation(p0,p,px);
    int or1=orientation(p1,p,px);
    // Bug Fix, Shai, Jan, 8, 2001
    // 'or' is a keyword in C++, changed to 'orient'
    int orient=or0*or1;                        
		
    if (orient < 0) 
    { //one is a left_turn the other right_turn
      return (or0 == LEFT_TURN); //left_turn
    }
    else 
    { //both are either left or right turns (or one is colinear)
      if (orientation(p0,p,p1)==RIGHT_TURN)
      {
        if ((or1 == 0) && (or0 == RIGHT_TURN)) 
          // the case where cvx and cv1 are colinear
          return false; 
        else 
          return true;
      }
      else
      {
        return false;
      }
    }
  }

  Comparison_result Comparison_result_from_int(int res) const 
  {
    if (res < 0) return SMALLER;
    if (res > 0) return LARGER;
    return EQUAL;
  }
  
  Comparison_result compare_x(const Point_2 &p1, const Point_2 &p2) const
  {
    OP_CONSUME(opCompare, 1);
    return Comparison_result_from_int(

#if (__LEDA__ >= 380)
      Pm_leda_segment_traits_2::Point_2::cmp_x(p1,p2)
#else // backward compatability to LEDA
      compare(p1.xcoord(),p2.xcoord())
#endif

                                      );
    /*static double d;
      d = p1.xcoordD() - p2.xcoordD();
      if (abs(d) >= 0.000001)
      if (d > 0)
      return LARGER;
      else
      return SMALLER;
      return Comparison_result_from_int(compare(p1.xcoord(), p2.xcoord()));
    */
  }
  
  Comparison_result compare_y(const Point_2 &p1, const Point_2 &p2) const
  {
    OP_CONSUME(opCompare, 1);
    return Comparison_result_from_int(

#if (__LEDA__ >= 380)
      Pm_leda_segment_traits_2::Point_2::cmp_y(p1,p2)
#else // backward compatability to LEDA   
      compare(p1.ycoord(),p2.ycoord()) 
#endif

                                      );

    /*static double d;
      d = p1.ycoordD() - p2.ycoordD();
      if (abs(d) >= 0.000001)
      if (d > 0)
      return LARGER;
      else
      return SMALLER;
      return Comparison_result_from_int(compare(p1.ycoord(), p2.ycoord()));
    */
  }
  
public:
  Point_2 point_to_left(const Point_2& p) const  
  {return Point_2(p.xcoord()-1, p.ycoord());}
	
  Point_2 point_to_right(const Point_2& p) const 
  {return Point_2(p.xcoord()+1, p.ycoord());}
	
  bool curve_is_same(const X_curve_2 &cv1, const X_curve_2 &cv2) const
  {
    OP_CONSUME(opCompare, 4);
    return (cv1 == cv2) != 0;
  }
	
  bool is_point_on_curve(const X_curve_2 &cv, const Point_2& p) const //check
  {
    //return (cv.sqr_dist(p) <= 0);
    OP_CONSUME(opSegmentX, 1);
    return cv.contains(p);
  }
	
public:
  bool is_left(const Point_2 &p1, const Point_2 &p2) const 
  { 
    OP_CONSUME(opCompare, 1); 
    return (compare_x(p1, p2) == SMALLER); 
  }
	
  bool is_right(const Point_2 &p1, const Point_2 &p2) const 
  { 
    OP_CONSUME(opCompare, 1);  
    return (compare_x(p1, p2) == LARGER);
  }
	
  bool is_same_x(const Point_2 &p1, const Point_2 &p2) const 
  { 
    OP_CONSUME(opCompare, 1);  
    return (compare_x(p1, p2) == EQUAL);
  }
	
  bool is_lower(const Point_2 &p1, const Point_2 &p2) const 
  { 
    OP_CONSUME(opCompare, 1);  
    return (compare_y(p1, p2) == SMALLER);
  }
	
  bool is_higher(const Point_2 &p1, const Point_2 &p2) const 
  { 
    OP_CONSUME(opCompare, 1); 
    return (compare_y(p1, p2) == LARGER); 
  }
	
  bool is_same_y(const Point_2 &p1, const Point_2 &p2) const 
  { 
    OP_CONSUME(opCompare, 1); 
    return (compare_y(p1, p2) == EQUAL); 
  }
	
  bool is_same(const Point_2 &p1, const Point_2 &p2) const
  { 
    OP_CONSUME(opCompare, 2); 
    return (p1 == p2); 
  }
	
  const Point_2& leftmost(const Point_2 &p1, const Point_2 &p2) const
  { return (is_left(p1, p2) ? p1 : p2); }
	
  const Point_2& rightmost(const Point_2 &p1, const Point_2 &p2) const
  { return (is_right(p1, p2) ? p1 : p2); }
	
  const Point_2& lowest(const Point_2 &p1, const Point_2 &p2) const
  { return (is_lower(p1, p2) ? p1 : p2); }
	
  const Point_2& highest(const Point_2 &p1, const Point_2 &p2) const
  { return (is_higher(p1, p2) ? p1 : p2); }
	
  Curve_status curve_get_status(const X_curve_2 &cv) const
  {
    if (curve_is_vertical(cv)) 
    {
      if ( is_higher(cv.target(), cv.source()) )
        return CURVE_VERTICAL_UP;
      else
        return CURVE_VERTICAL_DOWN;
    }
    else
    {
      if ( is_right(cv.target(), cv.source()) )
        return CURVE_RIGHT;
      else
        return CURVE_LEFT;
    }
  }
	
};

CGAL_END_NAMESPACE

#endif // CGAL_PM_LEDA_SEGMENT_TRAITS_2_H
// EOF
