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
// release       : $$
// release_date  : $$
//
// file          : include/CGAL/Arr_conic_traits_2.h
// package       : Arrangement (2.62)
// author(s)     : Ron Wein <wein@post.tau.ac.il>
// 
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_ARR_CONIC_TRAITS_2_H
#define CGAL_ARR_CONIC_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include "CGAL/Arrangement_2/Conic_arc_2.h"
#include <list>

CGAL_BEGIN_NAMESPACE

// ----------------------------------------------------------------------------
// Arrangement traits for conic arcs.
//

template <class Kernel_>
class Arr_conic_traits_2 
{
 public:
  typedef Kernel_                       Kernel;
  typedef typename Kernel::FT           NT;

  // Categories:
//#define HAS_LEFT_NOT
#if !defined(HAS_LEFT_NOT)
  typedef Tag_true                        Has_left_category;
#else
  typedef Tag_false                       Has_left_category;
#endif
    
  // The difference between Curve_2 and X_curve_2 is semantical only,
  // NOT syntactical.
  typedef Conic_arc_2<Kernel>        Curve_2;
  typedef Curve_2                    X_curve_2;
 
  // Using typename to please compiler (e.g., CC with IRIX64 on mips)
  typedef typename Curve_2::Kernel    R;
  typedef typename Curve_2::Point_2   Point_2;
  typedef typename Curve_2::Circle_2  Circle_2;
  typedef typename Curve_2::Segment_2 Segment_2;

  // For backward compatibility:
  typedef Curve_2                    Curve;
  typedef X_curve_2                  X_curve;
  typedef Point_2                    Point;
  typedef Circle_2                   Circle;
  typedef Segment_2                  Segment;

  enum Curve_point_status
  {
    UNDER_CURVE = -1,
    ABOVE_CURVE = 1,
    ON_CURVE = 2,
    CURVE_NOT_IN_RANGE = 0,
  };

#ifdef CGAL_CONIC_ARC_USE_CACHING
 private:

  typedef typename Curve_2::Intersections Intersections;
  mutable std::list<Intersections> inter_list; // For caching intersections.
#endif

 public:

  // Constructor.
  Arr_conic_traits_2()
  {}

  ////////// Planar Map methods: //////////

  // Compare the co-ordinates of two given points.
  Comparison_result compare_x(const Point_2& p0, const Point_2& p1) const
  {
    return (p0.compare_x(p1));
  }

  Comparison_result compare_xy(const Point_2& p0, const Point_2& p1) const
  {
    Comparison_result  x_res = p0.compare_x(p1);
    
    if (x_res != EQUAL)
      return (x_res);

    return (p0.compare_y(p1));
  }

  // Check whether the given curve is a vertical segment.
  bool curve_is_vertical(const X_curve_2& curve) const
  {
    return (curve.is_vertical_segment());
  }

  // Check whether the x co-ordinate of the given point is contained in the
  // x-range of the given x-monotone curve.
  bool curve_is_in_x_range (const X_curve_2& curve, const Point_2& p) const
  {
    CGAL_precondition(is_x_monotone(curve));

    if (curve.is_vertical_segment())
    {
      // Check if the vertical segment's x co-ordinate is the same as p's.
      return (compare_x (curve.source(), p) == EQUAL);
    }
    else
    {
      // Since the curve is x-monotone, if the point x co-ordinate is to the
      // left (or to the right of both curve's source and target points), then
      // the point is obviously not in the curve's x range.
      return (compare_x(p, curve.source()) != compare_x(p, curve.target()));
    }
  }

  // Decide wether curve1 is above, below or equal to curve2 at the
  // x co-ordinate of the given point.
  Comparison_result curve_compare_at_x (const X_curve_2& curve1, 
				        const X_curve_2& curve2, 
				        const Point_2& p) const
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    // Get the points on curve1 with the same x co-ordinate as p.
    int    n1;
    Point_2  ps1[2];

    if (curve1.is_vertical_segment())
    {
      // curve1 is a vertical segment.
      if (compare_x (curve1.source(), p) != EQUAL)
      {
	// The curve's x co-ordinate is different than p's.
	return (EQUAL);
      }
      else
      {
	// Both the source and target of the curve has the same x co-ordinate
	// as p. Make sure that ps1[0] has a smaller y value than ps1[1].
	n1 = 2;
	if (_compare_y (curve1.source(), curve1.target()) == SMALLER)
	{
	  ps1[0] = curve1.source();
	  ps1[1] = curve1.target();
	}
	else
	{
	  ps1[0] = curve1.target();
	  ps1[1] = curve1.source();
	}
      }
    }
    else
    {
      // Find all points on curve1 with the same x co-ordinate as p.
      n1 = curve1.get_points_at_x (p, ps1);

      if (n1 == 0)    // p is not in the x-range of curve1.
	return (EQUAL);

      CGAL_assertion(n1 == 1);
    }

    // Get the points on curve2 with the same x co-ordinate as p.
    int    n2;
    Point_2  ps2[2];

    if (curve2.is_vertical_segment())
    {
      // curve2 is a vertical segment.
      if (compare_x (curve2.source(), p) != EQUAL)
      {
	// The curve's x co-ordinate is different than p's.
	return (EQUAL);
      }
      else
      {
	// Both the source and target of the curve has the same x co-ordinate
	// as p. Make sure that ps2[0] has a smaller y value than ps2[1].
	n2 = 2;
	if (_compare_y (curve2.source(), curve2.target()) == SMALLER)
	{
	  ps2[0] = curve2.source();
	  ps2[1] = curve2.target();
	}
	else
	{
	  ps2[0] = curve2.target();
	  ps2[1] = curve2.source();
	}
      }
    }
    else
    {
      // Find all points on curve2 with the same x co-ordinate as p.
      n2 = curve2.get_points_at_x (p, ps2);

      if (n2 == 0)    // p is not in the x-range of curve2.
	return (EQUAL);

      CGAL_assertion(n2 == 1);
    }

    // Deal with vertical segments:
    if (n1 == 2)
    {
      // Check if the vertical segment curve1 contains ps2[0] or ps2[1].
      if (_compare_y (ps1[0], ps2[0]) != LARGER && 
	  _compare_y (ps1[1], ps2[0]) != SMALLER)
      {
	return (EQUAL);
      }

      if (n2 == 2)
      {
	if (_compare_y (ps1[0], ps2[1]) != LARGER && 
	    _compare_y (ps1[1], ps2[1]) != SMALLER)
	{
	  return (EQUAL);
	}
      }
    }
    else if (n2 == 2)
    {
      // Check if the vertical segment curve2 contains ps1[0].
      if (_compare_y (ps2[0], ps1[0]) != LARGER && 
	  _compare_y (ps2[1], ps1[0]) != SMALLER)
      {
	return (EQUAL);
      }
    }
    
    // None of the curves is a vertical segments and both have exactly
    // one point with the given x co-ordinate:
    // Compare the y co-ordinates of these two points.
    return (_compare_y (ps1[0], ps2[0]));
  }

  // Decide wether curve1 is above, below or equal to curve2 immediately to
  // the left of the x co-ordinate of the given point.
  Comparison_result curve_compare_at_x_left (const X_curve_2& curve1, 
					     const X_curve_2& curve2,
					     const Point_2& p) const
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    // If the curves are the same, they are equal to the left of p:
    if (curve_is_same(curve1,curve2))
      return (EQUAL);

    // Check the case of vertical segments: we assume that a vertical segment
    // does not contain a point immediately to the left of p.
    if (curve1.is_vertical_segment() || curve2.is_vertical_segment())
      return (EQUAL);

    // Get the points on curve1 with the same x co-ordinate as p.
    int    n1;
    Point_2  ps1[2];

    if (curve1.contains_point(p))
    {
      ps1[0] = p;
      n1 = 1;
    }
    else
    {
      n1 = curve1.get_points_at_x (p, ps1);
    }

    if (n1 == 0)    // p is not in the x-range of curve1.
      return (EQUAL);
    
    // Get the points on curve2 with the same x co-ordinate as p.
    int    n2;
    Point_2  ps2[2];

    if (curve2.contains_point(p))
    {
      ps2[0] = p;
      n2 = 1;
    }
    else
    {
      n2 = curve2.get_points_at_x (p, ps2);
    }

    if (n2 == 0)    // p is not in the x-range of curve2.
      return (EQUAL);

    // Check that both curves are defined to the left of p.
    // If not, return EQUAL.
    if ( (compare_x(curve1.source(), p) != SMALLER) &&
         (compare_x(curve1.target(), p) != SMALLER) )
      return (EQUAL);

    if ( (compare_x(curve2.source(), p) != SMALLER) &&
         (compare_x(curve2.target(), p) != SMALLER) )
      return (EQUAL);

    // If we reached here, make sure that there is only one point per curve.
    // Only then compare the y co-ordinates of the two points.
    CGAL_assertion(n1 == 1 && n2 == 1);

    Comparison_result result = _compare_y (ps1[0], ps2[0]);

    // In case the two curves do not intersect at the x co-ordinate of p,
    // just return the comparison result at p (since both curves are
    // continuous).    
    if (result != EQUAL)
      return (result);

    // The two curves intersect at ps1[0] = ps2[0] - proceed from here.
    return (_curve_compare_at_intersection_left (curve1, curve2, ps1[0]));
  }  

  // Decide wether curve1 is above, below or equal to curve2 immediately to
  // the right of the x co-ordinate of the given point.
  Comparison_result curve_compare_at_x_right (const X_curve_2& curve1, 
					      const X_curve_2& curve2,
					      const Point_2& p) const
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    // If the two curves are the same, they are equal to the right of p:
    if (curve_is_same(curve1,curve2))
      return (EQUAL);

    // Check the case of vertical segments: we assume that a vertical segment
    // does not contain a point immediately to the right of p.
    if (curve1.is_vertical_segment() || curve2.is_vertical_segment())
      return (EQUAL);

    // Get the points on curve1 with the same x co-ordinate as p.
    int    n1;
    Point_2  ps1[2];

    if (curve1.contains_point(p))
    {
      ps1[0] = p;
      n1 = 1;
    }
    else
    {
      n1 = curve1.get_points_at_x (p, ps1);
    }

    if (n1 == 0)    // p is not in the x-range of curve1.
      return (EQUAL);
    
    // Get the points on curve2 with the same x co-ordinate as p.
    int    n2;
    Point_2  ps2[2];

    if (curve2.contains_point(p))
    {
      ps2[0] = p;
      n2 = 1;
    }
    else
    {
      n2 = curve2.get_points_at_x (p, ps2);
    }

    if (n2 == 0)    // p is not in the x-range of curve2.
      return (EQUAL);

    // Check that both curves are defined to the right of p.
    // If not, return EQUAL.
    if ( (compare_x(curve1.source(), p) != LARGER) &&
         (compare_x(curve1.target(), p) != LARGER) )
      return (EQUAL);

    if ( (compare_x(curve2.source(), p) != LARGER) &&
         (compare_x(curve2.target(), p) != LARGER) )
      return (EQUAL);

    // If we reached here, make sure that there is only one point per curve.
    // Only then compare the y co-ordinates of the two points.
    CGAL_assertion(n1 == 1 && n2 == 1);

    Comparison_result result = _compare_y (ps1[0], ps2[0]);

    // In case the two curves do not intersect at the x co-ordinate of p,
    // just return the comparison result at p (since both curves are
    // continuous).    
    if (result != EQUAL)
      return (result);

    // The two curves intersect at ps1[0] = ps2[0] - proceed from here:
    return (_curve_compare_at_intersection_right (curve1, curve2, ps1[0]));
  }  

  // Check whether the given point is above, under or on the given curve.
  Curve_point_status curve_get_point_status (const X_curve_2& curve,
					     const Point_2& p) const
  {
    CGAL_precondition(is_x_monotone(curve));

    // A special treatment for vertical segments:
    if (curve.is_vertical_segment())
    {
      if (compare_x (curve.source(), p) != EQUAL)
	return (CURVE_NOT_IN_RANGE);

      // In case p has the same x c-ordinate of the vertical segment, compare
      // it to the segment endpoints to determine its position.
      if (_compare_y (curve.source(), p) == SMALLER &&
	  _compare_y (curve.target(), p) == SMALLER)
      {
	return (ABOVE_CURVE);
      }
      else if (_compare_y (curve.source(), p) == LARGER &&
	       _compare_y (curve.target(), p) == LARGER)
      {
	return (UNDER_CURVE);
      }
      else
      {
	return (ON_CURVE);
      }
    }

    // Since the curve is x-monotone, if the point x co-ordinate is to the
    // left (or to the right of both curve's source and target points), then
    // the point is obviously not in the curve's x range.
    if (compare_x(p, curve.source()) == compare_x(p, curve.target()))
    {
      return (CURVE_NOT_IN_RANGE);
    }

    // Check whether the point is exactly on the curve.
    if (curve.contains_point(p))
      return (ON_CURVE);

    // Get the points on the arc with the same x co-ordinate as p.
    int    n;
    Point_2  ps[2];

    n = curve.get_points_at_x (p, ps);

    // Make sure there is exactly one point.
    CGAL_assertion(n == 1);

    // Compare p with the a point of the curve with the same x co-ordinate.
    int result = _compare_y (p, ps[0]);

    if (result == SMALLER)
      return (UNDER_CURVE);
    else if (result == LARGER)
      return (ABOVE_CURVE);

    // We should never reach here:
    CGAL_assertion(false);
    return (ON_CURVE);
  }

  // Check whether the given curve in between c1 and c2, when going in the
  // clockwise direction from p from c1 to c2.
  // Notice that all three curves share the same end-point p.
  bool curve_is_between_cw (const X_curve_2& curve,
			    const X_curve_2& c1, const X_curve_2& c2,
			    const Point_2& p) const
  {
    CGAL_precondition(is_x_monotone(curve));
    CGAL_precondition(is_x_monotone(c1));
    CGAL_precondition(is_x_monotone(c2));

    CGAL_precondition(curve.source().equals(p) || curve.target().equals(p));
    CGAL_precondition(c1.source().equals(p) || c1.target().equals(p));
    CGAL_precondition(c2.source().equals(p) || c2.target().equals(p));

    // Make sure each one of the curves has p as one of its end-points.
    // Let p1_P, p2_P and pv_P point to the other end-points of c1, c2 and cv.
    const Point_2* p1_P = NULL;
    const Point_2* p2_P = NULL;
    const Point_2* pv_P = NULL;
    int            p1_flip = 1, p2_flip = 1, pv_flip = 1;

    if (curve.source().equals(p))
      pv_P = &(curve.target());
    else
    {
      pv_P = &(curve.source());
      pv_flip = -1;
    }

    if (c1.source().equals(p))
      p1_P = &(c1.target());
    else
    {
      p1_P = &(c1.source());
      p1_flip = -1;
    }

    if (c2.source().equals(p))
      p2_P = &(c2.target());
    else
    {
      p2_P = &(c2.source());
      p2_flip = -1;
    }

    // Make sure no two curves overlap.
    if (c1.has_same_base_conic (c2) &&
	(p1_flip * c1.conic().orientation() == 
	 p2_flip * c2.conic().orientation()) &&
	(c1.contains_point(*p2_P) ||
	 c2.contains_point(*p1_P)))
      return (false);

    if (c1.has_same_base_conic(curve) &&
	(p1_flip * c1.conic().orientation() == 
	 pv_flip * curve.conic().orientation()) &&
	(c1.contains_point(*pv_P) ||
	 curve.contains_point(*p1_P)))
      return (false);

    if (curve.has_same_base_conic(c2) &&
	(pv_flip * curve.conic().orientation() == 
	 p2_flip * c2.conic().orientation()) &&
	(curve.contains_point(*p2_P) ||
	 c2.contains_point(*pv_P)))
      return (false);

    // Decide whether each arc is defined to the left or to the right of p.
    bool c1_is_left = (compare_x(*p1_P, p) == SMALLER);
    bool c2_is_left = (compare_x(*p2_P, p) == SMALLER);
    bool curve_is_left = (compare_x(*pv_P, p) == SMALLER);

    // Special treatment for vertical segments:
    // For every curve we define an indicator, which is 0 if it is not a
    // vertical segment, and -1 / 1 if it is a vertical segment going up / down
    // resp.
    int c1_vertical = (c1.is_vertical_segment()) ?
      (_compare_y(*p1_P, p) == LARGER ? 1 : -1) : 0;
    int c2_vertical = (c2.is_vertical_segment()) ?
      (_compare_y(*p2_P, p) == LARGER ? 1 : -1) : 0;
    int curve_vertical = (curve.is_vertical_segment()) ?
      (_compare_y(*pv_P, p) == LARGER ? 1 : -1) : 0;

    if (c1_vertical != 0)
    {
      // c1 is a vertical segment:
      if (c2_vertical != 0)
      {
	if (curve_vertical != 0)
	  return (false);

	// Both c1 and c2 are vertical segments:
	if ((c1_vertical == 1) && (c2_vertical == -1))
	  return (!curve_is_left);
	else if ((c1_vertical == -1) && (c2_vertical == 1))
	  return (curve_is_left);
	else
	  return (false);
      }
      
      if (curve_vertical != 0)
      {
	// Both c1 and curve are vertical segments:
	if ((c1_vertical == 1) && (curve_vertical == -1))
	  return (c2_is_left);
	else if ((c1_vertical == -1) && (curve_vertical == 1))
	  return (!c2_is_left);
	else
	  return (false);
      }
      
      if (c1_vertical == 1)
      {
	// c1 is a vertical segment going up:
	if (c2_is_left)
	{
	  return
            ((!curve_is_left) ||
             (curve_is_left &&
              _curve_compare_at_intersection_left (c2, curve, p) == LARGER));
	}
	else
	{
	  return
            (!curve_is_left &&
             _curve_compare_at_intersection_right (c2, curve, p) == SMALLER);
	}
      }
      else
      {
	// c1 is a vertical segment going down:
	if (c2_is_left)
	{
	  return
            (curve_is_left &&
             _curve_compare_at_intersection_left (c2, curve, p) == LARGER);
	}
	else
	{
	  return
            ((curve_is_left) ||
             (curve_vertical == 1) ||
             (!curve_is_left &&
              _curve_compare_at_intersection_right (c2, curve, p) == SMALLER));
	}
      }
    }

    if (c2_vertical != 0)
    {
      if (curve_vertical != 0)
      {
	// Both c2 and curve are vertical segments:
	if ((c2_vertical == 1) && (curve_vertical == -1))
	  return (!c1_is_left);
	else if ((c2_vertical == -1) && (curve_vertical == 1))
	  return (c1_is_left);
	else
	  return (false);
      }
  
      // Only c2 is a vertical segment:
      if (c2_vertical == 1)
      {
	// c2 is a vertical segment going up:
	if (c1_is_left)
	{
	  return
            (curve_is_left &&
             _curve_compare_at_intersection_left (c1, curve, p) == SMALLER);
	}
	else
	{
	  return
            ((curve_is_left) ||
             (curve_vertical == -1) ||
             (!curve_is_left &&
              _curve_compare_at_intersection_right (c1, curve, p) == LARGER));
	}
      }
      else
      {
	// c2 is a vertical segment going down:
	if (c1_is_left)
	{
	  return
            ((!curve_is_left) ||
             (curve_vertical == 1) ||
             (curve_is_left &&
              _curve_compare_at_intersection_left (c1, curve, p) == SMALLER));
	}
	else
	{
	  return
            ((!curve_is_left) &&
             _curve_compare_at_intersection_right (c1, curve, p) == LARGER);
	}
      }
    }

    if (curve_vertical != 0)
    {
      // Only curve is a vertical segment:
      if (curve_vertical == 1)
      {
	// curve is a vertical segment going up:
	if (c1_is_left && !c2_is_left)
	  return (true);
	else if (c1_is_left && c2_is_left)
	  return (_curve_compare_at_intersection_left(c1, c2, p) == LARGER);
	else if (!c1_is_left && !c2_is_left)
	  return (_curve_compare_at_intersection_right(c1, c2, p) == SMALLER);
      }
      else
      {
	// curve is a vertical segment going down:
	if (!c1_is_left && c2_is_left)
	  return (true);
	else if (c1_is_left && c2_is_left)
	  return (_curve_compare_at_intersection_left(c1, c2, p) == LARGER);
	else if (!c1_is_left && !c2_is_left)
	  return (_curve_compare_at_intersection_right(c1, c2, p) == SMALLER);
      }

      return (false);
    }

    // None of the three curves is a vertical segment:
    // Check the following 4 cases:
    if (c1_is_left && c2_is_left)
    {
      CGAL_assertion(!c1.has_same_base_conic(c2) || 
		     c1.facing() != c2.facing());

      // Case 1: Both c1 and c2 are defined to the left of p.
      if (_curve_compare_at_intersection_left (c1, c2, p) == LARGER)
      { 
	// c1 is above c2:
        return
          (!curve_is_left ||
           !(_curve_compare_at_intersection_left (c2, curve, p) == SMALLER &&
             _curve_compare_at_intersection_left (c1, curve, p) == LARGER));
      }
      else
      { 
	// c2 is above c1:
        return
          (curve_is_left &&
           _curve_compare_at_intersection_left (c1, curve, p) == SMALLER &&
           _curve_compare_at_intersection_left (c2, curve, p) == LARGER);
      }
    }
    else if (!c1_is_left && !c2_is_left) 
    {
      CGAL_assertion(!c1.has_same_base_conic(c2) || 
		     c1.facing() != c2.facing());
      
      // Case 2: Both c1 and c2 are defined to the right of p.
      if (_curve_compare_at_intersection_right (c1, c2, p) == LARGER)
      {
	// c1 is above c2:
        return
          (!curve_is_left &&
           _curve_compare_at_intersection_right (c2, curve, p) == SMALLER &&
           _curve_compare_at_intersection_right (c1, curve, p) == LARGER);
      }
      else
      { 
	// c2 is above c1:
        return
          (curve_is_left ||
           !(_curve_compare_at_intersection_right (c1, curve, p) == SMALLER &&
             _curve_compare_at_intersection_right (c2, curve, p) == LARGER));
      }
    }
    else if (c1_is_left && !c2_is_left)
    {
      // Case 3: c1 is defined to the left and c2 is to the right of p.
      if (curve_is_left)
        return (_curve_compare_at_intersection_left(c1, curve, p) == SMALLER);
      else
        return (_curve_compare_at_intersection_right(c2, curve, p) == SMALLER);
    }
    else if (!c1_is_left && c2_is_left)
    {
      // Case 4: c1 is defined to the right and c2 is to the left of p. 
      if (curve_is_left)
        return (_curve_compare_at_intersection_left (c2, curve, p) == LARGER);
      else
        return (_curve_compare_at_intersection_right(c1, curve, p) == LARGER);
    }

    // We should never reach here.
    CGAL_assertion(false);
    return false;
  }

  // Cehck whether the two points are identical.
  bool point_is_same (const Point_2& p1, const Point_2& p2) const
  {
    return (p1.equals(p2));
  }
   
  // Check whether the two curves are identical.
  bool curve_is_same (const X_curve_2& curve1, const X_curve_2& curve2) const
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    // Check whether all arc features are the same.
    if (curve1.conic().orientation() == curve2.conic().orientation())
    {
      // Same orientation:
      return (curve1.has_same_base_conic(curve2) &&
	      curve1.source().equals(curve2.source()) &&
	      curve1.target().equals(curve2.target()));
    }
    else
    {
      // Check the flip case:
      return (curve1.has_same_base_conic(curve2) &&
	      curve1.source().equals(curve2.target()) &&
	      curve1.target().equals(curve2.source()));
    }
  }

  // Get the source and target vertex of the curve.
  Point_2 curve_source(const X_curve_2& curve) const
  {
    return (curve.source());
  }

  Point_2 curve_target(const X_curve_2& curve) const
  {
    return (curve.target());
  }

  // Return a point to the left or to the right of p.
  Point_2 point_to_left (const Point_2& p) const
  {
    NT   x = CGAL::to_double(p.x()) - 1;
    return (Point_2(x , p.y(),
                    Point_2::User_defined));
  }

  Point_2 point_to_right (const Point_2& p) const
  {
    NT   x = CGAL::to_double(p.x()) + 1;
    return (Point_2(x , p.y(),
                    Point_2::User_defined));
  }

#if defined(HAS_LEFT_NOT)
  // Reflect a point in y.
  Point_2 point_reflect_in_y (const Point_2& p) const
  {
    return (p.reflect_in_y());
  }
      
  // Reflect a curve in y.
  X_curve_2 curve_reflect_in_y (const X_curve_2& curve) const
  {
    return (curve.reflect_in_y());
  }

  // Reflect a point in x and y.
  Point_2 point_reflect_in_x_and_y (const Point_2& p) const
  {
    return (p.reflect_in_x_and_y());
  }
      
  // Reflect a curve in x and y.
  X_curve_2 curve_reflect_in_x_and_y (const X_curve_2& curve) const
  {
    return (curve.reflect_in_x_and_y());
  }
#endif
    
  ////////// Arrangement methods: //////////

  // Change the orientation of the curve (swap the source and the target).
  X_curve_2 curve_flip (const X_curve_2& curve) const
  {
    // Flip the arc.
    return (curve.flip());
  }

  // Check whether the curve is x-monotone.
  bool is_x_monotone (const Curve_2& curve) const
  {
    return (curve.is_x_monotone());
  }

  // Cut the curve to several x-monotone sub-curves.
  void make_x_monotone (const Curve_2& curve, 
			std::list<X_curve_2>& x_curves) const
  {
    CGAL_precondition(!is_x_monotone(curve));

    // Clear the output list.
    x_curves.clear();

    // Find the points of vertical tangency and act accordingly.
    int    n;
    Point_2  ps[2];

    n = curve.vertical_tangency_points (ps);

    CGAL_assertion (n > 0);

    // Split the conic arc into x-monotone sub-curves. 
    if (curve.is_full_conic())
    {
      // Make sure we have two vertical tangency points.
      CGAL_assertion(n == 2);

      // In case the curve is a full conic, split it to two x-monotone curves,
      // one going from ps[0] to ps[1], and the other from ps[1] to ps[0].
      x_curves.push_back (X_curve_2 (curve, ps[0], ps[1], false));
      x_curves.push_back (X_curve_2 (curve, ps[1], ps[0], false));
    }
    else
    {
      X_curve_2    sub_curve1;
      X_curve_2    sub_curve2;
      X_curve_2    sub_curve3;

      if (n == 1)
      {
	// Split the arc into two x-monotone sub-curves: one going from the
	// arc source to ps[0], and the other from ps[0] to the target. 
	_curve_split (curve, 
	              sub_curve1, sub_curve2, 
                      ps[0]);

	x_curves.push_back(sub_curve1);
	x_curves.push_back(sub_curve2);
      }
      else if (n == 2)
      {
	// Split the arc into three x-monotone sub-curves: one going from the
	// arc source to ps[0], one from ps[0] to ps[1], and the last one
	// from ps[1] to the target.
	// Notice that ps[0] and ps[1] might switch places.
	X_curve_2    temp;

	_curve_split (curve, 
	              sub_curve1, sub_curve2, 
                      ps[0]);

	if (sub_curve2.contains_point(ps[1]))
	{
	  temp = sub_curve2;
	  _curve_split (temp,
		        sub_curve2, sub_curve3,
		        ps[1]);
	}
	else if (sub_curve1.contains_point(ps[1]))
	{
	  // Actually we switch between ps[0] and ps[1].
	  temp = sub_curve1;
	  sub_curve3 = sub_curve2;
	  _curve_split (temp,
		        sub_curve1, sub_curve2,
		        ps[1]);
	}
	else
	{
	  // We should never reach here:
	  CGAL_assertion(false);
	}

	x_curves.push_back(sub_curve1);
	x_curves.push_back(sub_curve2);
	x_curves.push_back(sub_curve3);		    	
      }
      else
      {
	// We should never reach here:
	CGAL_assertion(false);
      }
    }

    return;
  }

  // Split the given curve into two sub-curves at the given point.
  void curve_split (const X_curve_2& curve, 
		    X_curve_2& sub_curve1, X_curve_2& sub_curve2, 
                    const Point_2& p) const 
  {
    CGAL_precondition(is_x_monotone(curve));

    // Make sure the point is on the curve and is not an end-point.
    CGAL_precondition(curve.contains_point(p));
    CGAL_precondition(! p.equals(curve.source()));
    CGAL_precondition(! p.equals(curve.target()));

    // Split the curve.
    _curve_split (curve,
		  sub_curve1, sub_curve2,
		  p);
    return;
  }

  // Check whether the intersection point between the two given curves is
  // lexicographically strictly to right of the given point.
  bool do_intersect_to_right(const X_curve_2& curve1, const X_curve_2& curve2,
                             const Point_2& p) const 
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    // Deal with overlapping curves:
    int     n_ovlps;
    X_curve_2 ovlp_arcs[2];

    n_ovlps = curve1.overlaps (curve2, ovlp_arcs);
    CGAL_assertion (n_ovlps < 2);

    if (n_ovlps == 1)
    {
      // Check if at least one end-point of the overlapping curve is to the
      // right of p.
      return (ovlp_arcs[0].source().compare_lex_xy(p) == LARGER ||
	      ovlp_arcs[0].target().compare_lex_xy(p) == LARGER);
    }

    // In case there are no overlaps and the base conics are the same,
    // there cannot be any intersection points, unless the two x-monotone
    // curves share an end point.
    if (curve1.has_same_base_conic(curve2))
    {
      if ((curve1.source().equals(curve2.source()) ||
	   curve1.source().equals(curve2.target())) &&
	  curve1.source().compare_lex_xy(p) == LARGER)
      {
	return (true);
      }

      if ((curve1.target().equals(curve2.source()) ||
	   curve1.target().equals(curve2.target())) &&
	  curve1.target().compare_lex_xy(p) == LARGER)
      {
	return (true);
      }
      
      // No overlaps at all: the two curves do not intersect.
      return (false);
    }

    // Find the intersection points and decide accordingly.
    int   n;
    Point_2 ps[4];
  
#ifdef CGAL_CONIC_ARC_USE_CACHING
    n = curve1.intersections_with (curve2, ps, &inter_list);
#else
    n = curve1.intersections_with (curve2, ps);    
#endif

    for (int i = 0; i < n; i++)
    {
      if (ps[i].compare_lex_xy(p) == LARGER)
	return (true);
    }

    return (false);
  }

  // RWRW: {
  // THE NEXT FUNCTION IS NOT REALLY USED AND WILL BE OVERRIDEN BY THE WRAPPER.
  // Check whether the intersection point between the two given curves is
  // lexicographically strictly to left of the given point.
  bool do_intersect_to_left (const X_curve_2& curve1, const X_curve_2& curve2,
                             const Point_2& p) const 
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    // Deal with overlapping curves:
    int     n_ovlps;
    X_curve_2 ovlp_arcs[2];

    n_ovlps = curve1.overlaps (curve2, ovlp_arcs);
    CGAL_assertion (n_ovlps < 2);

    if (n_ovlps == 1)
    {
      // Check if at least one end-point of the overlapping curve is to the
      // left of p.
      return (ovlp_arcs[0].source().compare_lex_xy(p) == SMALLER ||
	      ovlp_arcs[0].target().compare_lex_xy(p) == SMALLER);
    }

    // In case there are no overlaps and the base conics are the same,
    // there cannot be any intersection points, unless the two x-monotone
    // curves share an end point.
    if (curve1.has_same_base_conic(curve2))
    {
      if ((curve1.source().equals(curve2.source()) ||
	   curve1.source().equals(curve2.target())) &&
	  curve1.source().compare_lex_xy(p) == SMALLER)
      {
	return (true);
      }

      if ((curve1.target().equals(curve2.source()) ||
	   curve1.target().equals(curve2.target())) &&
	  curve1.target().compare_lex_xy(p) == SMALLER)
      {
	return (true);
      }
      
      // No overlaps at all: the two curves do not intersect.
      return (false);
    }

    // Find the intersection points and decide accordingly.
    int   n;
    Point_2 ps[4];
  
#ifdef CGAL_CONIC_ARC_USE_CACHING
    n = curve1.intersections_with (curve2, ps, &inter_list);
#else
    n = curve1.intersections_with (curve2, ps);    
#endif

    for (int i = 0; i < n; i++)
    {
      if (ps[i].compare_lex_xy(p) == SMALLER)
	return (true);
    }

    return (false);
  }
  // RWRW: }

  // Find the nearest intersection point between the two given curves to the
  // right of the given point.
  // In case of an overlap, p1 and p2 are the source and destination of the
  // overlapping curve. Otherwise p1=p2 is the calculated intersection point.
  bool nearest_intersection_to_right (const X_curve_2& curve1,
				      const X_curve_2& curve2,
				      const Point_2& p,
                                      Point_2& p1,
                                      Point_2& p2) const
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    // Deal with overlapping curves:
    int     n_ovlps;
    X_curve_2 ovlp_arcs[2];

    n_ovlps = curve1.overlaps (curve2, ovlp_arcs);
    CGAL_assertion (n_ovlps < 2);

    if (n_ovlps == 1)
    {
      Point_2  ovlp_source = ovlp_arcs[0].source();
      Point_2  ovlp_target = ovlp_arcs[0].target();

      if (ovlp_source.compare_lex_xy(p) == LARGER &&
	  ovlp_target.compare_lex_xy(p) == LARGER)
      {
	// The entire overlapping arc is to the right of p:
	p1 = ovlp_source;
	p2 = ovlp_target;
	return (true);
      }
      else if (ovlp_source.compare_lex_xy(p) != LARGER &&
	       ovlp_target.compare_lex_xy(p) == LARGER)
      {
	// The source is to the left of p, and the traget is to its right.
	p1 = p;
	p2 = ovlp_target;
	return (true);
      }
      else if (ovlp_source.compare_lex_xy(p) == LARGER &&
	       ovlp_target.compare_lex_xy(p) != LARGER)
      {
	// The source is to the right of p, and the traget is to its left.
	p1 = ovlp_source;
	p2 = p;
	return (true);
      }
      else
      {
	// The entire overlapping arc is to the left of p:
	return (false);
      }
    }

    // In case there are no overlaps and the base conics are the same,
    // there cannot be any intersection points, unless the two x-monotone
    // curves share an end point.
    if (curve1.has_same_base_conic(curve2))
    {
      const Point_2    *nearest_end_P = NULL;

      if ((curve1.source().equals(curve2.source()) ||
	   curve1.source().equals(curve2.target())) &&
	  curve1.source().compare_lex_xy(p) == LARGER)
      {
	nearest_end_P = &(curve1.source());
      }

      if ((curve1.target().equals(curve2.source()) ||
	   curve1.target().equals(curve2.target())) &&
	  curve1.target().compare_lex_xy(p) == LARGER)
      {
	if (nearest_end_P == NULL ||
	    nearest_end_P->compare_lex_xy (curve1.target()) == LARGER)
	{
	  nearest_end_P = &(curve1.target());
	}
      }

      if (nearest_end_P != NULL)
      {
	// A common end point was found:
	p1 = p2 = *nearest_end_P;
	return (true);
      }
      else
      {
	// No intersection:
	return (false);
      }
    }

    // Find the intersection points and choose the one nearest to p.
    int          n;
    Point_2        ps[4];
    const Point_2 *nearest_inter_P = NULL;
  
#ifdef CGAL_CONIC_ARC_USE_CACHING
    n = curve1.intersections_with (curve2, ps, &inter_list);
#else
    n = curve1.intersections_with (curve2, ps);    
#endif

    for (int i = 0; i < n; i++)
    {
      // Check if the current point is to the right of p.
      if (ps[i].compare_lex_xy(p) == LARGER)
      {
	  // Compare with the nearest point so far.
	if (nearest_inter_P == NULL ||
	    nearest_inter_P->compare_lex_xy (ps[i]) == LARGER)	
	{
	    nearest_inter_P = &(ps[i]);
	}
      }
    }

    if (nearest_inter_P != NULL)
    {
      // Return the nearest intersection point.
      p1 = p2 = *nearest_inter_P;
      return (true);
    }
    
    // No intersection found.
    return (false);
  }

  // RWRW: {
  // THE NEXT FUNCTION IS NOT REALLY USED AND WILL BE OVERRIDEN BY THE WRAPPER.
  // Find the nearest intersection point between the two given curves to the
  // left of the given point.
  // In case of an overlap, p1 and p2 are the source and destination of the
  // overlapping curve. Otherwise p1=p2 is the calculated intersection point.
  bool nearest_intersection_to_left (const X_curve_2& curve1,
				     const X_curve_2& curve2,
				     const Point_2& p,
				     Point_2& p1,
				     Point_2& p2) const
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    // Deal with overlapping curves:
    int     n_ovlps;
    X_curve_2 ovlp_arcs[2];

    n_ovlps = curve1.overlaps (curve2, ovlp_arcs);
    CGAL_assertion (n_ovlps < 2);

    if (n_ovlps == 1)
    {
      Point_2  ovlp_source = ovlp_arcs[0].source();
      Point_2  ovlp_target = ovlp_arcs[0].target();

      if (ovlp_source.compare_lex_xy(p) == SMALLER &&
	  ovlp_target.compare_lex_xy(p) == SMALLER)
      {
	// The entire overlapping arc is to the left of p:
	p1 = ovlp_source;
	p2 = ovlp_target;
	return (true);
      }
      else if (ovlp_source.compare_lex_xy(p) != SMALLER &&
	       ovlp_target.compare_lex_xy(p) == SMALLER)
      {
	// The source is to the right of p, and the traget is to its left.
	p1 = p;
	p2 = ovlp_target;
	return (true);
      }
      else if (ovlp_source.compare_lex_xy(p) == SMALLER &&
	       ovlp_target.compare_lex_xy(p) != SMALLER)
      {
	// The source is to the left of p, and the traget is to its right.
	p1 = ovlp_source;
	p2 = p;
	return (true);
      }
      else
      {
	// The entire overlapping arc is to the right of p:
	return (false);
      }
    }

    // In case there are no overlaps and the base conics are the same,
    // there cannot be any intersection points, unless the two x-monotone
    // curves share an end point.
    if (curve1.has_same_base_conic(curve2))
    {
      const Point_2    *nearest_end_P = NULL;

      if ((curve1.source().equals(curve2.source()) ||
	   curve1.source().equals(curve2.target())) &&
	  curve1.source().compare_lex_xy(p) == SMALLER)
      {
	nearest_end_P = &(curve1.source());
      }

      if ((curve1.target().equals(curve2.source()) ||
	   curve1.target().equals(curve2.target())) &&
	  curve1.target().compare_lex_xy(p) == SMALLER)
      {
	if (nearest_end_P == NULL ||
	    nearest_end_P->compare_lex_xy (curve1.target()) == SMALLER)
	{
	  nearest_end_P = &(curve1.target());
	}
      }

      if (nearest_end_P != NULL)
      {
	// A common end point was found:
	p1 = p2 = *nearest_end_P;
	return (true);
      }
      else
      {
	// No intersection:
	return (false);
      }
    }

    // Find the intersection points and choose the one nearest to p.
    int          n;
    Point_2        ps[4];
    const Point_2 *nearest_inter_P = NULL;
  
#ifdef CGAL_CONIC_ARC_USE_CACHING
    n = curve1.intersections_with (curve2, ps, &inter_list);
#else
    n = curve1.intersections_with (curve2, ps);    
#endif

    for (int i = 0; i < n; i++)
    {
      // Check if the current point is to the right of p.
      if (ps[i].compare_lex_xy(p) == SMALLER)
      {
	  // Compare with the nearest point so far.
	if (nearest_inter_P == NULL ||
	    nearest_inter_P->compare_lex_xy (ps[i]) == SMALLER)	
	{
	    nearest_inter_P = &(ps[i]);
	}
      }
    }

    if (nearest_inter_P != NULL)
    {
      // Return the nearest intersection point.
      p1 = p2 = *nearest_inter_P;
      return (true);
    }
    
    // No intersection found.
    return (false);
  }
  // RWRW: }

  // Check whether two curves overlap.
  bool curves_overlap (const X_curve_2& curve1, const X_curve_2& curve2) const
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    X_curve_2 ovlp_arcs[2];

    return (curve1.overlaps (curve2, ovlp_arcs) > 0);
  }

 private:

  ////////// Private auxiliary methods: //////////

  // Compare two points by their y coordinate:
  Comparison_result _compare_y(const Point_2& p0, const Point_2& p1) const
  {
    return (p0.compare_y(p1));
  }

  // Split the given curve into two sub-curves at the given point.
  // Since this is a private function, there are no preconditions.
  void _curve_split (const X_curve_2& curve, 
		     X_curve_2& sub_curve1, X_curve_2& sub_curve2, 
		     const Point_2& p) const
  {
    CGAL_precondition(! p.equals(curve.source()));
    CGAL_precondition(! p.equals(curve.target()));

    // Split the curve to source->p and p->target.
    sub_curve1 = X_curve_2 (curve, curve.source(), p, false);
    sub_curve2 = X_curve_2 (curve, p, curve.target(), false);
    
    CGAL_assertion(sub_curve1.is_x_monotone());
    CGAL_assertion(sub_curve2.is_x_monotone());

    return;
  }

  // Decide wether curve1 is above, below or equal to curve2 immediately to
  // the left of p_int, which is assumed to be an intersection of both curves.
  // Furthermore, the two curves are assumed to be defined to p_int's left.
  Comparison_result _curve_compare_at_intersection_left 
                                            (const X_curve_2& curve1, 
					     const X_curve_2& curve2,
					     const Point_2& p_int) const
  {
    // In case one arc is facing upwards and another facing downwards, it is
    // very clear that the one facing upward is above the one facing downwards.
    if (curve1.has_same_base_conic(curve2))
    {
      Comparison_result fac1 = curve1.facing();
      Comparison_result fac2 = curve2.facing();

      if (fac1 == LARGER && fac2 == SMALLER)
	return (LARGER);
      else if (fac1 == SMALLER && fac2 == LARGER)
	return (SMALLER);
    }

    // Otherwise, the two curves do intersect at p_int = ps1[0] = ps2[0]:
    // make a decision based on their partial derivatives.
    //const Point_2& p_int = ps1[0];
    const NT       _zero = 0;

    NT    slope1_numer, slope1_denom;
    NT    slope2_numer, slope2_denom;

    curve1.derive_by_x_at (p_int, 1, slope1_numer, slope1_denom);
    curve2.derive_by_x_at (p_int, 1, slope2_numer, slope2_denom);

    /*const*/ bool  is_vertical_slope1 = (slope1_denom == _zero);
    /*const*/ bool  is_vertical_slope2 = (slope2_denom == _zero);

    // RWRW - CHECK THIS !!!
    if (! is_vertical_slope1 &&
	p_int.is_approximate() &&
	eps_compare<APNT>(TO_APNT(slope1_denom), 0) == EQUAL)
    {
      is_vertical_slope1 = true;
    }
    if (! is_vertical_slope2 &&
	p_int.is_approximate() &&
	eps_compare<APNT>(TO_APNT(slope2_denom), 0) == EQUAL)
    {
      is_vertical_slope2 = true;
    }
    // TO HERE ...

    if (!is_vertical_slope1 && !is_vertical_slope2)
    {
      // The two curves have derivatives at p_int: use it to determine which
      // one is above the other (the one with a smaller slope in above). 
      
      // RWRW - CHECK THIS !!!
      Comparison_result  slope_res = 
	(!p_int.is_approximate()) ?
	CGAL::compare(slope2_numer*slope1_denom, slope1_numer*slope2_denom) :
	eps_compare<APNT>(TO_APNT(slope2_numer*slope1_denom), 
			  TO_APNT(slope1_numer*slope2_denom));
            
      if (slope_res != EQUAL)
	return (slope_res);

      // Use the second order derivative.
      curve1.derive_by_x_at (p_int, 2, slope1_numer, slope1_denom);
      curve2.derive_by_x_at (p_int, 2, slope2_numer, slope2_denom);

      slope_res = CGAL::compare (slope2_numer*slope1_denom, 
				 slope1_numer*slope2_denom);

      CGAL_assertion(slope_res != EQUAL);
	
      return ((slope_res == LARGER) ? SMALLER : LARGER);
    }
    else if (!is_vertical_slope2)
    {
      // The first curve has a vertical slope at p_int: check whether it is
      // facing upwards or downwards and decide accordingly.
      Comparison_result fac1 = curve1.facing();

      CGAL_assertion(fac1 != EQUAL);
      return (fac1);
    }
    else if (!is_vertical_slope1)
    {
      // The second curve has a vertical slope at p_int: check whether it is
      // facing upwards or downwards and decide accordingly.
      Comparison_result fac2 = curve2.facing();

      CGAL_assertion(fac2 != EQUAL);
      return ((fac2 == LARGER) ? SMALLER : LARGER);
    }
    else
    {
      // The two curves have vertical slopes at p_int: 
      // First check whether one is facing up and one down. In this case the
      // comparison result is trivial.
      Comparison_result fac1 = curve1.facing();
      Comparison_result fac2 = curve2.facing();

      if (fac1 == LARGER && fac2 == SMALLER)
	return (LARGER);
      else if (fac1 == SMALLER && fac2 == LARGER)
	return (SMALLER);

      // Compute the second order derivative by y and act according to it.
      curve1.derive_by_y_at (p_int, 2, slope1_numer, slope1_denom);
      curve2.derive_by_y_at (p_int, 2, slope2_numer, slope2_denom);

      Comparison_result  slope_res = CGAL::compare(slope2_numer*slope1_denom, 
						   slope1_numer*slope2_denom);

      CGAL_assertion(slope_res != EQUAL);

      if (fac1 == LARGER && fac2 == LARGER)
      {
	// Both are facing up.
	return ((slope_res == LARGER) ? SMALLER : LARGER);
      }
      else
      {
	// Both are facing down.
	return (slope_res);
      }
    }

    // We should never reach here:
    CGAL_assertion(false);
    return (EQUAL);
  }

  // Decide wether curve1 is above, below or equal to curve2 immediately to
  // the right of p_int, which is assumed to be an intersection of both curves.
  // Furthermore, the two curves are assumed to be defined to p_int's right.
  Comparison_result _curve_compare_at_intersection_right 
                                            (const X_curve_2& curve1, 
					     const X_curve_2& curve2,
					     const Point_2& p_int) const
  {
    // In case one arc is facing upwards and another facing downwards, it is
    // very clear that the one facing upward is above the one facing downwards.
    if (curve1.has_same_base_conic(curve2))
    {
      Comparison_result fac1 = curve1.facing();
      Comparison_result fac2 = curve2.facing();

      if (fac1 == LARGER && fac2 == SMALLER)
	return (LARGER);
      else if (fac1 == SMALLER && fac2 == LARGER)
	return (SMALLER);
    }

    // Otherwise, the two curves do intersect at p_int = ps1[0] = ps2[0]:
    // make a decision based on their partial derivatives.
    //const Point_2& p_int = ps1[0];
    const NT       _zero = 0;

    NT    slope1_numer, slope1_denom;
    NT    slope2_numer, slope2_denom;

    curve1.derive_by_x_at (p_int, 1, slope1_numer, slope1_denom);
    curve2.derive_by_x_at (p_int, 1, slope2_numer, slope2_denom);

    /*const*/ bool  is_vertical_slope1 = (slope1_denom == _zero);
    /*const*/ bool  is_vertical_slope2 = (slope2_denom == _zero);

    // RWRW - CHECK THIS !!!
    if (! is_vertical_slope1 &&
	p_int.is_approximate() &&
	eps_compare<APNT>(TO_APNT(slope1_denom), 0) == EQUAL)
    {
      is_vertical_slope1 = true;
    }
    if (! is_vertical_slope2 && 
	p_int.is_approximate() &&
	eps_compare<APNT>(TO_APNT(slope2_denom), 0) == EQUAL)
    {
      is_vertical_slope2 = true;
    }
    // TO HERE ...
 
    if (!is_vertical_slope1 && !is_vertical_slope2)
    {
      // The two curves have derivatives at p_int: use it to determine which
      // one is above the other (the one with a larger slope is below).
      // RWRW - CHECK THIS !!!
      Comparison_result  slope_res =
	(! p_int.is_approximate()) ?
	CGAL::compare (slope1_numer*slope2_denom, slope2_numer*slope1_denom) :
	eps_compare<APNT> (TO_APNT(slope1_numer*slope2_denom), 
			   TO_APNT(slope2_numer*slope1_denom));
      
      if (slope_res != EQUAL)
	return (slope_res);

      // Use the second order derivative.
      curve1.derive_by_x_at (p_int, 2, slope1_numer, slope1_denom);
      curve2.derive_by_x_at (p_int, 2, slope2_numer, slope2_denom);

      slope_res = CGAL::compare (slope1_numer*slope2_denom, 
				 slope2_numer*slope1_denom);

      CGAL_assertion(slope_res != EQUAL);
	
      return (slope_res);
    }
    else if (!is_vertical_slope2)
    {
      // The first curve has a vertical slope at p_int: check whether it is
      // facing upwards or downwards and decide accordingly.
      Comparison_result fac1 = curve1.facing();

      CGAL_assertion(fac1 != EQUAL);
      return (fac1);
    }
    else if (!is_vertical_slope1)
    {
      // The second curve has a vertical slope at p_int: check whether it is
      // facing upwards or downwards and decide accordingly.
      Comparison_result fac2 = curve2.facing();

      CGAL_assertion(fac2 != EQUAL);
      return ((fac2 == LARGER) ? SMALLER : LARGER);
    }
    else
    {
      // The two curves have vertical slopes at p_int: 
      // First check whether one is facing up and one down. In this case the
      // comparison result is trivial.
      Comparison_result fac1 = curve1.facing();
      Comparison_result fac2 = curve2.facing();

      if (fac1 == LARGER && fac2 == SMALLER)
	return (LARGER);
      else if (fac1 == SMALLER && fac2 == LARGER)
	return (SMALLER);

      // Compute the second order derivative by y and act according to it.
      curve1.derive_by_y_at (p_int, 2, slope1_numer, slope1_denom);
      curve2.derive_by_y_at (p_int, 2, slope2_numer, slope2_denom);

      Comparison_result  slope_res = CGAL::compare (slope1_numer*slope2_denom, 
						    slope2_numer*slope1_denom);

      CGAL_assertion(slope_res != EQUAL);

      if (fac1 == LARGER && fac2 == LARGER)
      {
	// Both are facing up.
	return ((slope_res == LARGER) ? SMALLER : LARGER);
      }
      else
      {
	// Both are facing down.
	return (slope_res);
      }
    }

    // We should never reach here:
    CGAL_assertion(false);
    return (EQUAL);
  }

};

CGAL_END_NAMESPACE

#endif
