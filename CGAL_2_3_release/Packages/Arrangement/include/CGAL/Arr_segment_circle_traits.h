// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-79 $
// release_date  : $CGAL_Date: 2001/07/03 $
//
// file          : include/CGAL/Arr_segment_circle_traits.h
// package       : Arrangement 
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Ron Wein <wein@post.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_ARR_CONIC_TRAITS_H
#define CGAL_ARR_CONIC_TRAITS_H

#include <CGAL/basic.h>
//#include <CGAL/Segment_2.h>
//#include <CGAL/Circle_2.h>
#include <list>

#include <CGAL/Segment_circle_2.h>

CGAL_BEGIN_NAMESPACE

// ----------------------------------------------------------------------------
// Arrangement traits for conic arcs.
//

template <class _NT>
class Arr_segment_circle_traits 
{
 public:
  
  typedef _NT                    NT;

  // The difference between Curve and X_curve is semantical only,
  // NOT syntactical.
  typedef Segment_circle_2<NT>   Curve;
  typedef Curve                  X_curve;
 
  // Using typename to please compiler (e.g., CC with IRIX64 on mips)
  typedef typename Curve::R       R;
  typedef typename Curve::Point   Point;
  typedef typename Curve::Segment Segment;
  typedef typename Curve::Circle  Circle;
  typedef typename Curve::Conic   Conic;

#ifndef __GNUC__
  enum Curve_point_status {
    UNDER_CURVE = -1,
    ABOVE_CURVE = 1,
    ON_CURVE = 2,
    CURVE_NOT_IN_RANGE = 0,
  };
#else
  //workaround for egcs, otherwise we get an ICE
  typedef int Curve_point_status;
  static const int UNDER_CURVE = -1;
  static const int ABOVE_CURVE = 1;
  static const int ON_CURVE = 2;
  static const int CURVE_NOT_IN_RANGE = 0;
#endif

  // Constructor.
  Arr_segment_circle_traits()
    {}

  ////////// Planar Map methods: //////////

  // Compare the co-ordinates of two given points.
  Comparison_result compare_x(const Point& p0, const Point& p1) const
  {
    return _compare_value(p0.x(),p1.x());
  }

  Comparison_result compare_y(const Point& p0, const Point& p1) const
  {
    return _compare_value(p0.y(),p1.y());
  }

  // Check whether the given curve is a vertical segment.
  bool curve_is_vertical(const X_curve& curve) const
  {
    return (curve.is_vertical_segment());
  }

  // Check whether the x co-ordinate of the given point is contained in the
  // x-range of the given x-monotone curve.
  bool curve_is_in_x_range (const X_curve& curve, const Point& p) const
  {
    CGAL_precondition(is_x_monotone(curve));

    if (curve.is_vertical_segment())
    {
      return (compare_x (curve.source(), p) == EQUAL);
    }
    else
    {
      // Find the number of points on the arc with the same x co-ordinate as p.
      // p is is the x-range of the arc only if there is at least one point.
      Point  ps[2];
      int    n = curve.get_points_at_x (p.x(), ps);

      return (n > 0);
    }
  }

  // Decide wether curve1 is above, below or equal to curve2 at the
  // x co-ordinate of the given point.
  Comparison_result curve_compare_at_x (const X_curve& curve1, 
				        const X_curve& curve2, 
				        const Point& p) const
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    // Get the points on curve1 with the same x co-ordinate as p.
    int    n1;
    Point  ps1[2];

    if (curve1.is_vertical_segment())
    {
      // Vertical segment.
      if (compare_x (curve1.source(), p) != EQUAL)
	return (EQUAL);

      n1 = 2;
      if (compare_y (curve1.source(), curve1.target()) == SMALLER)
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
    else
    {
      n1 = curve1.get_points_at_x (p.x(), ps1);

      if (n1 == 0)    // p is not in the x-range of curve1.
	return (EQUAL);
    }

    // Get the points on curve2 with the same x co-ordinate as p.
    int    n2;
    Point  ps2[2];

    if (curve2.is_vertical_segment())
    {
      // Vertical segment.
      if (compare_x (curve2.source(), p) != EQUAL)
	return (EQUAL);
      
      n2 = 2;
      if (compare_y (curve2.source(), curve2.target()) == SMALLER)
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
    else
    {
      n2 = curve2.get_points_at_x (p.x(), ps2);

      if (n2 == 0)    // p is not in the x-range of curve2.
	return (EQUAL);
    }

    // Deal with vertical segments:
    if (n1 == 2)
    {
      // Check if the vertical segment contains ps2[0].
      if (compare_y (ps1[0], ps2[0]) != LARGER && 
	  compare_y (ps1[1], ps2[0]) != SMALLER)
      {
	return (EQUAL);
      }

      if (n2 == 2)
      {
	if (compare_y (ps1[0], ps2[1]) != LARGER && 
	    compare_y (ps1[1], ps2[1]) != SMALLER)
	{
	  return (EQUAL);
	}
      }
    }
    else if (n2 == 2)
    {
      // Check if the vertical segment contains ps1[0].
      if (compare_y (ps2[0], ps1[0]) != LARGER && 
	  compare_y (ps2[1], ps1[0]) != SMALLER)
      {
	return (EQUAL);
      }
    }

    // Compare the y co-ordinates of the two points.
    return (compare_y (ps1[0], ps2[0]));
  }

  // Decide wether curve1 is above, below or equal to curve2 immediately to
  // the left of the x co-ordinate of the given point.
  Comparison_result curve_compare_at_x_left (const X_curve& curve1, 
					     const X_curve& curve2,
					     const Point& p) const
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
    Point  ps1[2];

    n1 = curve1.get_points_at_x (p.x(), ps1);

    if (n1 == 0)    // p is not in the x-range of curve1.
      return (EQUAL);
    
    // Get the points on curve2 with the same x co-ordinate as p.
    int    n2;
    Point  ps2[2];

    n2 = curve2.get_points_at_x (p.x(), ps2);

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

    Comparison_result result = compare_y (ps1[0], ps2[0]);

    // In case the two curves do not intersect at the x co-ordinate of p,
    // just return the comparison result at p (since both curves are
    // continuous).    
    if (result != EQUAL)
      return (result);

    // Otherwise, the two curves do intersect at p_int = ps1[0] = ps2[0]:
    // make a decision based on their partial derivatives.
    const Point& p_int = ps1[0];
    
    // In order to simplify the process, make sure the source is always to the
    // left of p.
    X_curve c1 = curve1;
    X_curve c2 = curve2;
    
    if (compare_x(c1.source(), p) != SMALLER)
      c1 = curve1.flip();
    if (compare_x(c2.source(), p) != SMALLER)
      c2 = curve2.flip();

    // Calculate the partial derivatives at the intersection point.
    NT       dC1_dx, dC1_dy;
    NT       dC2_dx, dC2_dy;
    const NT zero = 0;
    
    c1.partial_derivatives (p_int, dC1_dx, dC1_dy);
    c2.partial_derivatives (p_int, dC2_dx, dC2_dy);

    // Try to compare the slopes of the tangents to the two arcs.
    if (dC1_dy != zero && dC2_dy != zero)
    {
      Comparison_result  slope_result = _compare_value (dC1_dx/dC1_dy, 
      							dC2_dx/dC2_dy);
      if (slope_result != EQUAL)
      {
	return (slope_result);
      }
    }

    // If we reached here, the two curves intersect at p_int, and their
    // tangents at this point have the same slope. We shall therefore find
    // a point with extreme y-coordinate for each curve and compare the
    // y-values of these two extreme points.
    Point  p_extr1, p_extr2;
    Point  hpts[2];
    int    n_hpts;

    n_hpts = c1.horizontal_tangency_points(hpts);
    if (n_hpts == 0)
    {
      // No horizontal tangency points - the source must have an extreme
      // y-value.
      p_extr1 = c1.source();
    }
    else if (n_hpts == 1)
    {
      // Use the vertical tangency points only if it lies between the source
      // and p_int. Otherwise, use the source.
      if (compare_x(c1.source(), hpts[0]) == SMALLER &&
	  compare_x(p_int, hpts[0]) == LARGER)
      {
	p_extr1 = hpts[0];
      }
      else
      {
	p_extr1 = c1.source();
      }
    }
    else
    {
      // An x-monotone curve cannot have more than 1 horizontal tangencies.
      CGAL_assertion (false);
    }

    n_hpts = c2.horizontal_tangency_points(hpts);
    if (n_hpts == 0)
    {
      // No horizontal tangency points - the source must have an extreme
      // y-value.
      p_extr2 = c2.source();
    }
    else if (n_hpts == 1)
    {
      // Use the vertical tangency points only if it lies between the source
      // and p_int. Otherwise, use the source.
      if (compare_x(c2.source(), hpts[0]) == SMALLER &&
	  compare_x(p_int, hpts[0]) == LARGER)
      {
	p_extr2 = hpts[0];
      }
      else
      {
	p_extr2 = c2.source();
      }
    }
    else
    {
      // An x-monotone curve cannot have more than 1 horizontal tangencies.
      CGAL_assertion (false);
    }

    // Use the x-coordinate of the point next to x.
    result = compare_x (p_extr1, p_extr2);

    if (result == SMALLER)
    {
      n_hpts = c1.get_points_at_x (p_extr2.x(), hpts);
      CGAL_assertion(n_hpts == 1);
      p_extr1 = hpts[0];
    }
    else if (result == LARGER)
    {
      n_hpts = c2.get_points_at_x (p_extr1.x(), hpts);
      CGAL_assertion(n_hpts == 1);
      p_extr2 = hpts[0];
    }

    // Compare the y-coordinates of the two extreme points.
    Comparison_result extr_result = compare_y (p_extr1, p_extr2);

    if (extr_result != EQUAL)
      return (extr_result);

    // As a last resort:
    NT    x_mid = (p_extr1.x() + p_int.x()) / NT(2);
    Point p_mid1, p_mid2;

    n_hpts = c1.get_points_at_x (x_mid, hpts);
    CGAL_assertion(n_hpts == 1);
    p_mid1 = hpts[0];

    n_hpts = c2.get_points_at_x (x_mid, hpts);
    CGAL_assertion(n_hpts == 1);
    p_mid2 = hpts[0];

    Comparison_result mid_result = compare_y (p_mid1, p_mid2);

    if (mid_result != EQUAL)
      return (mid_result);

    // If we reached here, the two curves must be overlapping:
    int     n_ovlps;
    X_curve ovlp_arcs[2];

    n_ovlps = curve1.overlaps (curve2, ovlp_arcs);
    CGAL_assertion (n_ovlps == 1);
    return (EQUAL);
  }

  // Decide wether curve1 is above, below or equal to curve2 immediately to
  // the right of the x co-ordinate of the given point.
  Comparison_result curve_compare_at_x_right (const X_curve& curve1, 
					      const X_curve& curve2,
					      const Point& p) const
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
    Point  ps1[2];

    n1 = curve1.get_points_at_x (p.x(), ps1);

    if (n1 == 0)    // p is not in the x-range of curve1.
      return (EQUAL);
    
    // Get the points on curve2 with the same x co-ordinate as p.
    int    n2;
    Point  ps2[2];

    n2 = curve2.get_points_at_x (p.x(), ps2);

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

    Comparison_result result = compare_y (ps1[0], ps2[0]);

    // In case the two curves do not intersect at the x co-ordinate of p,
    // just return the comparison result at p (since both curves are
    // continuous).    
    if (result != EQUAL)
      return (result);

    // Otherwise, the two curves do intersect at p_int = ps1[0] = ps2[0]:
    // make a decision based on their partial derivatives.
    const Point& p_int = ps1[0];

    // In order to simplify the process, make sure the source is always to the
    // right of p.
    X_curve c1 = curve1;
    X_curve c2 = curve2;
    
    if (compare_x(c1.source(), p) != LARGER)
      c1 = curve1.flip();
    if (compare_x(c2.source(), p) != LARGER)
      c2 = curve2.flip();

    // Calculate the partial derivatives at the intersection point.
    NT       dC1_dx, dC1_dy;
    NT       dC2_dx, dC2_dy;
    const NT zero = 0;
    
    c1.partial_derivatives (p_int, dC1_dx, dC1_dy);
    c2.partial_derivatives (p_int, dC2_dx, dC2_dy);

    // Try to compare the slopes of the two arcs.
    if (dC1_dy != zero && dC2_dy != zero)
    {
      Comparison_result  slope_result = _compare_value (dC2_dx/dC2_dy,
      							dC1_dx/dC1_dy);
      if (slope_result != EQUAL)
      {
	return (slope_result);
      }
    }

    // If we reached here, the two curves intersect at p_int, and their
    // tangents at this point have the same slope. We shall therefore find
    // a point with extreme y-coordinate for each curve and compare the
    // y-values of these two extreme points.
    Point  p_extr1, p_extr2;
    Point  hpts[2];
    int    n_hpts;

    n_hpts = c1.horizontal_tangency_points(hpts);
    if (n_hpts == 0)
    {
      // No horizontal tangency points - the source must have an extreme
      // y-value.
      p_extr1 = c1.source();
    }
    else if (n_hpts == 1)
    {
      // Use the vertical tangency points only if it lies between p_int
      // and the source. Otherwise, use the source.
      if (compare_x(c1.source(), hpts[0]) == LARGER &&
	  compare_x(p_int, hpts[0]) == SMALLER)
      {
	p_extr1 = hpts[0];
      }
      else
      {
	p_extr1 = c1.source();
      }
    }
    else
    {
      // An x-monotone curve cannot have more than 1 horizontal tangencies.
      CGAL_assertion (false);
    }

    n_hpts = c2.horizontal_tangency_points(hpts);
    if (n_hpts == 0)
    {
      // No horizontal tangency points - the source must have an extreme
      // y-value.
      p_extr2 = c2.source();
    }
    else if (n_hpts == 1)
    {
      // Use the vertical tangency points only if it lies between p_int
      // and the source. Otherwise, use the source.
      if (compare_x(c2.source(), hpts[0]) == LARGER &&
	  compare_x(p_int, hpts[0]) == SMALLER)
      {
	p_extr2 = hpts[0];
      }
      else
      {
	p_extr2 = c2.source();
      }
    }
    else
    {
      // An x-monotone curve cannot have more than 1 horizontal tangencies.
      CGAL_assertion (false);
    }

    // Use the x-coordinate of the point next to x.
    result = compare_x (p_extr1, p_extr2);

    if (result == LARGER)
    {
      n_hpts = c1.get_points_at_x (p_extr2.x(), hpts);
      CGAL_assertion(n_hpts == 1);
      p_extr1 = hpts[0];
    }
    else if (result == SMALLER)
    {
      n_hpts = c2.get_points_at_x (p_extr1.x(), hpts);
      CGAL_assertion(n_hpts == 1);
      p_extr2 = hpts[0];
    }

    // Compare the y-coordinates of the two extreme points.
    Comparison_result extr_result = compare_y (p_extr1, p_extr2);

    if (extr_result != EQUAL)
      return (extr_result);

    // As a last resort:
    NT    x_mid = (p_extr1.x() + p_int.x()) / NT(2);
    Point p_mid1, p_mid2;

    n_hpts = c1.get_points_at_x (x_mid, hpts);
    CGAL_assertion(n_hpts == 1);
    p_mid1 = hpts[0];

    n_hpts = c2.get_points_at_x (x_mid, hpts);
    CGAL_assertion(n_hpts == 1);
    p_mid2 = hpts[0];

    Comparison_result mid_result = compare_y (p_mid1, p_mid2);

    if (mid_result != EQUAL)
      return (mid_result);

    // If we reached here, the two curves must be overlapping:
    int     n_ovlps;
    X_curve ovlp_arcs[2];

    n_ovlps = curve1.overlaps (curve2, ovlp_arcs);
    CGAL_assertion (n_ovlps == 1);
    return (EQUAL);
  }

  // Check whether the given point is above, under or on the given curve.
  Curve_point_status curve_get_point_status (const X_curve& curve,
					     const Point& p) const
  {
    CGAL_precondition(is_x_monotone(curve));

    // A special treatment for vertical segments:
    if (curve.is_vertical_segment())
    {
      if (compare_x (curve.source(), p) != EQUAL)
	return (CURVE_NOT_IN_RANGE);

      if (compare_y (curve.source(), p) == SMALLER &&
	  compare_y (curve.target(), p) == SMALLER)
      {
	return (ABOVE_CURVE);
      }
      else if (compare_y (curve.source(), p) == LARGER &&
	       compare_y (curve.target(), p) == LARGER)
      {
	return (UNDER_CURVE);
      }
      else
      {
	return (ON_CURVE);
      }
    }

    // Get the points on the arc with the same x co-ordinate as p.
    int    n;
    Point  ps[2];

    n = curve.get_points_at_x (p.x(), ps);

    // Make sure there is at most one point.
    CGAL_assertion(n <= 1);

    if (n == 0)
    {
      // p is not in the x-range of the curve.
      return (CURVE_NOT_IN_RANGE);
    }
    else
    {
     // Compare p with the a point of the curve with the same x co-ordinate.
      int result = compare_y (p, ps[0]);

      if (result == SMALLER)
	return (UNDER_CURVE);
      else if (result == LARGER)
	return (ABOVE_CURVE);
      else if (result == EQUAL)
	return (ON_CURVE);
    }

    // We should never reach here:
    CGAL_assertion(false);
    return (ON_CURVE);
  }

  // Check whether the given curve in between c1 and c2, when going in the
  // clockwise direction from p from c1 to c2.
  // Notice that all three curves share the same end-point p.
  bool curve_is_between_cw (const X_curve& curve,
			    const X_curve& c1, const X_curve& c2,
			    const Point& p) const
  {
    CGAL_precondition(is_x_monotone(curve));
    CGAL_precondition(is_x_monotone(c1));
    CGAL_precondition(is_x_monotone(c2));

    // Make sure that p is the source of all curves (otherwise flip them).
    X_curve cv1 = c1;
    X_curve cv2 = c2;
    X_curve cvx = curve;

    if (cv1.source() != p)
      cv1 = c1.flip();
    if (cv2.source() != p)
      cv2 = c2.flip();
    if (cvx.source() != p)
      cvx = curve.flip();

    // Decide whether each arc is defined to the left or to the right of p.
    bool cv1_is_left = (compare_x(cv1.target(), p) == SMALLER);
    bool cv2_is_left = (compare_x(cv2.target(), p) == SMALLER);
    bool cvx_is_left = (compare_x(cvx.target(), p) == SMALLER);

    if (cv1.is_vertical_segment())
    {
      // cv1 is a vertical segment:
      bool cv1_is_up = (compare_y(cv1.target(), p) == LARGER);

      if (cv2.is_vertical_segment())
      {
	// Both cv1 and cv2 are vertical segments:
	bool cv2_is_up = (compare_y(cv2.target(), p) == LARGER);

	if (cv1_is_up && !cv2_is_up)
	  return (!cvx_is_left);
	else if (!cv1_is_up && cv2_is_up)
	  return (cvx_is_left);
	else
	  return (false);
      }
      
      if (cv1_is_up)
      {
	if (cv2_is_left)
	{
	  return ((!cvx_is_left) ||
		  (cvx_is_left &&
		   curve_compare_at_x_left (cv2, cvx, p) == LARGER));
	}
	else
	{
	  return (!cvx_is_left &&
		  curve_compare_at_x_right (cv2, cvx, p) == SMALLER);
	}
      }
      else
      {
	if (cv2_is_left)
	{
	  return (cvx_is_left &&
		  curve_compare_at_x_left (cv2, cvx, p) == LARGER);
	}
	else
	{
	  return ((cvx_is_left) ||
		  (!cvx_is_left &&
		   curve_compare_at_x_right (cv2, cvx, p) == SMALLER));
	}
      }
    }

    if (cv2.is_vertical_segment())
    {
      // Only cv2 is a vertical segment:
      bool cv2_is_up = (compare_y(cv2.target(), p) == LARGER);

      if (cv2_is_up)
      {
	if (cv1_is_left)
	{
	  return (cvx_is_left &&
		  curve_compare_at_x_left (cv1, cvx, p) == SMALLER);
	}
	else
	{
	  return ((cvx_is_left) ||
		  (!cvx_is_left &&
		   curve_compare_at_x_right (cv1, cvx, p) == LARGER));
	}
      }
      else
      {
	if (cv1_is_left)
	{
	  return ((!cvx_is_left) ||
		  (cvx_is_left &&
		   curve_compare_at_x_left (cv1, cvx, p) == SMALLER));
	}
	else
	{
	  return (cvx_is_left &&
		  curve_compare_at_x_right (cv1, cvx, p) == LARGER);
	}
      }
    }

    if (cvx.is_vertical_segment())
    {
      // Only cvx is a vertical segment:
      bool cvx_is_up = (compare_y(cvx.target(), p) == LARGER);

      if (cvx_is_up)
      {
	if (cv1_is_left && !cv2_is_left)
	  return (true);
	else if (cv1_is_left && cv2_is_left)
	  return (curve_compare_at_x_left(cv1, cv2, p) == LARGER);
	else if (!cv1_is_left && !cv2_is_left) 
	  return (curve_compare_at_x_right(cv1, cv2, p) == SMALLER);
      }
      else
      {
	if (!cv1_is_left && cv2_is_left)
	  return (true);
	else if (cv1_is_left && cv2_is_left)
	  return (curve_compare_at_x_left(cv1, cv2, p) == SMALLER);
	else if (!cv1_is_left && !cv2_is_left)
	  return (curve_compare_at_x_right(cv1, cv2, p) == LARGER);
      }

      return (false);
    }

    // None of the three curves is a vertical segment:
    // Check the following 4 cases:
    if (cv1_is_left && cv2_is_left)
    {
      // Case 1: Both c1 and c2 are defined to the left of p.
      if (curve_compare_at_x_left (cv1, cv2, p) == LARGER)
      { 
	// c1 is above c2:
        return (!(curve_compare_at_x_left (cv2, cvx, p) == SMALLER &&
                  curve_compare_at_x_left (cv1, cvx, p) == LARGER));
      }
      else
      { 
	// c2 is above c1:
        return (curve_compare_at_x_left (cv1, cvx, p) == SMALLER &&
                curve_compare_at_x_left (cv2, cvx, p) == LARGER);
      }
    }
    else if (!cv1_is_left && !cv2_is_left) 
    {
      // Case 2: Both c1 and c2 are defined to the right of p.
      if (curve_compare_at_x_right (cv1, cv2, p) == LARGER)
      {
	// c1 is above c2:
        return (curve_compare_at_x_right (cv2, cvx, p) == SMALLER &&
                curve_compare_at_x_right (cv1, cvx, p) == LARGER);
      }
      else
      { 
	// c2 is above c1:
        return (!(curve_compare_at_x_right (cv1, cvx, p) == SMALLER &&
                  curve_compare_at_x_right (cv2, cvx, p) == LARGER));
      }
    }
    else if (cv1_is_left && !cv2_is_left)
    {
      // Case 3: c1 is defined to the left and c2 is to the right of p.
      if (cvx_is_left)
        return (curve_compare_at_x_left(cv1, cvx, p) == SMALLER);
      else
        return (curve_compare_at_x_right(cv2, cvx, p) == SMALLER); 
    }
    else if (!cv1_is_left && cv2_is_left)
    {
      // Case 4: c1 is defined to the right and c2 is to the left of p. 
      if (cvx_is_left)
        return (curve_compare_at_x_left (cv2, cvx, p) == LARGER);
      else
        return (curve_compare_at_x_right(cv1, cvx, p) == LARGER); 
    }

    // We should never reach here.
    CGAL_assertion(false);
    return false;
  }

  // Check whether the two curves are identical.
  bool curve_is_same (const X_curve& curve1, const X_curve& curve2) const
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    // Check whether all arc features are the same.
    if (curve1.conic().orientation() == curve2.conic().orientation())
    {
      // Same orientation:
      return (curve1.conic() == curve2.conic() &&
	      curve1.source() == curve2.source() &&
	      curve1.target() == curve2.target());
    }
    else
    {
      // Check the flip case:
      return (curve1.conic() == curve2.conic() &&
	      curve1.source() == curve2.target() &&
	      curve1.target() == curve2.source());
    }
  }

  // Get the source and target vertex of the curve.
  Point curve_source(const X_curve& curve) const
  {
    return (curve.source());
  }

  Point curve_target(const X_curve& curve) const
  {
    return (curve.target());
  }

  // Return a point to the left or to the right of p.
  Point point_to_left (const Point& p) const
  {
    return (Point(p.x()-NT(1),p.y()));
  }

  Point point_to_right (const Point& p) const
  {
    return (Point(p.x()+NT(1),p.y()));
  }

  // Reflect a point in y.
  Point point_reflect_in_y (const Point& p) const
  {
    // Use hx(), hy(), hw() in order to support both Homogeneous and Cartesian.
    return (Point (-p.hx(), p.hy(), p.hw()));
  }
      
  // Reflect a curve in y.
  X_curve curve_reflect_in_y (const X_curve& curve) const
  {
    Conic  ref_conic (curve.conic().r(),
                      curve.conic().s(),
		      -curve.conic().t(),
		      -curve.conic().u(),
		      curve.conic().v(),
		      curve.conic().w());
    X_curve ref_arc (ref_conic,
		     point_reflect_in_y (curve.source()),
		     point_reflect_in_y (curve.target()));
    return (ref_arc);
  }

  // Reflect a point in x and y.
  Point point_reflect_in_x_and_y (const Point& p) const
  {
    // Use hx(), hy(), hw() in order to support both Homogeneous and Cartesian.
    return (Point (-p.hx(), -p.hy(), p.hw()));
  }
      
  // Reflect a curve in x and y.
  X_curve curve_reflect_in_x_and_y (const X_curve& curve) const
  {
    Conic  ref_conic (curve.conic().r(),
                      curve.conic().s(),
		      curve.conic().t(),
		      -curve.conic().u(),
		      -curve.conic().v(),
		      curve.conic().w());
    X_curve ref_arc (ref_conic,
		     point_reflect_in_x_and_y (curve.source()),
		     point_reflect_in_x_and_y (curve.target()));
    return (ref_arc);
  }

  
  ////////// Arrangement methods: //////////

  // Change the orientation of the curve (swap the source and the target).
  X_curve curve_flip (const X_curve& curve) const
  {
    CGAL_precondition(is_x_monotone(curve));

    // Flip the arc.
    return (curve.flip());
  }

  // Check whether the curve is x-monotone.
  bool is_x_monotone (const Curve& curve) const
  {
    return (curve.is_x_monotone());
  }

  // Cut the curve to several x-monotone sub-curves.
  void make_x_monotone (const Curve& curve, 
			std::list<X_curve>& x_curves) const
  {
    CGAL_precondition(!is_x_monotone(curve));

    // Clear the output list.
    x_curves.clear();

    // Find the points of vertical tangency and act accordingly.
    int    n;
    Point  ps[2];

    n = curve.vertical_tangency_points (ps);

    CGAL_assertion (n > 0);

    // Split the conic arc into x-monotone sub-curves. 
    if (curve.is_full_conic())
    {
      // Make sure we have two vertical tangency points.
      CGAL_assertion(n == 2);

      // In case the curve is a full conic, split it to two x-monotone curves,
      // one going from ps[0] to ps[1], and the other from ps[1] to ps[0].
      x_curves.push_back (X_curve (curve.conic(), ps[0], ps[1]));
      x_curves.push_back (X_curve (curve.conic(), ps[1], ps[0]));
    }
    else
    {
      X_curve    sub_curve1;
      X_curve    sub_curve2;
      X_curve    sub_curve3;

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
	X_curve    temp;

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
  void curve_split (const X_curve& curve, 
		    X_curve& sub_curve1, X_curve& sub_curve2, 
                    const Point& p) const 
  {
    CGAL_precondition(is_x_monotone(curve));

    // Make sure the point is on the curve and is not an end-point.
    CGAL_precondition(curve.contains_point(p));
    CGAL_precondition(curve_source(curve) != p);
    CGAL_precondition(curve_target(curve) != p);

    // Split the curve.
    _curve_split (curve,
		  sub_curve1, sub_curve2,
		  p);
    return;
  }

  // Check whether the intersection point between the two given curves is
  // lexicographically strictly to right of the given point.
  bool do_intersect_to_right(const X_curve& curve1, const X_curve& curve2,
                             const Point& p) const 
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    // Deal with overlapping curves:
    int     n_ovlps;
    X_curve ovlp_arcs[2];

    n_ovlps = curve1.overlaps (curve2, ovlp_arcs);
    CGAL_assertion (n_ovlps < 2);

    if (n_ovlps == 1)
    {
      // Check if at least one end-point of the overlapping curve is to the
      // right of p.
      return (compare_x (ovlp_arcs[0].source(), p) == LARGER ||
	      compare_x (ovlp_arcs[0].target(), p) == LARGER);
    }

    // In case there are no overlaps and the base conics are the same,
    // there cannot be any intersection points, unless the two x-monotone
    // curves share an end point.
    if (curve1.conic() == curve2.conic())
    {
      if ((curve1.source() == curve2.source() ||
	   curve1.source() == curve2.target()) &&
	  compare_x (curve1.source(), p) == LARGER)
      {
	return (true);
      }

      if ((curve1.target() == curve2.source() ||
	   curve1.target() == curve2.target()) &&
	  compare_x (curve1.target(), p) == LARGER)
      {
	return (true);
      }
      
      // No overlaps at all: the two curves do not intersect.
      return (false);
    }

    // Find the intersection points and decide accordingly.
    int   n;
    Point ps[4];
  
    n = curve1.intersections_with (curve2, ps);

    for (int i = 0; i < n; i++)
    {
      if (compare_x (ps[i], p) == LARGER)
	return (true);
    }

    return (false);
  }

  // Find the nearest intersection point between the two given curves to the
  // right of the given point.
  // In case of an overlap, p1 and p2 are the source and destination of the
  // overlapping curve. Otherwise p1=p2 is the calculated intersection point.
  bool nearest_intersection_to_right (const X_curve& curve1,
				      const X_curve& curve2,
				      const Point& p,
                                      Point& p1,
                                      Point& p2) const
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    // Deal with overlapping curves:
    int     n_ovlps;
    X_curve ovlp_arcs[2];

    n_ovlps = curve1.overlaps (curve2, ovlp_arcs);
    CGAL_assertion (n_ovlps < 2);

    if (n_ovlps == 1)
    {
      Point  ovlp_source = ovlp_arcs[0].source();
      Point  ovlp_target = ovlp_arcs[0].target();

      if (compare_x (ovlp_source, p) == LARGER &&
	  compare_x (ovlp_target, p) == LARGER)
      {
	// The entire overlapping arc is to the right of p:
	p1 = ovlp_source;
	p2 = ovlp_target;
	return (true);
      }
      else if (compare_x (ovlp_source, p) != LARGER &&
	       compare_x (ovlp_target, p) == LARGER)
      {
	// The source is to the left of p, and the traget is to its right.
	p1 = p;
	p2 = ovlp_target;
	return (true);
      }
      else if (compare_x (ovlp_source, p) == LARGER &&
	       compare_x (ovlp_target, p) != LARGER)
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
    if (curve1.conic() == curve2.conic())
    {
      const Point    *nearest_end_P = NULL;

      if ((curve1.source() == curve2.source() ||
	   curve1.source() == curve2.target()) &&
	  compare_x (curve1.source(), p) == LARGER)
      {
	nearest_end_P = &(curve1.source());
      }

      if ((curve1.target() == curve2.source() ||
	   curve1.target() == curve2.target()) &&
	  compare_x (curve1.target(), p) == LARGER)
      {
	if (nearest_end_P == NULL ||
	    compare_x (*nearest_end_P, curve1.target()) == LARGER)
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
    Point        ps[4];
    const Point *nearest_inter_P = NULL;
  
    n = curve1.intersections_with (curve2, ps);

    for (int i = 0; i < n; i++)
    {
      if (compare_x (ps[i], p) == LARGER)
      {
	if (nearest_inter_P == NULL)
	{
	  // The first point to the right so far:
	  nearest_inter_P = &(ps[i]);
	}
	else
	{
	  // Compare with the nearest point so far.
	  if (compare_x (ps[i], *nearest_inter_P) == SMALLER)
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

  // Check whether two curves overlap.
  bool curves_overlap (const X_curve& curve1, const X_curve& curve2) const
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    X_curve ovlp_arcs[2];

    return (curve1.overlaps (curve2, ovlp_arcs) > 0);
  }

 private:

  ////////// Private auxiliary methods: //////////

  // Compare two values.
  Comparison_result _compare_value (const NT& a, const NT& b) const
  {
    return CGAL::compare(a,b);
  }

  // Split the given curve into two sub-curves at the given point.
  // Since this is a private function, there are no preconditions.
  void _curve_split (const X_curve& curve, 
		     X_curve& sub_curve1, X_curve& sub_curve2, 
		     const Point& p) const
  {
    // Split the curve to source->p and p->target.
    sub_curve1 = X_curve (curve.conic(), curve.source(), p);
    sub_curve2 = X_curve (curve.conic(), p, curve.target());

    return;
  }

};

CGAL_END_NAMESPACE

#endif






