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
// release       : $CGAL_Revision: CGAL-2.4-I-5 $
// release_date  : $CGAL_Date: 2001/08/31 $
//
// file          : include/CGAL/Arr_segment_circle_traits.h
// package       : Arrangement (2.19)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Ron Wein <wein@post.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_ARR_CONIC_TRAITS_2_H
#define CGAL_ARR_CONIC_TRAITS_2_H

#include <CGAL/basic.h>
#include <list>

#include <CGAL/Segment_circle_2.h>
#include <CGAL/tags.h>

CGAL_BEGIN_NAMESPACE

// ----------------------------------------------------------------------------
// Arrangement traits for conic arcs.
//

template <class _NT>
class Arr_segment_circle_traits 
{
 public:
  typedef Tag_false                 Has_left_category;
  
  typedef _NT                       NT;

  // The difference between Curve and X_curve is semantical only,
  // NOT syntactical.
  typedef Segment_circle_2<NT>      Curve_2;
  typedef Curve_2                   X_curve_2;
 
  // Using typename to please compiler (e.g., CC with IRIX64 on mips)
  typedef typename Curve_2::R       R;
  typedef typename Curve_2::Point   Point_2;
  typedef typename Curve_2::Segment Segment_2;
  typedef typename Curve_2::Circle  Circle_2;
  typedef typename Curve_2::Conic   Conic_2;

  // Obsolete, for backward compatibility
  typedef Point_2                   Point;
  typedef X_curve_2                 X_curve;
  typedef Curve_2                   Curve;
  typedef Segment_2                 Segment;
  typedef Circle_2                  Circle;
  typedef Conic_2                   Conic;

  // Constructor.
  Arr_segment_circle_traits()
    {}

  ////////// Planar Map methods: //////////

  // Compare the x co-ordinates of two given points.
  Comparison_result compare_x(const Point_2 & p0, const Point_2 & p1) const
  {
    return _compare_value(p0.x(),p1.x());
  }

  // Compare the two points lexicographically (by x, then by y).
  Comparison_result compare_xy(const Point_2 & p0, const Point_2 & p1) const
  {
    Comparison_result x_res = _compare_value(p0.x(),p1.x());

    if (x_res != EQUAL)
      return (x_res);

    return _compare_value(p0.y(),p1.y());
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
      return (compare_x (curve.source(), p) == EQUAL);
    }
    else
    {
      // Find the number of points on the arc with the same x co-ordinate as p.
      // p is is the x-range of the arc only if there is at least one point.
      Point_2  ps[2];
      int    n = curve.get_points_at_x (p.x(), ps);

      return (n > 0);
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
    CGAL_precondition(curve_is_in_x_range(curve1, p));
    CGAL_precondition(curve_is_in_x_range(curve2, p));

    // Get the points on curve1 with the same x co-ordinate as p.
    int      n1;
    Point_2  ps1[2];

    if (curve1.is_vertical_segment())
    {
      // Vertical segment.
      if (compare_x (curve1.source(), p) != EQUAL)
	return (EQUAL);

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
    else
    {
      n1 = curve1.get_points_at_x (p.x(), ps1);

      CGAL_assertion(n1 == 1);
    }

    // Get the points on curve2 with the same x co-ordinate as p.
    int      n2;
    Point_2  ps2[2];

    if (curve2.is_vertical_segment())
    {
      // Vertical segment.
      if (compare_x (curve2.source(), p) != EQUAL)
	return (EQUAL);
      
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
    else
    {
      n2 = curve2.get_points_at_x (p.x(), ps2);

      CGAL_assertion(n2 == 1);
    }

    // Deal with vertical segments:
    if (n1 == 2)
    {
      // Check if the vertical segment contains ps2[0].
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
      // Check if the vertical segment contains ps1[0].
      if (_compare_y (ps2[0], ps1[0]) != LARGER && 
	  _compare_y (ps2[1], ps1[0]) != SMALLER)
      {
	return (EQUAL);
      }
    }

    // Compare the y co-ordinates of the two points.
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

    // The two curve must not be vertical segments.
    CGAL_precondition(! curve1.is_vertical_segment());
    CGAL_precondition(! curve2.is_vertical_segment());

    // Check that both curves are defined to the left of p.
    CGAL_precondition((compare_x(curve1.source(), p) == SMALLER) ||
		      (compare_x(curve1.target(), p) == SMALLER));

    CGAL_precondition((compare_x(curve2.source(), p) == SMALLER) ||
		      (compare_x(curve2.target(), p) == SMALLER));

    // Get the points on curve1 with the same x co-ordinate as p.
    int      n1;
    Point_2  ps1[2];

    if (curve1.contains_point(p))
    {
      ps1[0] = p;
      n1 = 1;
    }
    else
    {
      n1 = curve1.get_points_at_x (p.x(), ps1);
    }

    // Make sure that there is exactly one point.
    CGAL_assertion(n1 == 1);
    
    // Get the points on curve2 with the same x co-ordinate as p.
    int      n2;
    Point_2  ps2[2];

    if (curve2.contains_point(p))
    {
      ps2[0] = p;
      n2 = 1;
    }
    else
    {
      n2 = curve2.get_points_at_x (p.x(), ps2);
    }

    // Make sure that there is exactly one point.
    CGAL_assertion(n2 == 1);

    // Make sure the two curves intersect at x(p).
    CGAL_precondition(_compare_y (ps1[0], ps2[0]) == EQUAL);

    // If the curves are the same, they are equal to the left of p:
    if (curve_is_same(curve1,curve2))
      return (EQUAL);

    // Otherwise, the two curves do intersect at p_int = ps1[0] = ps2[0]:
    // make a decision based on their partial derivatives.
    const Point_2& p_int = ps1[0];
    
    // In order to simplify the process, make sure the source is always to the
    // left of p.
    X_curve_2 c1 = curve1;
    X_curve_2 c2 = curve2;
    
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
    Point_2  p_extr1, p_extr2;
    Point_2  hpts[2];
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
    Comparison_result result = compare_x (p_extr1, p_extr2);

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
    Comparison_result extr_result = _compare_y (p_extr1, p_extr2);

    if (extr_result != EQUAL)
      return (extr_result);

    // As a last resort:
    NT    x_mid = (p_extr1.x() + p_int.x()) / NT(2);
    Point_2 p_mid1, p_mid2;

    n_hpts = c1.get_points_at_x (x_mid, hpts);
    CGAL_assertion(n_hpts == 1);
    p_mid1 = hpts[0];

    n_hpts = c2.get_points_at_x (x_mid, hpts);
    CGAL_assertion(n_hpts == 1);
    p_mid2 = hpts[0];

    Comparison_result mid_result = _compare_y (p_mid1, p_mid2);

    if (mid_result != EQUAL)
      return (mid_result);

    // If we reached here, the two curves must be overlapping:
    int     n_ovlps;
    X_curve_2 ovlp_arcs[2];

    n_ovlps = curve1.overlaps (curve2, ovlp_arcs);
    CGAL_assertion (n_ovlps == 1);
    return (EQUAL);
  }

  // Decide wether curve1 is above, below or equal to curve2 immediately to
  // the right of the x co-ordinate of the given point.
  Comparison_result curve_compare_at_x_right (const X_curve_2& curve1, 
					      const X_curve_2& curve2,
					      const Point_2& p) const
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    // The two curve must not be vertical segments.
    CGAL_precondition(! curve1.is_vertical_segment());
    CGAL_precondition(! curve2.is_vertical_segment());

    // Check that both curves are defined to the right of p.
    CGAL_precondition((compare_x(curve1.source(), p) == LARGER) ||
		      (compare_x(curve1.target(), p) == LARGER));

    CGAL_precondition((compare_x(curve2.source(), p) == LARGER) ||
		      (compare_x(curve2.target(), p) == LARGER));

    // Get the points on curve1 with the same x co-ordinate as p.
    int      n1;
    Point_2  ps1[2];

    if (curve1.contains_point(p))
    {
      ps1[0] = p;
      n1 = 1;
    }
    else
    {
      n1 = curve1.get_points_at_x (p.x(), ps1);
    }

    // Make sure we have a single point.
    CGAL_assertion(n1 == 1);
    
    // Get the points on curve2 with the same x co-ordinate as p.
    int      n2;
    Point_2  ps2[2];

    if (curve2.contains_point(p))
    {
      ps2[0] = p;
      n2 = 1;
    }
    else
    {
      n2 = curve2.get_points_at_x (p.x(), ps2);
    }

    // Make sure we have a single point.
    CGAL_assertion(n2 == 1);

    // Make sure the two curves intersect at x(p).
    CGAL_precondition(_compare_y (ps1[0], ps2[0]) == EQUAL);

    // If the two curves are the same, they are equal to the right of p:
    if (curve_is_same(curve1,curve2))
      return (EQUAL);

    // Otherwise, the two curves do intersect at p_int = ps1[0] = ps2[0]:
    // make a decision based on their partial derivatives.
    const Point_2& p_int = ps1[0];

    // In order to simplify the process, make sure the source is always to the
    // right of p.
    X_curve_2 c1 = curve1;
    X_curve_2 c2 = curve2;
    
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
    Point_2  p_extr1, p_extr2;
    Point_2  hpts[2];
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
    Comparison_result result = compare_x (p_extr1, p_extr2);

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
    Comparison_result extr_result = _compare_y (p_extr1, p_extr2);

    if (extr_result != EQUAL)
      return (extr_result);

    // As a last resort:
    NT    x_mid = (p_extr1.x() + p_int.x()) / NT(2);
    Point_2 p_mid1, p_mid2;

    n_hpts = c1.get_points_at_x (x_mid, hpts);
    CGAL_assertion(n_hpts == 1);
    p_mid1 = hpts[0];

    n_hpts = c2.get_points_at_x (x_mid, hpts);
    CGAL_assertion(n_hpts == 1);
    p_mid2 = hpts[0];

    Comparison_result mid_result = _compare_y (p_mid1, p_mid2);

    if (mid_result != EQUAL)
      return (mid_result);

    // If we reached here, the two curves must be overlapping:
    int     n_ovlps;
    X_curve_2 ovlp_arcs[2];

    n_ovlps = curve1.overlaps (curve2, ovlp_arcs);
    CGAL_assertion (n_ovlps == 1);
    return (EQUAL);
  }

  // Check whether the given point is above, under or on the given curve.
  Comparison_result curve_get_point_status (const X_curve_2& curve,
					    const Point_2& p) const
  {
    CGAL_precondition(is_x_monotone(curve));
    CGAL_precondition(curve_is_in_x_range(curve, p));

    // A special treatment for vertical segments:
    if (curve.is_vertical_segment())
    {
      Comparison_result res1 = _compare_y (curve.source(), p);
      Comparison_result res2 = _compare_y (curve.target(), p);

      if (res1 == res2)
	return (res1);
      else
	return (EQUAL);
    }

    // Get the points on the arc with the same x co-ordinate as p.
    int    n;
    Point_2  ps[2];

    n = curve.get_points_at_x (p.x(), ps);

    // Make sure there is at most one point.
    CGAL_assertion(n == 1);

    // Compare p with the a point of the curve with the same x co-ordinate.
    return (_compare_y (ps[0], p));
  }

  // Check whether the two curves are identical.
  bool curve_is_same (const Point_2 & p, const Point_2 & q) const
  {
    return (p == q);
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
  Point_2 curve_source(const X_curve_2& curve) const
  {
    return (curve.source());
  }

  Point_2 curve_target(const X_curve_2& curve) const
  {
    return (curve.target());
  }

  // Reflect a point in y.
  Point_2 point_reflect_in_y (const Point_2& p) const
  {
    // Use hx(), hy(), hw() in order to support both Homogeneous and Cartesian.
    return (Point_2 (-p.hx(), p.hy(), p.hw()));
  }
      
  // Reflect a curve in y.
  X_curve_2 curve_reflect_in_y (const X_curve_2& curve) const
  {
    Conic  ref_conic (curve.conic().r(),
                      curve.conic().s(),
		      -curve.conic().t(),
		      -curve.conic().u(),
		      curve.conic().v(),
		      curve.conic().w());
    X_curve_2 ref_arc (ref_conic,
		     point_reflect_in_y (curve.source()),
		     point_reflect_in_y (curve.target()));
    return (ref_arc);
  }

  // Reflect a point in x and y.
  Point_2 point_reflect_in_x_and_y (const Point_2& p) const
  {
    // Use hx(), hy(), hw() in order to support both Homogeneous and Cartesian.
    return (Point_2 (-p.hx(), -p.hy(), p.hw()));
  }
      
  // Reflect a curve in x and y.
  X_curve_2 curve_reflect_in_x_and_y (const X_curve_2& curve) const
  {
    Conic  ref_conic (curve.conic().r(),
                      curve.conic().s(),
		      curve.conic().t(),
		      -curve.conic().u(),
		      -curve.conic().v(),
		      curve.conic().w());
    X_curve_2 ref_arc (ref_conic,
		     point_reflect_in_x_and_y (curve.source()),
		     point_reflect_in_x_and_y (curve.target()));
    return (ref_arc);
  }

  
  ////////// Arrangement methods: //////////

  // Change the orientation of the curve (swap the source and the target).
  X_curve_2 curve_flip (const X_curve_2& curve) const
  {
    CGAL_precondition(is_x_monotone(curve));

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
    if (is_x_monotone(curve))
    {
      x_curves.clear();
      x_curves.push_back(X_curve(curve));
      return;
    } 

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
      x_curves.push_back (X_curve_2 (curve.conic(), ps[0], ps[1]));
      x_curves.push_back (X_curve_2 (curve.conic(), ps[1], ps[0]));
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
    CGAL_precondition(curve_source(curve) != p);
    CGAL_precondition(curve_target(curve) != p);

    // Split the curve.
    _curve_split (curve,
		  sub_curve1, sub_curve2,
		  p);
    return;
  }

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

      if (compare_lexicographically_xy (ovlp_source, p) == LARGER &&
	  compare_lexicographically_xy (ovlp_target, p) == LARGER)
      {
	// The entire overlapping arc is to the right of p:
	p1 = ovlp_source;
	p2 = ovlp_target;
	return (true);
      }
      else if (compare_lexicographically_xy (ovlp_source, p) != LARGER &&
	       compare_lexicographically_xy (ovlp_target, p) == LARGER)
      {
	// The source is to the left of p, and the traget is to its right.
	p1 = p;
	p2 = ovlp_target;
	return (true);
      }
      else if (compare_lexicographically_xy (ovlp_source, p) == LARGER &&
	       compare_lexicographically_xy (ovlp_target, p) != LARGER)
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
      const Point_2    *nearest_end_P = NULL;

      if ((curve1.source() == curve2.source() ||
	   curve1.source() == curve2.target()) &&
	  compare_lexicographically_xy (curve1.source(), p) == LARGER)
      {
	nearest_end_P = &(curve1.source());
      }

      if ((curve1.target() == curve2.source() ||
	   curve1.target() == curve2.target()) &&
	  compare_lexicographically_xy (curve1.target(), p) == LARGER)
      {
	if (nearest_end_P == NULL ||
	    compare_lexicographically_xy (*nearest_end_P,
					  curve1.target()) == LARGER)
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
  
    n = curve1.intersections_with (curve2, ps);

    for (int i = 0; i < n; i++)
    {
      if (compare_lexicographically_xy (ps[i], p) == LARGER)
      {
	if (nearest_inter_P == NULL)
	{
	  // The first point to the right so far:
	  nearest_inter_P = &(ps[i]);
	}
	else
	{
	  // Compare with the nearest point so far.
	  if (compare_lexicographically_xy (ps[i], 
					    *nearest_inter_P) == SMALLER)
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
  bool curves_overlap (const X_curve_2& curve1, const X_curve_2& curve2) const
  {
    CGAL_precondition(is_x_monotone(curve1));
    CGAL_precondition(is_x_monotone(curve2));

    X_curve_2 ovlp_arcs[2];

    return (curve1.overlaps (curve2, ovlp_arcs) > 0);
  }

  // Check if the two points are the same.
  bool point_is_same(const Point_2 & p1, const Point_2 & p2) const
  { 
    return (compare_xy(p1, p2) == EQUAL);
  }
  
private:

  ////////// Private auxiliary methods: //////////

  // Compare two values.
  Comparison_result _compare_value (const NT& a, const NT& b) const
  {
    return CGAL::compare(a,b);
  }

  // Compare the y coordinates of the two points.
  Comparison_result _compare_y(const Point_2 & p0, const Point_2 & p1) const
  {
    return _compare_value(p0.y(),p1.y());
  }

  // Split the given curve into two sub-curves at the given point.
  // Since this is a private function, there are no preconditions.
  void _curve_split (const X_curve_2& curve, 
		     X_curve_2& sub_curve1, X_curve_2& sub_curve2, 
		     const Point_2 & p) const
  {
    // Split the curve to source->p and p->target.
    sub_curve1 = X_curve_2 (curve.conic(), curve.source(), p);
    sub_curve2 = X_curve_2 (curve.conic(), p, curve.target());

    return;
  }

};

CGAL_END_NAMESPACE

#endif
