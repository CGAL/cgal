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
// file          : include/CGAL/Segment_circle_2.h
// package       : Arrangement 
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Ron Wein <wein@post.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_SEGMENT_CIRCLE_2_H
#define CGAL_SEGMENT_CIRCLE_2_H

// Segment_circle_2.h
//
// A modified version that specializes in segments and circular arcs.
// This class does NOT support general conic arcs

#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H

#ifndef CGAL_CONIC_2_H
#include <CGAL/Conic_2.h>
#endif // CGAL_CONIC_2_H

#include <fstream>
#include <list>


CGAL_BEGIN_NAMESPACE

// ----------------------------------------------------------------------------
// Solve a quadratic equation.
// The function returns the number of distinct solutions.
// The roots area must be at least of size 2.
// 
template <class NT>
static int _solve_quadratic_eq (const NT& a, const NT& b, const NT& c,
				NT* roots)
{
  static const NT _zero = NT(0);
  static const NT _two  = NT(2);
  static const NT _four = NT(4);

  if (a == _zero)
  {
    // We have a linear equation.
    if (b != _zero)
    {
      roots[0] = -c / b;
      return (1);
    }
    else
      return (0);
  }

  const NT        disc = b*b - _four*a*c;

  if (disc < _zero)
  {
    // Negative discriminant - no real solution.
    return (0);
  }
  else if (disc == _zero)
  {
    // Zero discriminant - one real solution (with multiplicty 2).
    roots[0] = roots[1] = -b / (_two * a);
    return (1);
  }
  else
  {
    // Positive discriminant - two distinct real solutions.
    const NT      sqrt_disc = CGAL::sqrt(disc);

    roots[0] = (sqrt_disc - b) / (_two * a);
    roots[1] = -(sqrt_disc + b) / (_two * a);
    return (2);
  }
}

CGAL_END_NAMESPACE

// ----------------------------------------------------------------------------
// Representation of a conic arc which is either a segment (a curve of
// degree 1), or a circular arc (of degree 2).
//

CGAL_BEGIN_NAMESPACE

template <class NT>
class Segment_circle_2
{

 public:

  typedef Cartesian<NT>        R;
  typedef Point_2<R>           Point;
  typedef Conic_2<R>           Conic;
  typedef Circle_2<R>          Circle;
  typedef Segment_2<R>         Segment;

 private:

  Conic    _conic;              // The conic that contains the arc.
  int      _deg;                // The degree of the conic (either 1 or 2).
  Point    _source;             // The source of the arc. 
  Point    _target;             // The target of the arc.
  bool     _is_full;            // Indicated whether the arc is a full conic.
      
 public:

  // Default constructor.
  Segment_circle_2 () :
    _deg(0),
    _is_full(false)
  {}

  // Copy constructor.
  Segment_circle_2 (const Segment_circle_2<NT>& arc) :
    _conic(arc._conic),
    _deg(arc._deg),
    _source(arc._source),
    _target(arc._target),
    _is_full(arc._is_full)
  {}

  // Constuct a segment arc from a segment.
  Segment_circle_2 (const Segment& segment) :
    _deg(1),
    _source(segment.source()),
    _target(segment.target()),
    _is_full(false)
  {
    static const NT _zero = 0;
    static const NT _one = 1;

    // The supporting conic is r=s=t=0, and u*x + v*y + w = 0 should hold
    // for both the source (x1,y1) and the target (x2, y2).
    if (_source.x() == _target.x())
    {
      // The supporting conic is a vertical line, of the form x = CONST.
      _conic.set (_zero, _zero, _zero,    // r = s = t = 0
                  _one,                   // u = 1
		  _zero,                  // v = 0
		  -_source.x());          // w = -CONST
    }
    else
    {
      // The supporting line is y = A*x + B, where:
      //
      //       y2 - y1              x2*y1 - x1*y2 
      //  A = ---------        B = ---------------
      //       x2 - x1                 x2 - x1
      //
      const NT    A = (_target.y() - _source.y()) /
                      (_target.x() - _source.x());
      const NT    B = (_target.x()*_source.y() - _source.x()*_target.y()) /
                      (_target.x() - _source.x());

      // Now we can set:
      _conic.set (_zero, _zero, _zero,    // r = s = t = 0
		  A,                      // u = A
		  -_one,                  // v = -1
		  B);                     // w = B
    }        
  }

  // Construct a circular arc from a circle and two points on that circle
  // (the orientation of the circle is assumed to be positive: i.e. we move
  // in a clockwise direction from the source to the target).
  // The source and the target must be on the conic boundary and must
  // not be the same.
  Segment_circle_2 (const Circle& circle,
		    const Point& source, const Point& target) :
    _deg(2),
    _source(source),
    _target(target),
    _is_full(false)
  {
    // Make sure the circle contains the two endpoints on its boundary.
    CGAL_precondition(circle.has_on_boundary(_source));
    CGAL_precondition(circle.has_on_boundary(_target));

    // Make sure that the source and the taget are not the same.
    CGAL_precondition(_source != _target);      

    // Produce the correponding conic: if the circle centre is (x0,y0)
    // and it radius is r, that its equation is:
    //   x^2 + y^2 - 2*x0*x - 2*y0*y + (x0^2 + y0^2 - r^2) = 0
    static const NT _zero = 0;
    static const NT _one = 1;
    static const NT _minus_two = -2;
    const NT    x0 = circle.center().x();
    const NT    y0 = circle.center().y();
    const NT    r_squared = circle.squared_radius();

    _conic.set (_one, _one,                  // r = s = 1
		_zero,                       // t = 0
		_minus_two*x0,
		_minus_two*y0,
		x0*x0 + y0*y0 - r_squared);
  }

  // Construct an arc which is basically a full circle.
  Segment_circle_2 (const Circle& circle) :
    _deg(2),
    _is_full(true)
  {
    // Produce the correponding conic: if the circle centre is (x0,y0)
    // and it radius is r, that its equation is:
    //   x^2 + y^2 - 2*x0*x - 2*y0*y + (x0^2 + y0^2 - r^2) = 0
    static const NT _zero = 0;
    static const NT _one = 1;
    static const NT _minus_two = -2;
    const NT    x0 = circle.center().x();
    const NT    y0 = circle.center().y();
    const NT    r_squared = circle.squared_radius();

    _conic.set (_one, _one,                    // r = s = 1
		_zero,                         // t = 0
		_minus_two*x0,
		_minus_two*y0,
		x0*x0 + y0*y0 - r_squared);

    // Set a fictitious source and destination.
    _source = Point(x0 + CGAL::sqrt(r_squared), y0);
    _target = _source;
  }

  // Construct a conic arc.
  // The conic must either be a circle or a line.
  // The source and the target must be on the conic boundary and must
  // not be the same.
  Segment_circle_2 (const Conic& conic,
	            const Point& source, const Point& target) :
    _conic(conic),
    _source(source),
    _target(target),
    _is_full(false)
  {
    // Make sure the conic contains the two endpoints on its boundary.
    CGAL_precondition(_conic.has_on_boundary(_source));
    CGAL_precondition(_conic.has_on_boundary(_target));

    // Make sure that the source and the taget are not the same.
    CGAL_precondition(_source != _target);      

    // Check whether we have a segment.
    static const NT _zero = 0;

    if (conic.r() == _zero && conic.s() == _zero && conic.t() == _zero &&
	(conic.u() != _zero || conic.v() != _zero))
    {
      _deg = 1;
      return;
    }

    // If the conic is of degree 2, it must be a circle.
    CGAL_precondition(_conic.is_ellipse());      
    CGAL_precondition(_conic.r() != NT(0));      
    CGAL_precondition(_conic.r() == _conic.s());      
    CGAL_precondition(_conic.t() == NT(0));
    _deg = 2;
  }

  // Destructor.
  virtual ~Segment_circle_2 ()
  {}

  // Get the arc's base conic.
  const Conic& conic () const
  {
    return (_conic);
  }

  // Get the arc's source.
  const Point& source () const
  {
    return (_source);
  }

  // Get the arc's target.
  const Point& target () const
  {
    return (_target);
  }

  // Check whether the curve is a segment.
  bool is_segment() const
  {
    return (_deg == 1);
  }

  // Check whether the curve is a circle.
  bool is_circle() const
  {
    return (_deg == 2);
  }

  // Get a segment if the arc is indeed one.
  Segment segment() const
  {
    CGAL_precondition(is_segment());

    return (Segment (_source, _target));
  }

  // Get a segment if the arc is indeed one.
  Circle circle() const
  {
    CGAL_precondition(is_circle());

    // Create the appropriate circle.
    static const NT _zero = 0;
    static const NT _two = 2;
    NT              x0, y0, r2;

    if (_conic.r() > _zero)
    {
      // Positive orientation. The conic has the form:
      //  x^2 + y^2 - (2*x0)*x - (2*y0)*y + (x0^2 + y0^2 - r^2) = 0 
      x0 = -(_conic.u() / _two);
      y0 = -(_conic.v() / _two);
      r2 = x0*x0 + y0*y0 - _conic.w();
    }
    else
    {
      // Negative orientation:
      //  - x^2 - y^2 + (2*x0)*x + (2*y0)*y + (r^2 - x0^2 - y0^2) = 0 
      x0 = _conic.u() / _two;
      y0 = _conic.v() / _two;
      r2 = x0*x0 + y0*y0 + _conic.w();
    }

    return (Circle (Point(x0, y0), r2));
  }

  // Check whether the arc is a full conic (i.e. a full circle).
  bool is_full_conic () const
  {
    return (_is_full);
  }

  // Check whether the curve is a vertical segment.
  bool is_vertical_segment () const
  {
    // A vertical segment is contained in the degenerate conic: u*x + w = 0.
    static const NT _zero = 0;

    return (_deg == 1 && _conic.v() == _zero);
  }

  // Check whether the curve is a horizontal segment.
  bool is_horizontal_segment () const
  {
    // A vertical segment is contained in the degenerate conic: v*y + w = 0.
    static const NT _zero = 0;

    return (_deg == 1 && _conic.u() == _zero);
  }

  // Check whether the given point is on the arc.
  bool contains_point (const Point& p) const
  {
    // Check whether the conic contains the point (x,y).
    if (! _conic.has_on_boundary(p))
      return (false);
      
    // If the point is on the conic, make sure it is between the two endpoints.
    else
      return (is_full_conic() || _is_between_endpoints(p));
  }

  // Calculate the vertical tangency points of the arc (ps should be allocated
  // to the size of 2).
  // The function returns the number of vertical tangency points.
  int vertical_tangency_points (Point* vpts) const
  {
    // No vertical tangency points for segments:
    if (_deg < 2)
      return (0);

    // Calculate the vertical tangency points of the conic.
    Point  ps[2];
    int    n;

    n = _conic_vertical_tangency_points (ps);

    // Return only the points that are contained in the arc interior.
    int    m = 0;
    
    for (int i = 0; i < n; i++)
    {
      if (is_full_conic() || _is_strictly_between_endpoints(ps[i]))
      {
	vpts[m] = ps[i];
	m++;
      }
    }

    // Return the number of vertical tangency points found.
    return (m);
  }

  // Calculate the horizontal tangency points of the arc (ps should be
  // allocated to the size of 2).
  // The function return the number of vertical tangency points.
  int horizontal_tangency_points (Point* hpts) const
  {
    // No vertical tangency points for segments:
    if (_deg < 2)
      return (0);

    // Calculate the vertical tangency points of the conic.
    Point  ps[2];
    int    n;

    n = _conic_horizontal_tangency_points (ps);

    // Return only the points that are contained in the arc interior.
    int    m = 0;
    
    for (int i = 0; i < n; i++)
    {
      if (is_full_conic() || _is_strictly_between_endpoints(ps[i]))
      {
	hpts[m] = ps[i];
	m++;
      }
    }

    // Return the number of vertical tangency points found.
    return (m);
  }

  // Check whether the arc is x-monotone.
  bool is_x_monotone() const 
  {
    // Any segment (including vertical segments) is considered x-monotone:
    if (_deg < 2)
      return (true);

    // Check the number of vertical tangency points.
    Point   ps[2];
    
    return (vertical_tangency_points (ps) == 0);
  }

  // Find all points on the arc with a given x-coordinate: ps should be
  // allocated to the size of 2.
  // The function return the number of points found.
  int get_points_at_x (const NT& x,
                       Point *ps) const
  {
    // Make sure the conic is not a vertical segment.
    CGAL_precondition(!is_vertical_segment());

    // Get the y co-ordinates of the points on the conic.
    NT    ys[2];
    int   n;

    n = _conic_get_y_coordinates (x, ys);
    
    // Find all the points that are contained in the arc.
    int   m = 0;
    
    for (int i = 0; i < n; i++)
    {
      ps[m] = Point (x, ys[i]);

      if (is_full_conic() || _is_between_endpoints(ps[m]))
	m++;
    }

    // Return the number of points on the arc.
    return (m);
  }
  
  // Find all points on the arc with a given y-coordinate: ps should be
  // allocated to the size of 2.
  // The function return the number of points found.
  int get_points_at_y (const NT& y,
                       Point *ps) const
  {
    // Make sure the conic is not a horizontal segment.
    CGAL_precondition(!is_horizontal_segment());

    // Get the x co-ordinates of the points on the conic.
    NT    xs[2];
    int   n;

    n = _conic_get_x_coordinates (y, xs);
    
    // Find all the points that are contained in the arc.
    int   m = 0;
    
    for (int i = 0; i < n; i++)
    {
      ps[m] = Point (xs[i], y);

      if (is_full_conic() || _is_between_endpoints(ps[m]))
	m++;
    }

    // Return the number of points on the arc.
    return (m);
  }
  
  // Return a flipped conic arc: exchange the arc's source and destination.
  Segment_circle_2 flip () const
  {
    Conic           opp_conic (-_conic.r(), -_conic.s(), -_conic.t(),
			       -_conic.u(), -_conic.v(), -_conic.w());

    return (Segment_circle_2<NT> (opp_conic, _target, _source));
  }

  // Get the partial derivatives of the arc at a given point.
  void partial_derivatives (const Point& p, 
		            NT& dC_dx, NT& dC_dy) const
  {
    // Make sure p is contained in the arc.
    CGAL_precondition(contains_point(p));

    // Calulate the partial derivatives of the conic C at p=(x,y), which are:
    //
    //    dC                            dC
    //   ---- = 2rx + ty + u           ---- = 2sy + tx + v
    //    dx                            dy
    //
    static const NT _two = 2;
    
    dC_dx = _two*_conic.r()*p.x() + _conic.t()*p.y() + _conic.u();
    dC_dy = _two*_conic.s()*p.y() + _conic.t()*p.x() + _conic.v();

    return;
  }

  // Calculate the intersection points between the arc and the given arc.
  // ps must be allocated at the size of 4.
  // The function returns the number of actual intersection point.
  int intersections_with (const Segment_circle_2<NT>& arc,
			  Point* ps) const
  {
    // For simplicity, assume that (*this) has the higher degree.
    if (_deg == 1 && arc._deg == 2)
    {
      return (arc.intersections_with (*this,
				      ps));
    }

    // Check the case when one of the two arcs is a vertical segment.
    const Segment_circle_2<NT>* vertical_P = NULL;
    const Segment_circle_2<NT>* other_P = NULL;
    
    if (is_vertical_segment())
    {
      vertical_P = this;
      other_P = &arc;
    }
    if (arc.is_vertical_segment())
    {
      // Both arcs are vertical segments: there should be no intersections.
      if (vertical_P != NULL)
	return (0);
      
      vertical_P = &arc;
      other_P = this;
    }

    if (vertical_P != NULL)
    {
      // Find all points on the other arc that have the segment's x value.
      Point  xps[2];
      int    n_ys;

      n_ys = other_P->get_points_at_x (vertical_P->source().x(), xps);

      // Make sure those points are on the vertical segment.
      int      n = 0;
      int      j;

      for (j = 0; j < n_ys; j++)
      {
	// Store this point only if it is contained on the other arc.
	if (vertical_P->contains_point(xps[j]))
	{
	  ps[n] = xps[j];
	  n++;
	}
      }
      return (n);
    }

    // Solve a quadratic equation to find the x co-ordinates of the potential
    // intersection points:
    //
    //  (r*E^2 + s*B^2)*x^2 + (u*E^2 + 2*s*B*C - v*B*E)*x + 
    //                      + (w*E^2 + s*C^2 - v*C*E) = 0
    //
    NT       B, C, E;
    NT       xs[2];
    int      n_roots;

    if (arc._deg == 1)
    {
      B = arc._conic.u();
      C = arc._conic.w();
      E = arc._conic.v();
    }
    else if (_conic.s() == arc._conic.s())
    {
      // Both conics have the same orientation.
      B = _conic.u() - arc._conic.u();
      C = _conic.w() - arc._conic.w();
      E = _conic.v() - arc._conic.v();
    }
    else
    {
      // The two conics have opposite orientations.
      B = _conic.u() + arc._conic.u();
      C = _conic.w() + arc._conic.w();
      E = _conic.v() + arc._conic.v();      
    }

    n_roots = _solve_quadratic_eq
      (_conic.r()*E*E + _conic.s()*B*B,
       _conic.u()*E*E + 2*_conic.s()*B*C - _conic.v()*B*E,
       _conic.w()*E*E + _conic.s()*C*C - _conic.v()*C*E,
       xs);

    // Go over all roots, and return only those located on both arcs.
    int      n_ys;
    Point    xps[2];
    int      n = 0;
    int      i, j;

    for (i = 0; i < n_roots; i++)
    {
      // Find all points on our arc that have xs[i] as their x co-ordinate.
      n_ys = get_points_at_x (xs[i], xps);
      
      for (j = 0; j < n_ys; j++)
      {
	// Store this point only if it is contained on the other arc.
	if (arc.contains_point(xps[j]))
	{
	  ps[n] = xps[j];
	  n++;
	}
      }
    }

    // Return the number of intersection points.
    return (n);
  }

  // Check whether the two arcs overlap.
  // The function computes the number of overlapping arcs (2 at most), and
  // returns their number (0 means there is not overlap).
  int overlaps (const Segment_circle_2<NT>& arc,
		Segment_circle_2<NT>* ovlp_arcs) const
  {
    // Two arcs can overlap only if their base conics are identical.
    if (_conic != arc._conic)
      return (0);

    // In case one of the arcs is a full conic, return the whole other conic.
    if (arc.is_full_conic())
    {
      ovlp_arcs[0] = *this;
      return (1);
    }
    else if (is_full_conic())
    {
      ovlp_arcs[0] = arc;
      return (1);
    }

    // In case the other arc has an opposite orientation, switch its source
    // and target.
    const Point *arc_sourceP;
    const Point *arc_targetP;

    if (_conic.orientation() == arc._conic.orientation())
    {
      arc_sourceP = &(arc._source);
      arc_targetP = &(arc._target);
    }
    else
    {
      arc_sourceP = &(arc._target);
      arc_targetP = &(arc._source);
    }

    if (_is_strictly_between_endpoints(*arc_sourceP))
    {
      if (_is_strictly_between_endpoints(*arc_targetP))
      {
	// Check the next special case (when there are 2 overlapping arcs):
	if (arc._is_strictly_between_endpoints(_source) &&
            arc._is_strictly_between_endpoints(_target))
	{
	  ovlp_arcs[0] = Segment_circle_2<NT>(_conic, _source, *arc_targetP);
	  ovlp_arcs[1] = Segment_circle_2<NT>(_conic, *arc_sourceP, _target);
	  return (2);
	}

	// Case 1 - *this:     +----------->     
        //            arc:       +=====>
	ovlp_arcs[0] = Segment_circle_2<NT>(_conic, *arc_sourceP,*arc_targetP);
	return (1);
      }
      else
      {
	// Case 2 - *this:     +----------->     
        //            arc:               +=====>
	ovlp_arcs[0] = Segment_circle_2<NT>(_conic, *arc_sourceP, _target);
	return (1);
      }
    }
    else if (_is_strictly_between_endpoints(*arc_targetP))
    {
      // Case 3 - *this:     +----------->     
      //            arc:   +=====>
      ovlp_arcs[0] = Segment_circle_2<NT>(_conic, _source, *arc_targetP);
      return (1);
    }
    else if (arc._is_strictly_between_endpoints(_source) &&
             arc._is_strictly_between_endpoints(_target))
    {
      // Case 3 - *this:     +----------->     
      //            arc:   +================>
      ovlp_arcs[0] = *this;
      return (1);
    }
    
    // If we reached here, there are no overlaps:
    return (0);
  }
		
  private:

  // Check whether the given point is between the source and the target.
  // The point is assumed to be on the conic's boundary.
  bool _is_between_endpoints (const Point& p) const
  {
    if (p == _source || p == _target)
      return (true);
    else
      return (_is_strictly_between_endpoints(p));
  }
  
  // Check whether the given point is strictly between the source and the
  // target (but not any of them).
  // The point is assumed to be on the conic's boundary.
  bool _is_strictly_between_endpoints (const Point& p) const
  {
    // In case this is a full conic, any point on its boundary is between
    // its end points.
    if (is_full_conic())
      return (true);

    // In case the conic is a line segment (i.e. of the form ux + vy + w = 0),
    // make sure that p = (x,y) satisfies the segment equation (where (x1,y1)
    // is the source, (x2,y2) is the target and 0 < lambda < 1:
    //  (x,y) = (x1,y1) + lambda*(x2-x1,y2-y1)
    static const NT _zero = 0;
    static const NT _one = 1;

    if (_deg == 1)
    {
      NT   lambda;

      if (_source.x() != _target.x())
	lambda = (p.x() - _source.x()) / (_target.x() - _source.x());
      else
	lambda = (p.y() - _source.y()) / (_target.y() - _source.y());

      return (lambda > _zero && lambda < _one);
    }
    
    // Otherwise, make a decision based on the conic's orientation and whether
    // (source,p,target) is a right or a left turn.
    if (_conic.orientation() == 1)
      return (leftturn<R>(_source, p, _target));
    else
      return (rightturn<R>(_source, p, _target));
  }

  // Find the y-coordinates of the conic at a given x-coordinate.
  int _conic_get_y_coordinates (const NT& x,
                                NT *ys) const
  {
    // Solve the quadratic equation for a given x and find the y values:
    //  s*y^2 + (t*x + v)*y + (r*x^2 + u*x + w) = 0
    return (_solve_quadratic_eq (_conic.s(),
				 x*_conic.t() + _conic.v(),
				 x*(x*_conic.r() + _conic.u()) + _conic.w(),
				 ys));
  }

  // Find the x-coordinates of the conic at a given y-coordinate.
  int _conic_get_x_coordinates (const NT& y,
                                NT *xs) const
  {
    // Solve the quadratic equation for a given y and find the x values:
    //  r*x^2 + (t*y + u)*x + (s*y^2 + v*y + w) = 0
    return (_solve_quadratic_eq (_conic.r(),
				 y*_conic.t() + _conic.u(),
				 y*(y*_conic.s() + _conic.v()) + _conic.w(),
				 xs));
  }
  
  // Find the vertical tangency points of the conic.
  int _conic_vertical_tangency_points (Point* ps) const
  {
    // In case the base conic is of degree 1 (and not 2), the arc has no
    // vertical tangency points.
    static const NT _zero = 0;
    static const NT _two = 2;

    if (_deg == 1)
      return (0);

    // Find the vertical tangency points of the circle:
    NT              x0, y0, r;

    if (_conic.r() > _zero)
    {
      // Positive orientation:
      x0 = -(_conic.u() / _two);
      y0 = -(_conic.v() / _two);
      r = CGAL::sqrt(x0*x0 + y0*y0 - _conic.w());
    }
    else
    {
      // Negative orientation:
      x0 = _conic.u() / _two;
      y0 = _conic.v() / _two;
      r = CGAL::sqrt(x0*x0 + y0*y0 + _conic.w());
    }

    ps[0] = Point (x0 - r, y0);
    ps[1] = Point (x0 + r, y0);
    
    return (2);
  }

  // Find the horizontal tangency points of the conic.
  int _conic_horizontal_tangency_points (Point* ps) const
  {
    // In case the base conic is of degree 1 (and not 2), the arc has no
    // horizontal tangency points.
    static const NT _zero = 0;
    static const NT _two = 2;

    if (_deg == 1)
      return (0);

    // Find the horizontal tangency points of the circle:
    NT              x0, y0, r;

    if (_conic.r() > _zero)
    {
      // Positive orientation:
      x0 = -(_conic.u() / _two);
      y0 = -(_conic.v() / _two);
      r = CGAL::sqrt(x0*x0 + y0*y0 - _conic.w());
    }
    else
    {
      // Negative orientation:
      x0 = _conic.u() / _two;
      y0 = _conic.v() / _two;
      r = CGAL::sqrt(x0*x0 + y0*y0 + _conic.w());
    }

    ps[0] = Point (x0, y0 - r);
    ps[1] = Point (x0, y0 + r);
    
    return (2);
  }

};

#ifndef NO_OSTREAM_INSERT_CONIC_ARC_2
template <class NT>
std::ostream& operator<< (std::ostream& os, const Segment_circle_2<NT>& arc)
{
  Segment_circle_2<NT>::Conic conic = arc.conic();
  Segment_circle_2<NT>::Point source = arc.source(), target = arc.target();

  os << "{" << CGAL::to_double(conic.r()) << "*x^2 + "
     << CGAL::to_double(conic.s()) << "*y^2 + "
     << CGAL::to_double(conic.t()) << "*xy + " 
     << CGAL::to_double(conic.u()) << "*x + "
     << CGAL::to_double(conic.v()) << "*y + "
     << CGAL::to_double(conic.w()) << "}: "
     << "(" << CGAL::to_double(source.x()) << "," 
     << CGAL::to_double(source.y()) << ") -> "
     << "(" << CGAL::to_double(target.x()) << "," 
     << CGAL::to_double(target.y()) << ")";

  return (os);
}
#endif // NO_OSTREAM_INSERT_CONIC_ARC_2

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_CIRCLE_2_H
