// Copyright (c) 1999  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_HYPER_SEGMENT_2_H
#define CGAL_HYPER_SEGMENT_2_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Bbox_2.h>

#include <fstream>
class CGAL::Window_stream;

CGAL_BEGIN_NAMESPACE

/*!
 * Representation of a segment of a simplified canonical hyperbola, 
 * given by the equation:
 *    a*x^2 + b*y2 + c*x + d = 0
 *
 * We consider the upper part of the hyprobola (y > 0), thus the following
 * holds:
 *    y = sqrt(A*x^2 + B*x + C)
 *
 * The end-points of the hyperbolic segment are defined by x_min and x_max. 
 */
template <class Kernel_>
class Hyper_segment_2
{
protected:
  typedef Hyper_segment_2<Kernel_>      Self;
        
public:
  typedef Kernel_                       Kernel;
  typedef typename Kernel::FT           NT;
    
  typedef typename Kernel::Point_2      Point_2;

private:

  bool          _is_seg;      // Is this hyper-segment actaully a line segment.
  NT            _A;           // The coefficients of the equation of 
  NT            _B;           // the underlying canonical hyperbola:
  NT            _C;           //   y = sqrt (A*x^2 + B*x + C)

  Point_2       _source;      // The hyperbolic segment source.
  Point_2       _target;      // The hyperbolic segment target.

 public:

  /*!
   * Default constructor.
   */
  Hyper_segment_2 () :
    _is_seg(true)
  {}

  /*!
   * Construct a hyperbolic segment which lies on the canocial hyperbola:
   *    a*x^2 + b*y2 + c*x + d = 0
   * And bounded by x-min and x-max.
   * \param a The coefficient of x^2 in the equation of the hyperbola.
   * \param b The coefficient of y^2 in the equation of the hyperbola.
   * \param c The coefficient of x in the equation of the hyperbola.
   * \param d The free coefficient in the equation of the hyperbola.
   * \param x_min The x-coordinate of the leftmost end-point of the segment.
   * \param x_max The x-coordinate of the righttmost end-point of the segment.
   * \pre (a > 0) and (b < 0) -- to make sure this is indeed a hyperbola.
   *      Furthermore, x_max > x_min.
   */
  Hyper_segment_2 (const NT& a, const NT& b, const NT& c, const NT& d,
		   const NT& x_min, const NT& x_max) :
    _is_seg(false)
  {
    CGAL_precondition_code (static const NT  _zero = 0;);
    CGAL_precondition(_compare (a, _zero) == LARGER);
    CGAL_precondition(_compare (b, _zero) == SMALLER);
    CGAL_precondition(_compare (x_min, x_max) == SMALLER);

    // Set the normalized coefficients of the hyperbola.
    _A = -(a/b);
    _B = -(c/b);
    _C = -(d/b);

    // Set the end-points.
    bool  source_ok = _get_point_at_x (x_min, _source);
    bool  target_ok = _get_point_at_x (x_max, _target);

    if (! (source_ok && target_ok))
    {
      CGAL_assertion (source_ok);
      CGAL_assertion (target_ok);
    }

    // Check that both end-points lie on the same hyperbolic branch.
    /*
    CGAL_assertion_code (
      NT                xs[2];
      int               n = _solve_quadratic_equation (_A, _B, _C,
						       xs[0], xs[1]);
      Comparison_result res1;
      Comparison_result res2;
      int               i;

      CGAL_assertion (n == 2);
      for (i = 0; i < 2; i++)
      {
	res1 = _compare (x_min, xs[i]);
	res2 = _compare (x_max, xs[i]);

	CGAL_assertion (res1 == EQUAL || res2 == EQUAL || res1 == res2);
      }
    );
    */
  }

  /*!
   * Construct a degenerate hyperbolic segment defined by a line segment.
   * \param ps The source point of the segment.
   * \param pt The target point of the segment.
   * \pre The segment is not vertical -- that is, x(ps) does not equal x(pt).
   *      Furthermore, both points should lie above the x-axis.
   */
  Hyper_segment_2 (const Point_2& ps, const Point_2& pt) :
    _is_seg(true)
  {
    // Make sure that the two points do not define a vertical segment.
    Comparison_result comp_x = _compare (ps.x(), pt.x());

    CGAL_precondition(comp_x != EQUAL);

    // Make sure that both points are above the x-axis.
    CGAL_precondition_code (static const NT  _zero = 0;);
    CGAL_precondition(_compare (ps.y(), _zero) == LARGER);
    CGAL_precondition(_compare (pt.y(), _zero) == LARGER);

    // Find the line (y = a*x + b)  that connects the two points.
    const NT  denom = ps.x() - pt.x();
    const NT  a = (ps.y() - pt.y()) / denom;
    const NT  b = (ps.x()*pt.y() - pt.x()*ps.y()) / denom;

    // Set the underlying hyperbola to be the pair of lines:
    //  y^2 = (a*x + b)^2 = (a^2)*x^2 + (2ab)*x + b^2
    _A = a*a;
    _B = 2*a*b;
    _C = b*b;

    // Set the source and target point.
    if (comp_x == SMALLER)
    {
      _source = ps;
      _target = pt;
    }
    else
    {
      _source = pt;
      _target = ps;
    }
  }

  /*!
   * Access the hyperblic coeffcients.
   */
  const NT& A () const {return (_A);}
  const NT& B () const {return (_B);}
  const NT& C () const {return (_C);}

  /*!
   * Get the source of the hyper-segment.
   * \return The source point.
   */
  const Point_2& source () const
  {
    return (_source);
  }

  /*!
   * Get the target of the hyper-segment.
   * \return The target point.
   */
  const Point_2& target () const
  {
    return (_target);
  }

  /*! 
   * Check whether the hyper-segment is actually a line segment.
   * \return (true) if the underlying hyperbola is actually a line segment.
   */
  bool is_linear () const
  {
    return (_is_seg);
  }

  /*!
   * Check if the two hyper-segments are equal.
   * \param seg The compares hyper-segments.
   * \return (true) if seg equals (*this) or if it is its flipped version.
   */
  bool is_equal (const Self& seg) const
  {
    if (this == &seg)
      return (true);

    // The underlying hyperbola must be the same:
    if (! _has_same_base_hyperbola(seg))
      return (false);

    // Compare the source and target.
    if (_compare (_source.x(), seg._source.x()) &&
	_compare (_target.x(), seg._target.x()))
      return (true);
    else if (_compare (_source.x(), seg._target.x()) &&
	     _compare (_target.x(), seg._source.x()))
      return (true);
    
    return (false);
  }

  /*!
   * Get a bounding box for the hyper-segment.
   * \return A bounding box.
   */
  Bbox_2 bbox () const
  {
    // Use the source and target to find the exterme coordinates.
    double  x1 = CGAL::to_double(_source.x());
    double  y1 = CGAL::to_double(_source.y());
    double  x2 = CGAL::to_double(_target.x());
    double  y2 = CGAL::to_double(_target.y());
    double  x_min, x_max;
    double  y_min, y_max;

    if (x1 < x2)
    {
      x_min = x1;
      x_max = x2;
    }
    else
    {
      x_min = x2;
      x_max = x1;
    }

    if (y1 < y2)
    {
      y_min = y1;
      y_max = y2;
    }
    else
    {
      y_min = y2;
      y_max = y1;
    }

    // Return the resulting bounding box.
    return (Bbox_2 (x_min, y_min, x_max, y_max));
  }

  /*!
   * Check if the given point is in the x-range of the hyper-segment.
   * \param q The query point.
   * \return (true) if q is in the x-range of the hyper-segment.
   */
  bool point_is_in_x_range (const Point_2& q) const
  {
    Comparison_result res1 = _compare (q.x(), _source.x());
    Comparison_result res2 = _compare (q.x(), _target.x());
    
    return ((res1 == EQUAL) || (res2 == EQUAL) || (res1 != res2));
  }

  /*!
   * Check the position of the query point with respect to the hyper-segment.
   * \param q The query point.
   * \return SMALLER if q lies under the hyper-segment;
   *         LARGER if it lies above the hyper-segment;
   *         EQUAL if q lies on the hyper-segment.
   * \pre q is in the x-range of the hyper-segment.
   */
  Comparison_result point_position (const Point_2& q) const
  {
    CGAL_precondition(point_is_in_x_range (q)); 

    // Substitute q's x-coordinate into the equation of the hyperbola and
    // compare the result.
    return (_compare (q.y(),
		      CGAL::sqrt(_A*q.x()*q.x() + _B*q.x() + _C)));
  }

  /*!
   * Return a flipped hyper-segment.
   * \return The flipped hyper-segment.
   */
  Self flip () const
  {
    Self     flipped (*this);

    flipped._source = _target;
    flipped._target = _source;

    return (flipped);
  }

  /*!
   * Compare the y-coordinates of two hyper-segments at a given x-coordinate.
   * \param seg The other segment.
   * \param q The point (a placeholder for the x-coordinate).
   * \return SMALLER if (*this) is below seg at q;
   *         LARGER if it is above seg at q;
   *         EQUAL if the two hyper-segment intersect at this x-coordinate.
   * \pre q is in the x-range of both hyper-segments.
   */
  Comparison_result compare_y_at_x (const Self& seg,
				    const Point_2& q) const
  {
    CGAL_precondition(this->point_is_in_x_range(q));
    CGAL_precondition(seg.point_is_in_x_range(q));
  
    // Compare y1 = A1*x^2 + B1*x + C1 
    //     and y2 = A2*x^2 + B2*x + C2:
    return (_compare (_A*q.x()*q.x() + _B*q.x() + _C,
		      seg._A*q.x()*q.x() + seg._B*q.x() + seg._C));
  }

  /*!
   * Compare the slopes of two hyper-segments at their intersection point.
   * \param seg The other segment.
   * \param q A placeholder for the x-coordinate of the intersection point.
   * \return SMALLER if the slope of (*this) is less than seg's at q;
   *         LARGER if it is larger at q;
   *         EQUAL if the two slopes are equal (in case of an overlap).
   * \pre Both (*this) and seg are equal at q.x()
   */
  Comparison_result compare_slopes (const Self& seg,
				    const Point_2& q) const
  {
    // The underlying hyperbola is:
    //   y = sqrt(A*x^2 + B*x + C)
    //
    // If we derive, we obtain:
    //                 2A*x + B            2A*x + B
    //   y' = ------------------------- = ----------
    //         2*sqrt(A*x^2 + B*x + C)        y
    //
    NT                y2 = _A*q.x()*q.x() + _B*q.x() + _C;
    Comparison_result res = _compare (y2, 0);

    CGAL_precondition(_compare (y2,
				seg._A*q.x()*q.x() + seg._B*q.x() + seg._C) 
		      == EQUAL);

    CGAL_assertion(res != SMALLER);
    
    if (res == LARGER)
    {
      // The y-coordinate of the intersection point is greater than 0.
      // It is therefore sufficient to compare the numerators of the 
      // expressions of the two derivatives:
      res = _compare (2*_A*q.x() + _B,
		      2*seg._A*q.x() + seg._B);

      if (res != EQUAL)
	return (res);

      // Compare the second derivatives, where:
      //
      //          2A*y^2 - (2A*x + B)^2
      //   y'' = -----------------------
      //                2*y^3
      //
      return (_compare (2*_A*y2 - (2*_A*q.x() + _B)*(2*_A*q.x() + _B),
			2*seg._A*y2 - (2*seg._A*q.x() + seg._B)*
			(2*seg._A*q.x() + seg._B)));
    }
    else // (res == EQUAL), that is the y-coordinate is 0.
    {
      // Compare the second derivative by x, which in our case equals:
      return (_compare (2*_A*q.x() + _B,
			2*seg._A*q.x() + seg._B));
    }
  }

  /*!
   * Intersect the two hyper-segments.
   * \param seg The second hyper-segment.
   * \param p1 The output intersection point (if one exists).
   * \param p2 The second intersection point (in case of an overlap).
   * \return The number of intersection points computed:
   *         0 - No intersection.
   *         1 - A simple intersection returned as p1.
   *         2 - An overlapping segment, whose end-points are p1 and p2.
   */
  int intersect (const Self& seg,
		 Point_2& p1, Point_2& p2) const
  {
    // First check the case of overlapping hyper-segments.
    if (_has_same_base_hyperbola(seg))
    {
      bool bs = point_is_in_x_range (seg._source); 
      bool bt = point_is_in_x_range (seg._target); 
      
      if (bs && bt)
      {
	// (*this) contains the second segment:
	if (_compare (seg._source.x(), seg._target.x()) == SMALLER)
	{
	  p1 = seg._source;
	  p2 = seg._target;
	}
	else
	{
	  p2 = seg._source;
	  p1 = seg._target;
	}
	return (2);
      }
      else if (bs || bt)
      {
	// Only one of the end-points is contained in (*this):
	const Point_2& p_in  = bs ? seg._source : seg._target;
	const Point_2& p_out = bs ? seg._target : seg._source;
	bool to_right = (_compare(_source.x(), _target.x()) == SMALLER);
	Comparison_result res = _compare (p_out.x(), _source.x());

	if (_compare (p_in.x(), _source.x()) == EQUAL)
	{
	  // The two hyper-segments share just a common end-point.
	  if ((to_right && res == SMALLER) ||
	      (!to_right && res == LARGER))
	  {
	    p1 = p_in;
	    return (1);
	  }
	}
	else if (_compare (p_in.x(), _target.x()) == EQUAL)
	{
	  // The two hyper-segments share just a common end-point.
	  if ((to_right && res == LARGER) ||
	      (!to_right && res == SMALLER))
	  {
	    p1 = p_in;
	    return (1);
	  }
	}

	if (res == SMALLER)
	{
	  p1 = to_right ? _source : _target;
	  p2 = p_in;
	}
	else
	{
	  p1 = p_in;
	  p2 = to_right ? _target : _source;
	}
	return (2);
      }
      else
      {
	// Check if seg contains (*this).
	bs = seg.point_is_in_x_range (_source); 
	bt = seg.point_is_in_x_range (_target); 

	if (bs && bt)
	{
	  if (_compare (_source.x(), _target.x()) == SMALLER)
	  {
	    p1 = _source;
	    p2 = _target;
	  }
	  else
	  {
	    p2 = _source;
	    p1 = _target;
	  }
	  return (2); 
	}
      }

      // The two segments are disjoint:
      return (0);
    }

    // The equations of the two underlying hyperbola are given by:
    //   H1: y = A1*x^2 + B1*x + C1
    //   H2: y = A2*x^2 + B2*x + C2
    //
    // So the following must hold for the x-coordinates of the intersection
    // points:
    //   (A1 - A2)*x^2 + (B1 - B2)*x + (C1 - C2) = 0
    //
    NT     xs[2];
    bool   bs[2];
    int    n_xs;
    int    i;

    n_xs = _solve_quadratic_equation (_A - seg._A,
				      _B - seg._B,
				      _C - seg._C,
				      xs[0], xs[1]);

    // Check the legality of each of the x-coordinates.
    bs[0] = false;
    bs[1] = false;

    for (i = 0; i < n_xs; i++)
    {
      bs[i] = true;

      Comparison_result res1 = _compare(xs[i], _source.x());
      Comparison_result res2 = _compare(xs[i], _target.x());

      if (res1 != EQUAL && res2 != EQUAL && res1 == res2)
	bs[i] = false;
      else
      {
	res1 = _compare(xs[i], seg._source.x());
	res2 = _compare(xs[i], seg._target.x());
      
	if (res1 != EQUAL && res2 != EQUAL && res1 == res2)
	  bs[i] = false;
      }
    }

    // Make sure we return one point at most.
    if (!bs[0] && !bs[1])
      return (0);

    CGAL_assertion (!(bs[0] && bs[1]));

    if (bs[0])
      _get_point_at_x (xs[0], p1);
    else
      _get_point_at_x (xs[1], p1);

    return (1);
  }

  /*!
   * Split a hyper-segment at a given point.
   * \param p The split point.
   * \param seg1 The first resulting hyper-segment.
   * \param seg2 The second resulting hyper-segment.
   * \pre p lies in the interior of the curve (and is not and end-point).
   */
  void split (const Point_2& p,
	      Self& seg1, Self& seg2) const
  {
    CGAL_precondition(point_position(p) == EQUAL);
    CGAL_precondition(_compare(_source.x(), p.x()) != EQUAL);
    CGAL_precondition(_compare(_target.x(), p.x()) != EQUAL);

    seg1 = (*this);
    seg1._target = p;

    seg2 = (*this);
    seg2._source = p;

    return;
  }

  /*
   * Trim the hyper-segment at the given x-coodinates.
   * \param p1 A placeholder for the x-coordinates of the new source.
   * \param p2 A placeholder for the x-coordinates of the new target.
   * \return The trimmed hyper-segment.
   * \pre p1 and p2 are both in the x-range of the hyper-segment.
   */
  Self trim (const Point_2& p1,
	     const Point_2& p2) const
  {
    CGAL_precondition(point_is_in_x_range(p1));
    CGAL_precondition(point_is_in_x_range(p2));

    Self     trimmed (*this);

    _get_point_at_x (p1.x(), trimmed._source);
    _get_point_at_x (p2.x(), trimmed._target);

    return (trimmed);
  }

protected:

  /*!
   * Compare two values.
   */
  inline Comparison_result _compare (const NT& val1, const NT& val2) const
  {
    return (CGAL_NTS compare(val1, val2));
  }

  /*!
   * Check whether the underlying hyperbola contains a point with the given
   * x-coordinate and compute it if it does.
   * \param x The x-coordinate.
   * \param p The output point.
   * \return (true) if the underlying hyperbola contains a point with the 
   *         given x-coordinate.
   */
  bool _get_point_at_x (const NT& x,
                        Point_2& p) const
  {
    // Try to compute the y-coordinate.
    NT              val = _A*x*x + _B*x + _C;
    static const NT _zero = 0;

    if (_compare(val, _zero) == SMALLER)
      return (false);

    p = Point_2(x, CGAL::sqrt(val));
    return (true);
  }

  /*!
   * Check whether the two hyper-segments lie on the same hyperbola.
   * \param seg The compared hyper-segment.
   * \return (true) if the two hyper-segments lie on the same hyperbola.
   */
  bool _has_same_base_hyperbola (const Self& seg) const
  {
    return ((_compare (_A, seg._A) == EQUAL) &&
	    (_compare (_B, seg._B) == EQUAL) &&
	    (_compare (_C, seg._C) == EQUAL));
  }

  /*!
   * Solve the equation a*x^2 + b*x + c = 0.
   * \return The number of solutions found.
   */
  int _solve_quadratic_equation (const NT& a, const NT& b, const NT& c,
				 NT& x1, NT& x2) const
  {
    static const NT _zero = 0;
    static const NT _two = 2;
    static const NT _four = 4;

    if (_compare (a, _zero) == EQUAL)
    {
      // If this is a linear equation:
      if (_compare (b, _zero) == EQUAL)
	return (0);

      x1 = -c / b;
      return (1);
    }

    // Act according to the sign of (b^2 - 4ac):
    const NT          disc = b*b - _four*a*c;
    Comparison_result res = _compare (disc, _zero);

    if (res == SMALLER)
    {
      return (0);
    }
    else if (res == EQUAL)
    {
      x1 = (-b) / (_two*a);
      return (1);
    }
    
    const NT          sqrt_disc = CGAL::sqrt (disc);
    
    x1 = (sqrt_disc - b) / (_two*a);
    x2 = -(sqrt_disc + b) / (_two*a);
    return (2);
  }
};

#ifndef NO_OSTREAM_INSERT_CONIC_ARC_2
template <class Kernel>
std::ostream& operator<< (std::ostream& os, 
			  const Hyper_segment_2<Kernel>& seg)
{
  os << "{ y^2 = " << CGAL::to_double(seg.A()) << "*x^2 + "
     << CGAL::to_double(seg.B()) << "*x + " 
     << CGAL::to_double(seg.C()) << " } : "
     << "(" << CGAL::to_double(seg.source().x()) << "," 
     << CGAL::to_double(seg.source().y()) << ") -> "
     << "(" << CGAL::to_double(seg.target().x()) << "," 
     << CGAL::to_double(seg.target().y()) << ")";

  return (os);
}
#endif // NO_OSTREAM_INSERT_CONIC_ARC_2

CGAL_END_NAMESPACE

#endif
