// Copyright (c) 2005  Tel-Aviv University (Israel).
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

#ifndef CGAL_CONIC_ARC_2_H
#define CGAL_CONIC_ARC_2_H

/*! \file
 * Header file for the _Conic_arc_2<Int_kernel, Alg_kernel, Nt_traits> class.
 */

#include <CGAL/Arr_traits_2/Conic_point_2.h>
#include <CGAL/Bbox_2.h>

#include <ostream>

CGAL_BEGIN_NAMESPACE


/*!
 * Representation of a conic arc -- a bounded segment that lies on a conic
 * curve, the loci of all points satisfying the equation:
 *   r*x^2 + s*y^2 + t*xy + u*x + v*y +w = 0
 *
 * The class is templated with three parameters: 
 * Rat_kernel A kernel that provides the input objects or coefficients.
 *            Rat_kernel::FT should be an integral or a rational type.
 * Alg_kernel A geometric kernel, where Alg_kernel::FT is the number type
 *            for the coordinates of arrangement vertices, which are algebraic
 *            numbers of degree up to 4 (preferably it is CORE::Expr).
 * Nt_traits A traits class for performing various operations on the integer,
 *           rational and algebraic types. 
 */

template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_>
class _Conic_arc_2
{
public:

  typedef Rat_kernel_                                      Rat_kernel;
  typedef Alg_kernel_                                      Alg_kernel;
  typedef Nt_traits_                                       Nt_traits;

  typedef _Conic_arc_2<Rat_kernel, Alg_kernel, Nt_traits>  Self;
 
  typedef typename Rat_kernel::FT                          Rational;
  typedef typename Rat_kernel::Point_2                     Rat_point_2;
  typedef typename Rat_kernel::Segment_2                   Rat_segment_2;
  typedef typename Rat_kernel::Circle_2                    Rat_circle_2;

  typedef typename Alg_kernel::FT                          Algebraic;
  typedef typename Alg_kernel::Point_2                     Point_2;
  typedef _Conic_point_2<Alg_kernel>                       Conic_point_2;

protected:

  typedef typename Nt_traits::Integer                      Integer;

  Integer        _r;       //
  Integer        _s;       // The coefficients of the supporting conic curve:
  Integer        _t;       //
  Integer        _u;       //
  Integer        _v;       //   r*x^2 + s*y^2 + t*xy + u*x + v*y +w = 0 .
  Integer        _w;       //

  Orientation    _orient;  // The orientation of the conic.

  int            _info;    // Does the arc represent a full conic curve.
  Conic_point_2  _source;  // The source of the arc (if _info is 0).
  Conic_point_2  _target;  // The target of the arc (if _info is 0).

  // For arcs whose base is a hyperbola we store the axis (a*x + b*y + c = 0)
  // which separates the two bracnes of the hyperbola. We also store the side
  // (-1 or 1) that the arc occupies.
  struct Hyperbolic_arc_data
  {
    Algebraic     a;
    Algebraic     b;
    Algebraic     c;
    int           side;
  };

  Hyperbolic_arc_data *_hyper_data_P;

 public:

  /// \name Construction and destruction fucntions.
  //@{

  /*!
   * Default constructor.
   */
  _Conic_arc_2 () :
    _r(0), _s(0), _t(0), _u(0), _v(0), _w(0),
    _orient (COLLINEAR),
    _info (0),
    _hyper_data_P (NULL)
  {}

  /*!
   * Copy constructor.
   * \param arc The copied arc.
   */
  _Conic_arc_2 (const Self& arc) :
    _r(arc._r), _s(arc._s), _t(arc._t), _u(arc._u), _v(arc._v), _w(arc._w),
    _orient (arc._orient),
    _info (arc._info),
    _source(arc._source),
    _target(arc._target)
  {
    if (arc._hyper_data_P != NULL)
      _hyper_data_P = new Hyperbolic_arc_data (*(arc._hyper_data_P));
    else
      _hyper_data_P = NULL;
  }

  /*! 
   * Construct a conic arc which is the full conic:
   *   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
   * \pre The conic C must be an ellipse (so 4rs - t^2 > 0).
   */
  _Conic_arc_2 (const Rational& r, const Rational& s, const Rational& t,
                const Rational& u, const Rational& v, const Rational& w)
  {
    // Make sure the given curve is an ellipse.
    CGAL_precondition (CGAL_NTS compare(4*r*s - t*t, Rational(0)) == LARGER);

    // Set the arc to be the full conic (and compute the orientation).
    Rational    rat_coeffs [6];

    rat_coeffs[0] = r;
    rat_coeffs[1] = s;
    rat_coeffs[2] = t;
    rat_coeffs[3] = u;
    rat_coeffs[4] = v;
    rat_coeffs[5] = w;

    _set_full (rat_coeffs, true);
  }

  /*! 
   * Construct a conic arc which lies on the conic:
   *   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
   * \param orient The orientation of the arc (clockwise or counterclockwise).
   * \param source The source point.
   * \param target The target point.
   * \pre The source and the target must be on the conic boundary and must
   * not be the same.
   */
  _Conic_arc_2 (const Rational& r, const Rational& s, const Rational& t,
                const Rational& u, const Rational& v, const Rational& w,
                const Orientation& orient,
                const Point_2& source, const Point_2& target) :
    _orient (orient),
    _source (source),
    _target (target)
  {
    // Make sure that the source and the taget are not the same.
    CGAL_precondition (Alg_kernel().compare_xy_2_object() (source,
                                                           target) != EQUAL);

    // Set the arc properties (no need to compute the orientation).
    Rational    rat_coeffs [6];

    rat_coeffs[0] = r;
    rat_coeffs[1] = s;
    rat_coeffs[2] = t;
    rat_coeffs[3] = u;
    rat_coeffs[4] = v;
    rat_coeffs[5] = w;

    _set (rat_coeffs);
  }

  /*!
   * Construct a conic arc from the given line segment.
   * \param seg The line segment with rational endpoints.
   */
  _Conic_arc_2 (const Rat_segment_2& seg) :
    _orient (COLLINEAR)
  {
    // Set the source and target.
    Rat_kernel        ker;
    Rat_point_2       source = ker.construct_vertex_2_object() (seg, 0); 
    Rat_point_2       target = ker.construct_vertex_2_object() (seg, 1); 
    Rational          x1 = source.x();
    Rational          y1 = source.y();
    Rational          x2 = target.x();
    Rational          y2 = target.y();
    Nt_traits         nt_traits;

    _source = Point_2 (nt_traits.convert (x1), nt_traits.convert (y1));
    _target = Point_2 (nt_traits.convert (x2), nt_traits.convert (y2));
 
    // Make sure that the source and the taget are not the same.
    CGAL_precondition (Alg_kernel().compare_xy_2_object() (_source,
                                                           _target) != EQUAL);

    // The supporting conic is r=s=t=0, and u*x + v*y + w = 0 should hold
    // for both the source (x1,y1) and the target (x2, y2).
    const Rational    _zero (0);
    const Rational    _one (1);
    Rational          rat_coeffs [6];

    rat_coeffs[0] = _zero;
    rat_coeffs[1] = _zero;
    rat_coeffs[2] = _zero;

    if (CGAL_NTS compare (x1, x2) == EQUAL)
    {
      // The supporting conic is a vertical line, of the form x = CONST.
      rat_coeffs[3] = _one;
      rat_coeffs[4] = _zero;
      rat_coeffs[5] = -x1;
    }
    else
    {
      // The supporting line is A*x + B*y + C = 0, where:
      //
      //  A = y2 - y1,    B = x1 - x2,    C = x2*y1 - x1*y2 
      //
      rat_coeffs[3] = y2 - y1;
      rat_coeffs[4] = x1 - x2;
      rat_coeffs[5] = x2*y1 - x1*y2;
    }

    // Set the arc properties (no need to compute the orientation).
    _set (rat_coeffs);
  }

  /*!
   * Set a circular arc that corresponds to a full circle.
   * \param circ The circle (with rational center and rational squared radius).
   */ 
  _Conic_arc_2 (const Rat_circle_2& circ) :
    _orient (CLOCKWISE)
  {
    // Get the circle properties.
    Rat_kernel        ker;
    Rat_point_2       center = ker.construct_center_2_object() (circ); 
    Rational          x0 = center.x();
    Rational          y0 = center.y();
    Rational          R_sqr = ker.compute_squared_radius_2_object() (circ); 

    // Produce the correponding conic: if the circle center is (x0,y0)
    // and its squared radius is R^2, that its equation is:
    //   x^2 + y^2 - 2*x0*x - 2*y0*y + (x0^2 + y0^2 - R^2) = 0
    // Note that this equation describes a curve with a negative (clockwise) 
    // orientation.
    const Rational    _zero (0);
    const Rational    _one (1);
    const Rational    _minus_two (-2);
    Rational          rat_coeffs [6];

    rat_coeffs[0] = _one;
    rat_coeffs[1] = _one;
    rat_coeffs[2] = _zero;
    rat_coeffs[3] = _minus_two*x0;
    rat_coeffs[4] = _minus_two*y0;
    rat_coeffs[5] = x0*x0 + y0*y0 - R_sqr;

    // Set the arc to be the full conic (no need to compute the orientation).
    _set_full (rat_coeffs, false);
  }

  /*!
   * Set a circular arc that lies on the given circle:
   *   C: (x - x0)^2 + (y - y0)^2 = R^2
   * \param orient The orientation of the circle.
   * \param source The source point.
   * \param target The target point.
   * \pre The source and the target must be on the conic boundary and must
   *      not be the same.
   */  
  _Conic_arc_2 (const Rat_circle_2& circ,
                const Orientation& orient,
                const Point_2& source, const Point_2& target) :
    _orient(orient),
    _source(source),
    _target(target)
  {
    // Make sure that the source and the taget are not the same.
    CGAL_precondition (Alg_kernel().compare_xy_2_object() (source,
                                                           target) != EQUAL);
    CGAL_precondition (orient != COLLINEAR);

    // Get the circle properties.
    Rat_kernel        ker;
    Rat_point_2       center = ker.construct_center_2_object() (circ);
    Rational          x0 = center.x();
    Rational          y0 = center.y();
    Rational          R_sqr = ker.compute_squared_radius_2_object() (circ); 

    // Produce the correponding conic: if the circle center is (x0,y0)
    // and it squared radius is R^2, that its equation is:
    //   x^2 + y^2 - 2*x0*x - 2*y0*y + (x0^2 + y0^2 - R^2) = 0
    // Since this equation describes a curve with a negative (clockwise) 
    // orientation, we multiply it by -1 if necessary to obtain a positive
    // (counterclockwise) orientation.
    const Rational    _zero (0);
    Rational          rat_coeffs[6];

    if (_orient == COUNTERCLOCKWISE)
    {
      const Rational  _minus_one (-1);
      const Rational  _two (2);

      rat_coeffs[0] = _minus_one;
      rat_coeffs[1] = _minus_one;
      rat_coeffs[2] = _zero;
      rat_coeffs[3] = _two*x0;
      rat_coeffs[4] = _two*y0;
      rat_coeffs[5] = R_sqr - x0*x0 - y0*y0;
    }
    else
    {
      const Rational    _one (1);
      const Rational    _minus_two (-2);

      rat_coeffs[0] = _one;
      rat_coeffs[1] = _one;
      rat_coeffs[2] = _zero;
      rat_coeffs[3] = _minus_two*x0;
      rat_coeffs[4] = _minus_two*y0;
      rat_coeffs[5] = x0*x0 + y0*y0 - R_sqr;
    }

    // Set the arc properties (no need to compute the orientation).
    _set (rat_coeffs);
  }

  /*!
   * Construct a circular arc from the given three points.
   * \param p1 The arc source.
   * \param p2 A point in the interior of the arc.
   * \param p3 The arc target.
   * \pre The three points must not be collinear.
   */
  _Conic_arc_2 (const Rat_point_2& p1,
                const Rat_point_2& p2,
                const Rat_point_2& p3)
  {
    // Set the source and target.
    Rational          x1 = p1.x();
    Rational          y1 = p1.y();
    Rational          x2 = p2.x();
    Rational          y2 = p2.y();
    Rational          x3 = p3.x();
    Rational          y3 = p3.y();
    Nt_traits         nt_traits;

    _source = Point_2 (nt_traits.convert (x1), nt_traits.convert (y1));
    _target = Point_2 (nt_traits.convert (x3), nt_traits.convert (y3));
 
    // Make sure that the source and the taget are not the same.
    CGAL_precondition (Alg_kernel().compare_xy_2_object() (_source,
                                                           _target) != EQUAL);

    // Compute the lines: A1*x + B1*y + C1 = 0,
    //               and: A2*x + B2*y + C2 = 0,
    // where:
    const Rational  _two  = 2;

    const Rational  A1 = _two*(x1 - x2);
    const Rational  B1 = _two*(y1 - y2);
    const Rational  C1 = y2*y2 - y1*y1 + x2*x2 - x1*x1;

    const Rational  A2 = _two*(x2 - x3);
    const Rational  B2 = _two*(y2 - y3);
    const Rational  C2 = y3*y3 - y2*y2 + x3*x3 - x2*x2;

    // Compute the coordinates of the intersection point between the
    // two lines, given by (Nx / D, Ny / D), where:
    const Rational  Nx = B1*C2 - B2*C1;
    const Rational  Ny = A2*C1 - A1*C2;
    const Rational  D = A1*B2 - A2*B1;

    // Make sure the three points are not collinear.
    CGAL_precondition_code(
      const bool points_collinear = (CGAL::sign (D) == ZERO);
    );
    CGAL_precondition(!points_collinear);

    // The equation of the underlying circle is given by:
    Rational          rat_coeffs[6];
    
    rat_coeffs[0] = D*D;
    rat_coeffs[1] = D*D;
    rat_coeffs[2] = _zero;
    rat_coeffs[3] = -_two*D*Nx;
    rat_coeffs[4] = -_two*D*Ny;
    rat_coeffs[5] = 
      Nx*Nx + Ny*Ny - ((D*x2 - Nx)*(D*x2 - Nx) + (D*y2 - Ny)*(D*y2 - Ny));

    // Determine the orientation: If the mid-point forms a left-turn with
    // the source and the target points, the orientation is positive (going
    // counterclockwise).
    // Otherwise, it is negative (going clockwise).
    Alg_kernel                         ker;
    typename Alg_kernel::Orientation_2 orient_f = ker.orientation_2_object();


    Point_2  p_mid = Point_2 (nt_traits.convert (x2), nt_traits.convert (y2));
 
    if (orient_f(_source, p_mid, _target) == LEFT_TURN)
      _orient = COUNTERCLOCKWISE;
    else
      _orient = CLOCKWISE;

    // Set the arc properties (no need to compute the orientation).
    _set (rat_coeffs);
  }

  /*!
   * Destructor.
   */
  virtual ~_Conic_arc_2 ()
  {
    if (_hyper_data_P != NULL)
      delete _hyper_data_P;
  }

  /*!
   * Assignment operator.
   * \param arc The copied arc.
   */
  const Self& operator= (const Self& arc)
  {
    if (this == &arc)
      return (*this);

    // Free any existing data.
    if (_hyper_data_P != NULL)
      delete _hyper_data_P;

    // Copy the arc's attributes.
    _r = arc._r;
    _s = arc._s;
    _t = arc._t;
    _u = arc._u;
    _v = arc._v;
    _w = arc._w;

    _orient = arc._orient;
    _info = arc._info;
    _source = arc._source;
    _target = arc._target;

    // Duplicate the data for hyperbolic or circular arcs.
    if (arc._hyper_data_P != NULL)
      _hyper_data_P = new Hyperbolic_arc_data (*(arc._hyper_data_P));
    else
      _hyper_data_P = NULL;

    return (*this);
  }
  //@}

  /// \name Get the arc properties.
  //@{

  /*! 
   * Get the coefficients of the underlying conic.
   */
  const Integer& r () const {return (_r);}
  const Integer& s () const {return (_s);}
  const Integer& t () const {return (_t);}
  const Integer& u () const {return (_u);}
  const Integer& v () const {return (_v);}
  const Integer& w () const {return (_w);}

  /*!
   * Check whether the arc represents a full conic curve.
   */
  bool is_full_conic () const
  {
    return ((_info & 1) != 0);
  }

  /*!
   * Get the arc's source.
   * \return The source point.
   * \pre The arc does not represent a full conic curve.
   */
  const Point_2& source () const
  {
    CGAL_precondition (! is_full_conic());

    return (_source);
  }

  /*!
   * Get the arc's target.
   * \return The target point.
   * \pre The arc does not represent a full conic curve.
   */
  const Point_2& target () const
  {
    CGAL_precondition (! is_full_conic());

    return (_target);
  }

  /*!
   * Get the orientation of the arc.
   * \return The orientation.
   */
  Orientation orientation () const
  {
    return (_orient);
  }

  /*!
   * Get a bounding box for the conic arc.
   * \return The bounding box.
   */
  Bbox_2 bbox () const
  {
    double    x_min, y_min;
    double    x_max, y_max;

    if (is_full_conic())
    {
      // In case of a full conic (an ellipse or a circle), compute the
      // horizontal and vertical tangency points and use them to bound the arc.
      Point_2   tan_ps[2];
      int       n_tan_ps;

      n_tan_ps = vertical_tangency_points (tan_ps);
      CGAL_assertion (n_tan_ps == 2);

      if (CGAL::to_double(tan_ps[0].x()) < CGAL::to_double(tan_ps[1].x()))
      {
        x_min = CGAL::to_double(tan_ps[0].x());
        x_max = CGAL::to_double(tan_ps[1].x());
      }
      else
      {
        x_min = CGAL::to_double(tan_ps[1].x());
        x_max = CGAL::to_double(tan_ps[0].x());
      }

      n_tan_ps = horizontal_tangency_points (tan_ps);
      CGAL_assertion (n_tan_ps == 2);

      if (CGAL::to_double(tan_ps[0].y()) < CGAL::to_double(tan_ps[1].y()))
      {
        y_min = CGAL::to_double(tan_ps[0].y());
        y_max = CGAL::to_double(tan_ps[1].y());
      }
      else
      {
        y_min = CGAL::to_double(tan_ps[1].y());
        y_max = CGAL::to_double(tan_ps[0].y());
      }
    }
    else
    {
      // Use the source and target to initialize the exterme points.
      bool   source_left = 
        CGAL::to_double(_source.x()) < CGAL::to_double(_target.x());
      double x_min = source_left ?
        CGAL::to_double(_source.x()) : CGAL::to_double(_target.x());
      double x_max = source_left ?
        CGAL::to_double(_target.x()) : CGAL::to_double(_source.x());
      bool   source_down = 
        CGAL::to_double(_source.y()) < CGAL::to_double(_target.y());
      double y_min = source_down ?
        CGAL::to_double(_source.y()) : CGAL::to_double(_target.y());
      double y_max = source_down ?
        CGAL::to_double(_target.y()) : CGAL::to_double(_source.y());

      // Go over the vertical tangency points and try to update the x-points.
      Point_2    tan_ps[2];
      int        n_tan_ps;
      int        i;

      n_tan_ps = vertical_tangency_points (tan_ps);
      for (i = 0; i < n_tan_ps; i++)
      {
        if (CGAL::to_double(tan_ps[i].x()) < x_min)
          x_min = CGAL::to_double(tan_ps[i].x());
        if (CGAL::to_double(tan_ps[i].x()) > x_max)
          x_max = CGAL::to_double(tan_ps[i].x());
      }

      // Go over the horizontal tangency points and try to update the y-points.
      n_tan_ps = horizontal_tangency_points (tan_ps);
      for (i = 0; i < n_tan_ps; i++)
      {
        if (CGAL::to_double(tan_ps[i].y()) < y_min)
          y_min = CGAL::to_double(tan_ps[i].y());
        if (CGAL::to_double(tan_ps[i].y()) > y_max)
          y_max = CGAL::to_double(tan_ps[i].y());
      }
    }

    // Return the resulting bounding box.
    return (Bbox_2 (x_min, y_min, x_max, y_max));
  }
  //@}

  /// \name Compute points on the arc.
  //@{

  /*!
   * Calculate the vertical tangency points of the arc.
   * \param vpts The vertical tangency points.
   * \pre The vpts vector should be allocated at the size of 2.
   * \return The number of vertical tangency points.
   */
  int vertical_tangency_points (Point_2* vpts) const
  {
    // No vertical tangency points for line segments:
    if (_orient == COLLINEAR)
      return (0);

    // Calculate the vertical tangency points of the supporting conic.
    Point_2 ps[2];
    int     n;

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

  /*!
   * Calculate the horizontal tangency points of the arc.
   * \param hpts The horizontal tangency points.
   * \pre The hpts vector should be allocated at the size of 2.
   * \return The number of horizontal tangency points.
   */
  int horizontal_tangency_points (Point_2* hpts) const
  {
    // No horizontal tangency points for line segments:
    if (_orient == COLLINEAR)
      return (0);

    // Calculate the horizontal tangency points of the conic.
    Point_2    ps[2];
    int        n;

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

    // Return the number of horizontal tangency points found.
    return (m);
  }

  /*!
   * Find all points on the arc with a given x-coordinate.
   * \param p A placeholder for the x-coordinate.
   * \param ps The point on the arc at x(p).
   * \pre The vector ps should be allocated at the size of 2.
   * \return The number of points found.
   */
  int get_points_at_x (const Point_2& p,
                       Point_2 *ps) const
  {
    // Get the y coordinates of the points on the conic.
    Algebraic    ys[2];
    int          n;

    n = _conic_get_y_coordinates (p.x(), ys);
    
    // Find all the points that are contained in the arc.
    int   m = 0;
    
    for (int i = 0; i < n; i++)
    {
      ps[m] = Point_2 (p.x(), ys[i], _conic_id);

      if (is_full_conic() || _is_between_endpoints(ps[m]))
        m++;
    }

    // Return the number of points on the arc.
    return (m);
  }

  /*!
   * Find all points on the arc with a given y-coordinate.
   * \param p A placeholder for the y-coordinate.
   * \param ps The point on the arc at x(p).
   * \pre The vector ps should be allocated at the size of 2.
   * \return The number of points found.
   */
  int get_points_at_y (const Point_2& p,
                       Point_2 *ps) const
  {
    // Get the y coordinates of the points on the conic.
    Algebraic    xs[2];
    int          n;

    n = _conic_get_x_coordinates (p.y(), xs);
    
    // Find all the points that are contained in the arc.
    int   m = 0;
    
    for (int i = 0; i < n; i++)
    {
      ps[m] = Point_2 (xs[i], p.y(), _conic_id);

      if (is_full_conic() || _is_between_endpoints(ps[m]))
        m++;
    }

    // Return the number of points on the arc.
    return (m);
  }
  //@}

private:

  /// \name Auxiliary construction functions.
  //@{

  /*!
   * Set the properties of a conic arc (for the usage of the constructors).
   * \param rat_coeffs A vector of size 6, storing the rational coefficients
   *                   of x^2, y^2, xy, x, y and the free coefficient resp.
   */
  void _set (const Rational* rat_coeffs)
  {
    // Convert the coefficients vector to an equivalent vector of integer
    // coefficients.
    Nt_traits         nt_traits;
    Integer           int_coeffs[6];

    nt_traits.convert_coefficients (rat_coeffs, rat_coeffs + 6,
                                    int_coeffs);

    // Check the orientation of conic curve, and negate the conic coefficients
    // if its given orientation.
    typename Rat_kernel::Conic_2   temp_conic (rat_coeffs[0], rat_coeffs[1],
                                               rat_coeffs[2], rat_coeffs[3],
                                               rat_coeffs[4], rat_coeffs[5]);

    if (_orient == temp_conic.orientation())
    {
      _r = int_coeffs[0];
      _s = int_coeffs[1];
      _t = int_coeffs[2];
      _u = int_coeffs[3];
      _v = int_coeffs[4];
      _w = int_coeffs[5];
    }
    else
    {
      _r = -int_coeffs[0];
      _s = -int_coeffs[1];
      _t = -int_coeffs[2];
      _u = -int_coeffs[3];
      _v = -int_coeffs[4];
      _w = -int_coeffs[5];
    }

    // In case the base conic is a hyperbola, build the hyperbolic data
    // (this happens when (4rs - t^2) < 0).
    if (CGAL::sign (4*_r*_s - _t*_t) == NEGATIVE)
      _build_hyperbolic_arc_data ();
    else
      hyper_data_P = NULL;

    // Mark that this arc is not a full conic curve.
    _info = 0;

    // In case of a non-degenerate parabola or a hyperbola, make sure 
    // the arc is not infinite.
    CGAL_precondition_code (
      if ((CGAL::sign (_r) != ZERO ||
           CGAL::sign (_s) != ZERO ||
           CGAL::sign (_t) != ZERO) &&
          CGAL::sign (4*_r*_s - _t*_t) != POSITIVE)
      {
        Alg_kernel       ker;
        Point_2          p_mid = ker.construct_midpoint_2_object() (_source,
                                                                    _target);
        Point_2          ps[2];

        bool  finite_at_x = (get_points_at_x(p_mid, ps) > 0);
        bool  finite_at_y = (get_points_at_y(p_mid, ps) > 0);
      
        CGAL_precondition(finite_at_x && finite_at_y);
      }
    );

    return;
  }

  /*!
   * Set the properties of a conic arc that is really a full curve
   * (that is, an ellipse).
   * \param rat_coeffs A vector of size 6, storing the rational coefficients
   *                   of x^2, y^2, xy, x, y and the free coefficient resp.
   * \param comp_orient Should we compute the orientation of the given curve.
   */
  void _set_full (const Rational* rat_coeffs,
                  const bool& comp_orient)
  {
    // Convert the coefficients vector to an equivalent vector of integer
    // coefficients.
    Nt_traits         nt_traits;
    Integer           int_coeffs[6];

    nt_traits.convert_coefficients (rat_coeffs, rat_coeffs + 6,
                                    int_coeffs);

    // Check the orientation of conic curve, and negate the conic coefficients
    // if its given orientation.
    typename Rat_kernel::Conic_2   temp_conic (rat_coeffs[0], rat_coeffs[1],
                                               rat_coeffs[2], rat_coeffs[3],
                                               rat_coeffs[4], rat_coeffs[5]);
    const Orientation              temp_orient = temp_conic.orientation();

    if (comp_orient)
      _orient = temp_orient;

    if (_orient == temp_orient)
    {
      _r = int_coeffs[0];
      _s = int_coeffs[1];
      _t = int_coeffs[2];
      _u = int_coeffs[3];
      _v = int_coeffs[4];
      _w = int_coeffs[5];
    }
    else
    {
      _r = -int_coeffs[0];
      _s = -int_coeffs[1];
      _t = -int_coeffs[2];
      _u = -int_coeffs[3];
      _v = -int_coeffs[4];
      _w = -int_coeffs[5];
    }

    // Make sure the conic is a non-degenerate ellipse:
    // The coefficients should satisfy (4rs - t^2) > 0.
    CGAL_precondition (CGAL::sign (4*_r*_s - _t*_t) == POSITIVE);

    _hyper_data_P = NULL;

    // Mark that this arc is a full conic curve.
    _info = 1;

    return;
  }

  /*!
   * Build the data for hyperbolic arc, contaning the characterization of the
   * hyperbolic branch the arc is placed on.
   */
  void _build_hyperbolic_arc_data ()
  {
    // Let phi be the rotation angle of the conic from its canonic form.
    // We can write:
    // 
    //                          t
    //  sin(2*phi) = -----------------------
    //                sqrt((r - s)^2 + t^2)
    // 
    //                        r - s
    //  cos(2*phi) = -----------------------
    //                sqrt((r - s)^2 + t^2)
    //
    Nt_traits        nt_traits;
    const int        or_fact = (_orient == CLOCKWISE) ? -1 : 1;
    const Algebraic  r = nt_traits.convert (or_fact * _r);
    const Algebraic  s = nt_traits.convert (or_fact * _s);
    const Algebraic  t = nt_traits.convert (or_fact * _t);
    const Algebraic  cos_2phi = (r - s) / nt_traits.sqrt((r-s)*(r-s) + t*t);
    const Algebraic  _one = 1;
    const Algebraic  _two = 2;
    Algebraic        sin_phi;
    Algebraic        cos_phi;

    // Calculate sin(phi) and cos(phi) according to the half-angle formulae:
    // 
    //  sin(phi)^2 = 0.5 * (1 - cos(2*phi))
    //  cos(phi)^2 = 0.5 * (1 + cos(2*phi))
    Sign             sign_t = CGAL::sign (t);

    if (sign_t == ZERO)
    {
      // sin(2*phi) == 0, so phi = 0 or phi = PI/2
      if (CGAL::sign (cos_2phi) == POSITIVE)
      {
        // phi = 0.
        sin_phi = _zero;
        cos_phi = _one;
      }
      else
      {
        // phi = PI.
        sin_phi = _zero;
        cos_phi = -_one;
      }
    }
    else if (sign_t == POSITIVE)
    {
      // sin(2*phi) > 0 so 0 < phi < PI/2.
      sin_phi = nt_traits.sqrt((_one + cos_2phi) / _two);
      cos_phi = nt_traits.sqrt((_one - cos_2phi) / _two);
    }
    else
    {
      // sin(2*phi) < 0 so PI/2 < phi < PI.
      sin_phi = nt_traits.sqrt((_one + cos_2phi) / _two);
      cos_phi = -nt_traits.sqrt((_one - cos_2phi) / _two);
    }
                               
    // Calculate the center (x0, y0) of the conic, given by the formulae:
    //
    //        t*v - 2*s*u                t*u - 2*r*v
    //  x0 = -------------   ,     y0 = -------------
    //        4*r*s - t^2                4*r*s - t^2
    //
    // The denominator (4*r*s - t^2) must be negative for hyperbolas.
    const Algebraic  u = nt_traits.convert (or_fact * _u);
    const Algebraic  v = nt_traits.convert (or_fact * _v);
    const Algebraic  det = 4*r*s - t*t;
    Algebraic        x0, y0;

    CGAL_assertion (CGAL::sign (det) == NEGATIVE);
    
    x0 = (t*v - _two*s*u) / det;
    y0 = (t*u - _two*r*v) / det;
    
    // The axis separating the two branches of the hyperbola is now given by:
    // 
    //  cos(phi)*x + sin(phi)*y - (cos(phi)*x0 + sin(phi)*y0) = 0
    //
    _hyper_data_P = new Hyperbolic_arc_data;

    _hyper_data_P->a = cos_phi;
    _hyper_data_P->b = sin_phi;
    _hyper_data_P->c = - (cos_phi*x0 + sin_phi*y0);

    // Make sure that the two endpoints are located on the same branch
    // of the hyperbola.
    _hyper_data_P->side = _hyperbolic_arc_side(_source);

    CGAL_assertion (_hyper_data_P->side = _hyperbolic_arc_side(_target));

    return;
  }
  //@}

protected:

  /// \name Auxiliary functions.
  //@{

  /*!
   * Find on which branch of the hyperbola is the given point located.
   * The point is assumed to be on the hyperbola.
   * \param p The query point.
   * \return The branch ID (either -1 or 1).
   */
  int _hyperbolic_arc_side (const Point_2& p) const
  {
    if (_hyper_data_P == NULL)
      return (0);

    Algebraic         val = (_hyper_data_P->a*p.x() + _hyper_data_P->b*p.y() + 
                             _hyper_data_P->c);
    Sign              sign_val = CGAL::sign (val);

    CGAL_assertion (sign_val != ZERO);
    return ((sign_val == POSITIVE) ? 1 : -1);
  }
 
  /*!
   * Check whether the given point is between the source and the target.
   * The point is assumed to be on the conic's boundary.
   * \param p The query point.
   * \return (true) if the point is between the two endpoints, 
   *         (false) if it is not.
   */
  bool _is_between_endpoints (const Point_2& p) const
  {
    CGAL_precondition (! is_full_conic());

    // Check if p is one of the endpoints.
    Alg_kernel                         ker;

    if (ker.equal_2_object() (p, _source) || 
        ker.equal_2_object() (p, _target))
    {
      return (true);
    }
    else
    {
      return (_is_strictly_between_endpoints(p));
    }
  }

  /*!
   * Check whether the given point is strictly between the source and the
   * target (but not any of them).
   * The point is assumed to be on the conic's boundary.
   * \param p The query point.
   * \return (true) if the point is strictly between the two endpoints, 
   *         (false) if it is not.
   */
  bool _is_strictly_between_endpoints (const Point_2& p) const
  {
    // In case this is a full conic, any point on its boundary is between
    // its end points.
    if (is_full_conic())
      return (true);

    // In case of a hyperbolic arc, make sure the point is located on the
    // same branch as the arc.
    if (_hyper_data_P != NULL)
    {
      if (_hyperbolic_arc_side(p) != _hyper_data_P->side)
        return (false);
    }

    // Act according to the conic degree.
    Alg_kernel                         ker;

    if (_orient == COLLINEAR)
    {
      Comparison_result  res1;
      Comparison_result  res2;

      if (ker.compare_x_2_object() (_source, _target) == EQUAL)
      {
        // In case of a vertical segment - just check whether the y coordinate
        // of p is between those of the source's and of the target's.
        res1 = ker.compare_y_2_object() (p, _source);
        res2 = ker.compare_y_2_object() (p, _target);
      }
      else
      {
        // Otherwise, since the segment is x-monotone, just check whether the
        // x coordinate of p is between those of the source's and of the 
        // target's.
        res1 = ker.compare_x_2_object() (p, _source);
        res2 = ker.compare_x_2_object() (p, _target);
      }

      return (res1 != EQUAL && res2 != EQUAL && res1 != res2);
    }
    else
    {
      // In case of a conic of degree 2, make a decision based on the conic's
      // orientation and whether (source,p,target) is a right or a left turn.
      if (_orient == COUNTERCLOCKWISE)
        return (ker.orientation_2_object()(_source, p, _target) == LEFT_TURN);
      else
        return (ker.orientation_2_object()(_source, p, _target) == RIGHT_TURN);
    }
  }
  
  /*!
   * Find the vertical tangency points of the undelying conic.
   * \param ps The output points of vertical tangency.
   *           This area must be allocated at the size of 2.
   * \return The number of vertical tangency points.
   */
  int _conic_vertical_tangency_points (Point_2* ps) const
  {
    // In case the base conic is of degree 1 (and not 2), the arc has no
    // vertical tangency points.
    if (CGAL::sign (_s) == ZERO)
      return (0);

    // We are interested in the x coordinates where the quadratic equation:
    //  s*y^2 + (t*x + v)*y + (r*x^2 + u*x + w) = 0
    // has a single solution (obviously if s = 0, there are no such points).
    // We therefore demand that the discriminant of this equation is zero:
    //  (t*x + v)^2 - 4*s*(r*x^2 + u*x + w) = 0
    const Integer _two = 2;
    const Integer _four = 4;
    Algebraic     xs[2];
    Algebraic    *xs_end;
    int           n_xs;
    Nt_traits     nt_traits;
    
    xs_end = nt_traits.solve_quadratic_equation (_t*_t - _four*_r*_s,
						 _two*_t*_v - _four*_s*_u,
						 _v*_v - _four*_s*_w,
						 xs);
    n_xs = xs_end - xs;

    // Find the y-coordinates of the vertical tangency points.
    Algebraic     ys[2];
    Algebraic    *ys_end;
    int           n_ys;

    if (CGAL::sign (_t) == ZERO)
    {
      // The two vertical tangency points have the same y coordinate:
      ys[0] = nt_traits.convert (-_v) /nt_traits.convert (_two*_s);
      n_ys = 1;
    }
    else
    {
      ys_end =
        nt_traits.solve_quadratic_equation (_four*_r*_s*_s - _s*_t*_t,
                                            _four*_r*_s*_v - _two*_s*_t*_u,
                                            _r*_v*_v - _t*_u*_v + _t*_t*_w,
                                            ys);
      n_ys = ys_end - ys;
    }

    // Pair the x and y coordinates and obtain the vertical tangency points.
    int   n = 0;
    int   i, j;

    for (i = 0; i < n_xs; i++)
    {
      if (n_ys == 1)
      {
        ps[n] = Point_2 (xs[i], ys[0]);
        n++;
      }
      else
      {
        for (j = 0; j < n_ys; j++)
        {
          if (CGAL::compare (ys[j], 
                             -(nt_traits.convert(_t) * xs[i] + 
                               nt_traits.convert(_v)) /
                             nt_traits.convert(_two*_s)) == EQUAL)
          {
            ps[n] = Point_2 (xs[i], ys[j]);
            n++;
            break;
          }
        }
      }
    }

    return (n);
  }

  /*!
   * Find the horizontal tangency points of the undelying conic.
   * \param ps The output points of horizontal tangency.
   *           This area must be allocated at the size of 2.
   * \return The number of horizontal tangency points.
   */
  int _conic_horizontal_tangency_points (Point_2* ps) const
  {
    const Integer _zero = 0;

    // In case the base conic is of degree 1 (and not 2), the arc has no
    // vertical tangency points.
    if (CGAL::sign (_r) == ZERO)
      return (0);

    // We are interested in the y coordinates were the quadratic equation:
    //  r*x^2 + (t*y + u)*x + (s*y^2 + v*y + w) = 0
    // has a single solution (obviously if r = 0, there are no such points).
    // We therefore demand that the discriminant of this equation is zero:
    //  (t*y + u)^2 - 4*r*(s*y^2 + v*y + w) = 0
    const Integer _two = 2;
    const Integer _four = 4;
    int           n;
    Algebraic     ys[2];
    Algebraic    *ys_end;
    Nt_traits     nt_traits;

    ys_end = nt_traits.solve_quadratic_equation (_t*_t - _four*_r*_s,
						 _two*_t*_u - _four*_r*_v,
						 _u*_u - _four*_r*_w,
						 ys);
    n = ys_end - ys;

    // Compute the x coordinates and construct the horizontal tangency points.
    Algebraic     x;
    int           i;

    for (i = 0; i < n; i++)
    {
      // Having computed y, x is the simgle solution to the quadratic equation
      // above, and since its discriminant is 0, x is simply given by:
      x = -(nt_traits.convert(_t)*ys[i] + nt_traits.convert(_u)) / 
        nt_traits.cnvert(_two*_r);

      ps[i] = Point_2 (x, ys[i]);
    }
      
    return (n);
  }

  /*!
   * Find the y coordinates of the underlying conic at a given x coordinate.
   * \param x The x coordinate.
   * \param ys The output y coordinates. 
   * \pre The vector ys must be allocated at the size of 2.
   * \return The number of y coordinates computed (either 0, 1 or 2).
   */
  int _conic_get_y_coordinates (const Algebraic& x,
                                Algebraic *ys) const
  {
    // Solve the quadratic equation for a given x and find the y values:
    //  s*y^2 + (t*x + v)*y + (r*x^2 + u*x + w) = 0
    Nt_traits     nt_traits;
    Algebraic     A = nt_traits.convert(_s);
    Algebraic     B = nt_traits.convert(_t)*x + nt_traits.convert(_v);
    Algebraic     C = (nt_traits.convert(_r)*x + 
                       nt_traits.convert(_u))*x + nt_traits.convert(_w);

    return (_solve_quadratic_equation (A, B, C, ys[0], ys[1]));
  }

  /*!
   * Find the x coordinates of the underlying conic at a given y coordinate.
   * \param y The y coordinate.
   * \param xs The output x coordinates. 
   * \pre The vector xs must be allocated at the size of 2.
   * \return The number of x coordinates computed (either 0, 1 or 2).
   */
  int _conic_get_x_coordinates (const Algebraic& y,
                                Algebraic *xs) const
  {
    // Solve the quadratic equation for a given y and find the x values:
    //  r*x^2 + (t*y + u)*x + (s*y^2 + v*y + w) = 0
    Nt_traits     nt_traits;
    Algebraic     A = nt_traits.convert(_r);
    Algebraic     B = nt_traits.convert(_t)*y + nt_traits.convert(_u);
    Algebraic     C = (nt_traits.convert(_s)*y + 
                       nt_traits.convert(_v))*y + nt_traits.convert(_w);

    return (_solve_quadratic_equation (A, B, C, xs[0], xs[1]));
  }

  /*!
   * Solve the given quadratic equation: Ax^2 + B*x + C = 0.
   * \param x_minus The root obtained from taking -sqrt(discriminant).
   * \param x_plus The root obtained from taking -sqrt(discriminant).
   * \return The number of disticnt solutions to the equation.
   */
  int _solve_quadratic_equation (const Algebraic& A,
                                 const Algebraic& B,
                                 const Algebraic& C,
                                 Algebraic& x_minus, Algebraic& x_plus) const
  {
    const Algebraic  disc = B*B - 4*A*C;
    Sign             sign_disc = CGAL::sign (disc);

    if (sign_disc == NEGATIVE)
    {
      // No real-valued solutions:
      return (0);
    }
    else if (sign_disc == ZERO)
    {
      // One distinct solution:
      x_minus = x_plus = -B / (2*A);
      return (1);
    }

    // Compute the two distinct solutions:
    Algebraic     _2A = 2*A;
    Nt_traits     nt_traits;
    Algebraic     sqrt_disc = nt_traits.sqrt (disc);

    x_minus = -(B + sqrt_disc) / _2A;
    x_plus = (sqrt_disc - B) / _2A;
    return (2);
  }
  //@}

};

/*!
 * Exporter for conic arcs.
 */
template <class Rat_kernel, class Alg_kernel, class Nt_traits>
std::ostream& 
operator<< (std::ostream& os, 
            const _Conic_arc_2<Rat_kernel, Alg_kernel, Nt_traits> & arc)
{
  os << "{" << CGAL::to_double(arc.r()) << "*x^2 + "
     << CGAL::to_double(arc.s()) << "*y^2 + "
     << CGAL::to_double(arc.t()) << "*xy + " 
     << CGAL::to_double(arc.u()) << "*x + "
     << CGAL::to_double(arc.v()) << "*y + "
     << CGAL::to_double(arc.w()) << "}";

  if (arc.is_full_conic())
  {
    os << " - Full curve";
  }
  else
  {
    os << " : (" << CGAL::to_double(arc.source().x()) << "," 
       << CGAL::to_double(arc.source().y()) << ") ";

    if (arc.orientation() == CLOCKWISE)
      os << "--cw-->";
    else if (arc.orientation() == COUNTERCLOCKWISE)
      os << "--ccw-->";
    else
      os << "--l-->";

    os << " (" << CGAL::to_double(arc.target().x()) << "," 
       << CGAL::to_double(arc.target().y()) << ")";
  }

  return (os);
}

CGAL_END_NAMESPACE

#endif
