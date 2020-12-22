// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_CONIC_ARC_2_H
#define CGAL_CONIC_ARC_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Header file for the _Conic_arc_2<Int_kernel, Alg_kernel, Nt_traits> class.
 */

#include <CGAL/Arr_geometry_traits/Conic_point_2.h>
#include <CGAL/Arr_geometry_traits/Conic_intersections_2.h>
#include <CGAL/Bbox_2.h>

#include <ostream>

namespace CGAL {


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

  typedef typename Nt_traits::Integer                      Integer;

  typedef typename Alg_kernel::FT                          Algebraic;
  typedef typename Alg_kernel::Point_2                     Point_2;
  typedef _Conic_point_2<Alg_kernel>                       Conic_point_2;

protected:

  Integer        _r;       //
  Integer        _s;       // The coefficients of the supporting conic curve:
  Integer        _t;       //
  Integer        _u;       //
  Integer        _v;       //   r*x^2 + s*y^2 + t*xy + u*x + v*y +w = 0 .
  Integer        _w;       //

  Orientation    _orient;  // The orientation of the conic.

  // Bit masks for the _info field.
  enum
  {
    IS_VALID = 1,
    IS_FULL_CONIC = 2
  };

  int            _info;    // Does the arc represent a full conic curve.
  Conic_point_2  _source;  // The source of the arc (if not a full curve).
  Conic_point_2  _target;  // The target of the arc (if not a full curve).

  /*! \struct
   * For arcs whose base is a hyperbola we store the axis (a*x + b*y + c = 0)
   * which separates the two bracnes of the hyperbola. We also store the side
   * (NEGATIVE or POSITIVE) that the arc occupies.
   * In case of line segments connecting two algebraic endpoints, we use this
   * structure two store the coefficients of the line supporting this segment.
   * In this case we set the side field to be ZERO.
   */
  struct Extra_data
  {
    Algebraic     a;
    Algebraic     b;
    Algebraic     c;
    Sign          side;
  };

  Extra_data    *_extra_data_P;  // The extra data stored with the arc
                                 // (may be nullptr).

public:

  /// \name Construction and destruction functions.
  //@{

  /*!
   * Default constructor.
   */
  _Conic_arc_2 () :
    _r(0), _s(0), _t(0), _u(0), _v(0), _w(0),
    _orient (COLLINEAR),
    _info (0),
    _extra_data_P (nullptr)
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
    if (arc._extra_data_P != nullptr)
      _extra_data_P = new Extra_data (*(arc._extra_data_P));
    else
      _extra_data_P = nullptr;
  }

  /*!
   * Construct a conic arc which is the full conic:
   *   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
   * \pre The conic C must be an ellipse (so 4rs - t^2 > 0).
   */
  _Conic_arc_2 (const Rational& r, const Rational& s, const Rational& t,
                const Rational& u, const Rational& v, const Rational& w) :
    _extra_data_P (nullptr)
  {
    // Make sure the given curve is an ellipse (4rs - t^2 should be positive).
    CGAL_precondition (CGAL::sign (4*r*s - t*t) == POSITIVE);

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
    _target (target),
    _extra_data_P (nullptr)
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
    _orient (COLLINEAR),
    _extra_data_P (nullptr)
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

    if (CGAL::compare (x1, x2) == EQUAL)
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
    _orient (CLOCKWISE),
    _extra_data_P (nullptr)
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
    _target(target),
    _extra_data_P (nullptr)
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
                const Rat_point_2& p3):
    _extra_data_P (nullptr)
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
    const bool  points_collinear = (CGAL::sign (D) == ZERO);

    if (points_collinear)
    {
      _info = 0;           // Inavlid arc.
      return;
    }

    // The equation of the underlying circle is given by:
    Rational          rat_coeffs[6];

    rat_coeffs[0] = D*D;
    rat_coeffs[1] = D*D;
    rat_coeffs[2] = 0;
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
   * Construct a conic arc from the given five points, specified by the
   * points p1, p2, p3, p4 and p5.
   * \param p1 The source point of the given arc.
   * \param p2,p3,p4 Points lying on the conic arc, between p1 and p5.
   * \param p5 The target point of the given arc.
   * \pre No three points are collinear.
   */
  _Conic_arc_2 (const Rat_point_2& p1,
                const Rat_point_2& p2,
                const Rat_point_2& p3,
                const Rat_point_2& p4,
                const Rat_point_2& p5) :
    _extra_data_P(nullptr)
  {
    // Make sure that no three points are collinear.
    Rat_kernel                         ker;
    typename Rat_kernel::Orientation_2 orient_f = ker.orientation_2_object();
    const bool                         point_collinear =
      (orient_f (p1, p2, p3) == COLLINEAR ||
       orient_f (p1, p2, p4) == COLLINEAR ||
       orient_f (p1, p2, p5) == COLLINEAR ||
       orient_f (p1, p3, p4) == COLLINEAR ||
       orient_f (p1, p3, p5) == COLLINEAR ||
       orient_f (p1, p4, p5) == COLLINEAR ||
       orient_f (p2, p3, p4) == COLLINEAR ||
       orient_f (p2, p3, p5) == COLLINEAR ||
       orient_f (p2, p4, p5) == COLLINEAR ||
       orient_f (p3, p4, p5) == COLLINEAR);

    if (point_collinear)
    {
      _info = 0;           // Inavlid arc.
      return;
    }

    // Set the source and target.
    Rational          x1 = p1.x();
    Rational          y1 = p1.y();
    Rational          x5 = p5.x();
    Rational          y5 = p5.y();
    Nt_traits         nt_traits;

    _source = Point_2 (nt_traits.convert (x1), nt_traits.convert (y1));
    _target = Point_2 (nt_traits.convert (x5), nt_traits.convert (y5));

    // Set a conic curve that passes through the five given point.
    typename Rat_kernel::Conic_2   temp_conic;
    Rational                       rat_coeffs [6];

    temp_conic.set (p1, p2, p3, p4, p5);

    // Get the conic coefficients.
    rat_coeffs[0] = temp_conic.r();
    rat_coeffs[1] = temp_conic.s();
    rat_coeffs[2] = temp_conic.t();
    rat_coeffs[3] = temp_conic.u();
    rat_coeffs[4] = temp_conic.v();
    rat_coeffs[5] = temp_conic.w();

    // Determine the orientation: If one of the midpoints forms a left-turn
    // with the source and the target points, the orientation is positive
    // (going counterclockwise).
    // Otherwise, it is negative (going clockwise).
    const Orientation                  turn = orient_f(p1, p2, p5);

    if (turn == LEFT_TURN)
    {
      _orient = COUNTERCLOCKWISE;
      CGAL_precondition (orient_f(p1, p3, p5) == LEFT_TURN &&
                         orient_f(p1, p4, p5) == LEFT_TURN);
    }
    else
    {
      _orient = CLOCKWISE;
      CGAL_precondition (orient_f(p1, p3, p5) != LEFT_TURN &&
                         orient_f(p1, p4, p5) != LEFT_TURN);
    }

    // Set the arc properties (no need to compute the orientation).
    _set (rat_coeffs);

    // Make sure that all midpoints are strictly between the
    // source and the target.
    Point_2  mp2 = Point_2 (nt_traits.convert (p2.x()),
                            nt_traits.convert (p2.y()));
    Point_2  mp3 = Point_2 (nt_traits.convert (p3.x()),
                            nt_traits.convert (p3.y()));
    Point_2  mp4 = Point_2 (nt_traits.convert (p4.x()),
                            nt_traits.convert (p4.y()));

    if (! _is_strictly_between_endpoints (mp2) ||
        ! _is_strictly_between_endpoints (mp3) ||
        ! _is_strictly_between_endpoints (mp4))
    {
      _info = 0;               // Inalvid arc.
    }
  }

  /*!
   * Construct a conic arc which lies on the conic:
   *   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
   * The source and the target are specified by the intersection of the
   * conic with:
   *   C_1: r_1*x^2 + s_1*y^2 + t_1*xy + u_1*x + v_1*y + w_1 = 0
   *   C_2: r_2*x^2 + s_2*y^2 + t_2*xy + u_2*x + v_2*y + w_2 = 0
   * The user must also specify the source and the target with approximated
   * coordinates. The actual intersection points that best fits the source
   * (or the target) will be selected.
   */
  _Conic_arc_2 (const Rational& r, const Rational& s, const Rational& t,
                const Rational& u, const Rational& v, const Rational& w,
                const Orientation& orient,
                const Point_2& app_source,
                const Rational& r_1, const Rational& s_1, const Rational& t_1,
                const Rational& u_1, const Rational& v_1, const Rational& w_1,
                const Point_2& app_target,
                const Rational& r_2, const Rational& s_2, const Rational& t_2,
                const Rational& u_2, const Rational& v_2, const Rational& w_2):
    _orient(orient),
    _extra_data_P(nullptr)
  {
    // Create the integer coefficients of the base conic.
    Rational          rat_coeffs [6];
    Nt_traits         nt_traits;
    Integer           base_coeffs[6];
    int               deg_base;

    rat_coeffs[0] = r;
    rat_coeffs[1] = s;
    rat_coeffs[2] = t;
    rat_coeffs[3] = u;
    rat_coeffs[4] = v;
    rat_coeffs[5] = w;

    nt_traits.convert_coefficients (rat_coeffs, rat_coeffs + 6,
                                    base_coeffs);

    if (CGAL::sign (base_coeffs[0]) == ZERO &&
        CGAL::sign (base_coeffs[1]) == ZERO &&
        CGAL::sign (base_coeffs[2]) == ZERO)
    {
      deg_base = 1;
    }
    else
    {
      deg_base = 2;
    }

    // Compute the endpoints.
    Rational          aux_rat_coeffs [6];
    Integer           aux_coeffs[6];
    int               deg_aux;
    Algebraic         xs[4];
    int               n_xs;
    Algebraic         ys[4];
    int               n_ys;
    int               i, j;
    Algebraic         val;
    bool              found;
    double            dx, dy;
    double            curr_dist;
    double            min_dist = -1;
    int               k;

    for (k = 1; k <= 2; k++)
    {
      // Get the integer coefficients of the k'th auxiliary conic curve.
      aux_rat_coeffs[0] = (k == 1) ? r_1 : r_2;
      aux_rat_coeffs[1] = (k == 1) ? s_1 : s_2;
      aux_rat_coeffs[2] = (k == 1) ? t_1 : t_2;
      aux_rat_coeffs[3] = (k == 1) ? u_1 : u_2;
      aux_rat_coeffs[4] = (k == 1) ? v_1 : v_2;
      aux_rat_coeffs[5] = (k == 1) ? w_1 : w_2;

      nt_traits.convert_coefficients (aux_rat_coeffs, aux_rat_coeffs + 6,
                                      aux_coeffs);

      if (CGAL::sign (aux_coeffs[0]) == ZERO &&
          CGAL::sign (aux_coeffs[1]) == ZERO &&
          CGAL::sign (aux_coeffs[2]) == ZERO)
      {
        deg_aux = 1;
      }
      else
      {
        deg_aux = 2;
      }

      // Compute the x- and y-coordinates of intersection points of the base
      // conic and the k'th auxiliary conic.
      n_xs = _compute_resultant_roots (nt_traits,
                                       base_coeffs[0], base_coeffs[1],
                                       base_coeffs[2],
                                       base_coeffs[3], base_coeffs[4],
                                       base_coeffs[5],
                                       deg_base,
                                       aux_coeffs[0], aux_coeffs[1],
                                       aux_coeffs[2],
                                       aux_coeffs[3], aux_coeffs[4],
                                       aux_coeffs[5],
                                       deg_aux,
                                       xs);

      n_ys = _compute_resultant_roots (nt_traits,
                                       base_coeffs[1], base_coeffs[0],
                                       base_coeffs[2],
                                       base_coeffs[4], base_coeffs[3],
                                       base_coeffs[5],
                                       deg_base,
                                       aux_coeffs[1], aux_coeffs[0],
                                       aux_coeffs[2],
                                       aux_coeffs[4], aux_coeffs[3],
                                       aux_coeffs[5],
                                       deg_aux,
                                       ys);

      // Find the intersection point which is nearest the given approximation
      // and set it as the endpoint.
      found = false;
      for (i = 0; i < n_xs; i++)
      {
        for (j = 0; j < n_ys; j++)
        {
          // Check if the point (xs[i], ys[j]) lies on both conics.
          val = nt_traits.convert(base_coeffs[0]) * xs[i]*xs[i] +
                nt_traits.convert(base_coeffs[1]) * ys[j]*ys[j] +
                nt_traits.convert(base_coeffs[2]) * xs[i]*ys[j] +
                nt_traits.convert(base_coeffs[3]) * xs[i] +
                nt_traits.convert(base_coeffs[4]) * ys[j] +
                nt_traits.convert(base_coeffs[5]);

          if (CGAL::sign (val) != ZERO)
            continue;

          val = nt_traits.convert(aux_coeffs[0]) * xs[i]*xs[i] +
                nt_traits.convert(aux_coeffs[1]) * ys[j]*ys[j] +
                nt_traits.convert(aux_coeffs[2]) * xs[i]*ys[j] +
                nt_traits.convert(aux_coeffs[3]) * xs[i] +
                nt_traits.convert(aux_coeffs[4]) * ys[j] +
                nt_traits.convert(aux_coeffs[5]);

          if (CGAL::sign (val) == ZERO)
          {
            // Compute the distance of (xs[i], ys[j]) from the approximated
            // endpoint.
            if (k == 1)
            {
              dx = CGAL::to_double (xs[i] - app_source.x());
              dy = CGAL::to_double (ys[j] - app_source.y());
            }
            else
            {
              dx = CGAL::to_double (xs[i] - app_target.x());
              dy = CGAL::to_double (ys[j] - app_target.y());
            }

            curr_dist = dx*dx + dy*dy;

            // Update the endpoint if (xs[i], ys[j]) is the nearest pair so
            // far.
            if (! found || curr_dist < min_dist)
            {
              if (k == 1)
                _source = Point_2 (xs[i], ys[j]);
              else
                _target = Point_2 (xs[i], ys[j]);

              min_dist = curr_dist;
              found = true;
            }
          }
        }
      }

      if (! found)
      {
        _info = 0;           // Invalid arc.
        return;
      }
    }

    // Make sure that the source and the target are not the same.
    if (Alg_kernel().compare_xy_2_object() (_source,
                                            _target) == EQUAL)
    {
      _info = 0;      // Invalid arc.
      return;
    }

    // Set the arc properties (no need to compute the orientation).
    _set (rat_coeffs);
  }

  /*!
   * Destructor.
   */
  virtual ~_Conic_arc_2 ()
  {
    if (_extra_data_P != nullptr)
      delete _extra_data_P;
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
    if (_extra_data_P != nullptr)
      delete _extra_data_P;

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

    // Duplicate the extra data, if necessary.
    if (arc._extra_data_P != nullptr)
      _extra_data_P = new Extra_data (*(arc._extra_data_P));
    else
      _extra_data_P = nullptr;

    return (*this);
  }
  //@}

  /// \name Get the arc properties.
  //@{

  /*!
   * Check if the arc is valid.
   */
  bool is_valid () const
  {
    return ((_info & IS_VALID) != 0);
  }

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
   * Check whether the arc is x-monotone.
   */
  bool is_x_monotone () const
  {
    // Check if the arc contains no vertical tangency points.
    Point_2      vtan_ps[2];
    return (vertical_tangency_points (vtan_ps) == 0);
  }

  /*!
   * Check whether the arc is y-monotone.
   */
  bool is_y_monotone () const
  {
    // Check if the arc contains no horizontal tangency points.
    Point_2      htan_ps[2];
    return (horizontal_tangency_points (htan_ps) == 0);
  }

  /*!
   * Check whether the arc represents a full conic curve.
   */
  bool is_full_conic () const
  {
    return ((_info & IS_FULL_CONIC) != 0);
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
    CGAL_precondition (is_valid());

    double    x_min = 0, y_min = 0;
    double    x_max = 0, y_max = 0;

    if (is_full_conic())
    {
      // In case of a full conic (an ellipse or a circle), compute the
      // horizontal and vertical tangency points and use them to bound the arc.
      Point_2   tan_ps[2];
      CGAL_assertion_code(int n_tan_ps);

      CGAL_assertion_code(n_tan_ps = vertical_tangency_points(tan_ps));
      CGAL_assertion(n_tan_ps == 2);

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

      CGAL_assertion_code(n_tan_ps = horizontal_tangency_points(tan_ps));
      CGAL_assertion(n_tan_ps == 2);

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
      x_min = source_left ?
        CGAL::to_double(_source.x()) : CGAL::to_double(_target.x());
      x_max = source_left ?
        CGAL::to_double(_target.x()) : CGAL::to_double(_source.x());

      bool   source_down =
        CGAL::to_double(_source.y()) < CGAL::to_double(_target.y());
      y_min = source_down ?
        CGAL::to_double(_source.y()) : CGAL::to_double(_target.y());
      y_max = source_down ?
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

  /// \name Modifying functions.
  //@{

  /*!
   * Set the source point of the conic arc.
   * \param ps The new source point.
   * \pre The arc is not a full conic curve.
   *      ps must lie on the supporting conic curve.
   */
  void set_source (const Point_2& ps)
  {
    CGAL_precondition (! is_full_conic());
    CGAL_precondition (_is_on_supporting_conic (ps));
    CGAL_precondition (Alg_kernel().orientation_2_object()
                       (_source, ps, _target) == _orient ||
                       Alg_kernel().orientation_2_object()
                       (ps, _source, _target) == _orient);

    _source = ps;
    return;
  }

  /*!
   * Set the target point of the conic arc.
   * \param pt The new source point.
   * \pre The arc is not a full conic curve.
   *      pt must lie on the supporting conic curve.
   */
  void set_target (const Point_2& pt)
  {
    CGAL_precondition (! is_full_conic());
    CGAL_precondition (_is_on_supporting_conic (pt));
    CGAL_precondition (Alg_kernel().orientation_2_object()
                       (_source, pt, _target) == _orient ||
                       Alg_kernel().orientation_2_object()
                       (_source, _target, pt) == _orient);

    _target = pt;
    return;
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
    CGAL_assertion (m <= 2);
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
    CGAL_assertion (m <= 2);
    return (m);
  }

  /*!
   * Find all points on the arc with a given x-coordinate.
   * \param p A placeholder for the x-coordinate.
   * \param ps The point on the arc at x(p).
   * \pre The vector ps should be allocated at the size of 2.
   * \return The number of points found.
   */
  int points_at_x (const Point_2& p,
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
      ps[m] = Point_2 (p.x(), ys[i]);

      if (is_full_conic() || _is_between_endpoints(ps[m]))
        m++;
    }

    // Return the number of points on the arc.
    CGAL_assertion (m <= 2);
    return (m);
  }

  /*!
   * Find all points on the arc with a given y-coordinate.
   * \param p A placeholder for the y-coordinate.
   * \param ps The point on the arc at x(p).
   * \pre The vector ps should be allocated at the size of 2.
   * \return The number of points found.
   */
  int points_at_y (const Point_2& p,
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
      ps[m] = Point_2 (xs[i], p.y());

      if (is_full_conic() || _is_between_endpoints(ps[m]))
        m++;
    }

    // Return the number of points on the arc.
    CGAL_assertion (m <= 2);
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
    _info = IS_VALID;

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

    // Make sure both endpoint lie on the supporting conic.
    if (! _is_on_supporting_conic (_source) ||
        ! _is_on_supporting_conic (_target))
    {
      _info = 0;          // Invalid arc.
      return;
    }

    _extra_data_P = nullptr;

    // Check whether we have a degree 2 curve.
    if ((CGAL::sign (_r) != ZERO ||
         CGAL::sign (_s) != ZERO ||
         CGAL::sign (_t) != ZERO))
    {
      if (_orient == COLLINEAR)
      {
        // We have a segment of a line pair with rational coefficients.
        // Compose the equation of the underlying line
        // (with algebraic coefficients).
        const Algebraic        x1 = _source.x(), y1 = _source.y();
        const Algebraic        x2 = _target.x(), y2 = _target.y();

        // The supporting line is A*x + B*y + C = 0, where:
        //
        //  A = y2 - y1,    B = x1 - x2,    C = x2*y1 - x1*y2
        //
        // We use the extra dat field to store the equation of this line.
        _extra_data_P = new Extra_data;
        _extra_data_P->a = y2 - y1;
        _extra_data_P->b = x1 - x2;
        _extra_data_P->c = x2*y1 - x1*y2;
        _extra_data_P->side = ZERO;

        // Make sure the midpoint is on the line pair (thus making sure that
        // the two points are not taken from different lines).
        Alg_kernel       ker;
        Point_2          p_mid = ker.construct_midpoint_2_object() (_source,
                                                                    _target);

        if (CGAL::sign ((nt_traits.convert(_r)*p_mid.x() +
                         nt_traits.convert(_t)*p_mid.y() +
                         nt_traits.convert(_u)) * p_mid.x() +
                        (nt_traits.convert(_s)*p_mid.y() +
                         nt_traits.convert(_v)) * p_mid.y() +
                        nt_traits.convert(_w)) != ZERO)
        {
          _info = 0;          // Invalid arc.
          return;
        }
      }
      else
      {
        // The sign of (4rs - t^2) detetmines the conic type:
        // - if it is possitive, the conic is an ellipse,
        // - if it is negative, the conic is a hyperbola,
        // - if it is zero, the conic is a parabola.
        CGAL::Sign   sign_conic = CGAL::sign (4*_r*_s - _t*_t);

        if (sign_conic == NEGATIVE)
          // Build the extra hyperbolic data
          _build_hyperbolic_arc_data ();

        if (sign_conic != POSITIVE)
        {
          // In case of a non-degenerate parabola or a hyperbola, make sure
          // the arc is not infinite.
          Alg_kernel       ker;
          Point_2          p_mid = ker.construct_midpoint_2_object() (_source,
                                                                      _target);
          Point_2          ps[2];

          bool  finite_at_x = (points_at_x(p_mid, ps) > 0);
          bool  finite_at_y = (points_at_y(p_mid, ps) > 0);

          if (! finite_at_x && ! finite_at_y)
          {
            _info = 0;          // Invalid arc.
            return;
          }
        }
      }
    }

    // Mark that this arc valid and is not a full conic curve.
    _info = IS_VALID;

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
    const bool  is_ellipse = (CGAL::sign (4*_r*_s - _t*_t) == POSITIVE);
    CGAL_assertion (is_ellipse);

    // We do not have to store any extra data with the arc.
    _extra_data_P = nullptr;

    // Mark that this arc is a full conic curve.
    if (is_ellipse)
      _info = IS_VALID | IS_FULL_CONIC;
    else
      _info = 0;

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
    const Algebraic  _zero = 0;
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
        // phi = PI/2.
        sin_phi = _one;
        cos_phi = _zero;
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
    // We store the equation of this line in the extra data structure and also
    // the sign (side of half-plane) our arc occupies with respect to the line.
    _extra_data_P = new Extra_data;

    _extra_data_P->a = cos_phi;
    _extra_data_P->b = sin_phi;
    _extra_data_P->c = - (cos_phi*x0 + sin_phi*y0);

    // Make sure that the two endpoints are located on the same branch
    // of the hyperbola.
    _extra_data_P->side = _sign_of_extra_data (_source.x(), _source.y());

    CGAL_assertion (_extra_data_P->side != ZERO);
    CGAL_assertion (_extra_data_P->side == _sign_of_extra_data(_target.x(),
                                                              _target.y()));

    return;
  }
  //@}

protected:

  /// \name Auxiliary functions.
  //@{

  /*!
   * Evaluate the sign of (a*x + b*y + c) stored with the extra data field
   * at a given point.
   * \param px The x-coordinate of query point.
   * \param py The y-coordinate of query point.
   * \return The sign of (a*x + b*y + c).
   */
  Sign _sign_of_extra_data (const Algebraic& px,
                            const Algebraic& py) const
  {
    CGAL_assertion (_extra_data_P != nullptr);

    if (_extra_data_P == nullptr)
      return (ZERO);

    Algebraic         val = (_extra_data_P->a*px + _extra_data_P->b*py +
                             _extra_data_P->c);

    return (CGAL::sign (val));
  }

  /*!
   * Check whether the given point lies on the supporting conic of the arc.
   * \param p The query point.
   * \return (true) if p lies on the supporting conic; (false) otherwise.
   */
  bool _is_on_supporting_conic (const Point_2& p) const
  {
    // Check whether p satisfies the conic equation.
    // The point must satisfy: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0.
    Nt_traits        nt_traits;
    const Algebraic  val = (nt_traits.convert(_r)*p.x() +
                            nt_traits.convert(_t)*p.y() +
                            nt_traits.convert(_u)) * p.x() +
                           (nt_traits.convert(_s)*p.y() +
                            nt_traits.convert(_v)) * p.y() +
                           nt_traits.convert(_w);

    return (CGAL::sign (val) == ZERO);
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

    // Check if we have extra data available.
    if (_extra_data_P != nullptr)
    {
      if (_extra_data_P->side != ZERO)
      {
        // In case of a hyperbolic arc, make sure the point is located on the
        // same branch as the arc.
        if (_sign_of_extra_data(p.x(), p.y()) != _extra_data_P->side)
          return (false);
      }
      else
      {
        // In case we have a segment of a line pair, make sure that p really
        // satisfies the equation of the line.
        if (_sign_of_extra_data(p.x(), p.y()) != ZERO)
          return (false);
      }
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

      // If p is not in the (open) x-range (or y-range) of the segment, it
      // cannot be contained in the segment.
      if (res1 == EQUAL || res2 == EQUAL || res1 == res2)
        return (false);

      // Perform an orientation test: This is crucial for segment of line
      // pairs, as we want to make sure that p lies on the same line as the
      // source and the target.
      return (ker.orientation_2_object()(_source, p, _target) == COLLINEAR);
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
    n_xs = static_cast<int>(xs_end - xs);

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
      n_ys = static_cast<int>(ys_end - ys);
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
          if (CGAL::compare (nt_traits.convert(_two*_s) * ys[j],
                             -(nt_traits.convert(_t) * xs[i] +
                               nt_traits.convert(_v))) == EQUAL)
          {
            ps[n] = Point_2 (xs[i], ys[j]);
            n++;
            break;
          }
        }
      }
    }

    CGAL_assertion (n <= 2);
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
    n = static_cast<int>(ys_end - ys);

    // Compute the x coordinates and construct the horizontal tangency points.
    Algebraic     x;
    int           i;

    for (i = 0; i < n; i++)
    {
      // Having computed y, x is the single solution to the quadratic equation
      // above, and since its discriminant is 0, x is simply given by:
      x = -(nt_traits.convert(_t)*ys[i] + nt_traits.convert(_u)) /
        nt_traits.convert(_two*_r);

      ps[i] = Point_2 (x, ys[i]);
    }

    CGAL_assertion (n <= 2);
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
    // Check if we actually have a linear equation.
    if (CGAL::sign(A) == ZERO)
    {
      if (CGAL::sign(B) == ZERO)
        return (0);

      x_minus = x_plus = -C / B;
      return (1);
    }

    // Compute the discriminant and act according to its sign.
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

} //namespace CGAL

#endif
