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

#ifndef CGAL_CONIC_ARC_2_CORE_H
#define CGAL_CONIC_ARC_2_CORE_H

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Conic_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Arrangement_2/Conic_arc_2_eq.h>

#include <list>
#include <ostream>

#define CGAL_CONIC_ARC_USE_CACHING

CGAL_BEGIN_NAMESPACE

/*!
 * A class that stores additional information with the point's coordinates:
 * Namely the conic IDs of the generating curves.
 */
template <class Kernel_>
class Point_ex_2 : public Kernel_::Point_2
{
public:

  typedef Kernel_                       Kernel;
  typedef typename Kernel::Point_2      Base;
  typedef Point_ex_2<Kernel>            Self;
    
  typedef typename Kernel::FT           CoNT;
 
private:

  int    _id1;       // The ID of the first generating conic (or 0).
  int    _id2;       // The ID of the second generating conic (or 0).

 public:

  // Constructors.
  Point_ex_2 () :
    Base(),
    _id1(0),
    _id2(0)
  {}

  Point_ex_2 (const CoNT& hx, const CoNT& hy, const CoNT& hz) :
    Base(hx,hy,hz),
    _id1(0),
    _id2(0)
  {}

  Point_ex_2 (const CoNT& hx, const CoNT& hy,
	      const int& id1 = 0, const int& id2 = 0) :
    Base(hx,hy),
    _id1(id1),
    _id2(id2)
  {}

  Point_ex_2 (const Base& p) :
    Base(p),
    _id1(0),
    _id2(0)
  {}

  Point_ex_2 (const Base& p,
	      const int& id1, const int& id2 = 0) :
    Base(p),
    _id1(id1),
    _id2(id2)
  {}

  // Set the generating conic IDs.
  void set_generating_conics (const int& id1,
			      const int& id2 = 0)
  {
    _id1 = id1;
    _id2 = id2;

    return;
  }

  // Check if the given conic generates the point.
  bool is_generating_conic_id (const int& id) const
  {
    return (id != 0 &&
	    (id == _id1 || id == _id2));
  }

  // Check whether two points are equal.
  bool equals (const Self& p) const
  {
    // If the two points are the same:
    if (this == &p)
      return (true);

    // Use the parent's equality operator.
    return ((*this) == p);
  }
  
  // Compare the x coordinates.
  Comparison_result compare_x (const Self& p) const
  {
    return (CGAL_NTS compare(x(), p.x()));
  }

  // Compare the y coordinates.
  Comparison_result compare_y (const Self& p) const
  {
    return (CGAL_NTS compare(y(), p.y()));
  }

  // Compare two points lexicographically.
  Comparison_result compare_lex_xy (const Self& p) const
  {
    Comparison_result   res = this->compare_x (p);

    if (res != EQUAL)
      return (res);
    
    return (this->compare_y (p));
  }
};

/*!
 * Representation of a conic arc -- a bounded segment that lies on a conic
 * curve, the loci of all points satisfying the equation:
 *   r*x^2 + s*y^2 + t*xy + u*x + v*y +w = 0
 *
 * The class is templated with two parameters: 
 * Int_kernel_ is a kernel that provides the input points or coefficients.
 *             Int_kernel_::FT must be of an integral type (e.g. CORE:BigInt).
 * Alg_kernel_ is a geometric kernel, where Kernel_::FT is the number type
 *             for the coordinates of points, which are algebraic numbers
 *             (preferably it is CORE::Expr).
 */

static int _conics_count = 0;

template <class Int_kernel_, class Alg_kernel_> class Arr_conic_traits_2;

template <class Int_kernel_, class Alg_kernel_>
class Conic_arc_2
{
protected:
  typedef Conic_arc_2<Int_kernel_, Alg_kernel_>  Self;
        
public:

  typedef Int_kernel_                  Int_kernel;
  typedef Alg_kernel_                  Alg_kernel;

  // Define objects from the integral kernel.
  typedef typename Int_kernel::FT        CfNT;
  typedef typename Int_kernel::Point_2   Int_point_2;
  typedef typename Int_kernel::Segment_2 Int_segment_2;
  typedef typename Int_kernel::Line_2    Int_line_2;
  typedef typename Int_kernel::Circle_2  Int_circle_2;

  // Define objects from the algebraic kernel.
  typedef typename Alg_kernel::FT        CoNT;   
  typedef Point_ex_2<Alg_kernel>         Point_2;
  
 protected:

  friend class Arr_conic_traits_2<Int_kernel, Alg_kernel>;

  enum
  {
    DEGREE_0 = 0,
    DEGREE_1 = 1,
    DEGREE_2 = 2,
    DEGREE_MASK = 1 + 2,
    FULL_CONIC = 4,
    X_MONOTONE = 8,
    X_MON_UNDEFINED = 8 + 16,
    FACING_UP = 32,
    FACING_DOWN = 64,
    FACING_MASK = 32 + 64,
    IS_VERTICAL_SEGMENT = 128,
    IS_VALID = 256
  };

  CfNT     _r;              //
  CfNT     _s;              // The coefficeint of the underlying conic curve:
  CfNT     _t;              //
  CfNT     _u;              //
  CfNT     _v;              // r*x^2 + s*y^2 + t*xy + u*x + v*y +w = 0
  CfNT     _w;              //

  Orientation _orient;      // The orientation of the conic.

  int         _conic_id;    // The id of the conic.
  Point_2     _source;      // The source of the arc. 
  Point_2     _target;      // The target of the arc.
  int         _info;        // A bit array with extra information:
                            // Bit 0 & 1 - The degree of the conic 
                            //             (either 1 or 2).
                            // Bit 2     - Whether the arc is a full conic.
                            // Bit 3 & 4 - Whether the arc is x-monotone
                            //             (00, 01 or 11 - undefined).
                            // Bit 5 & 6 - Indicate whether the arc is
                            //             facing upwards or downwards (for
                            //             x-monotone curves of degree 2).
                            // Bit 7     - Is the arc a vertical segment.
                            // Bit 8     - Is the arc valid.

  // For arcs whose base is a hyperbola we store the axis (a*x + b*y + c = 0)
  // which separates the two bracnes of the hyperbola. We also store the side
  // (-1 or 1) that the arc occupies.
  struct Hyperbolic_arc_data
  {
    CoNT     a;
    CoNT     b;
    CoNT     c;
    int      side;
  };

#ifdef CGAL_CONIC_ARC_USE_CACHING
  struct Intersections
  {
    int     id1;
    int     id2;
    int     n_points;
    Point_2 ps[4];
  };
#endif

  Hyperbolic_arc_data *_hyper_P;

  // Produce a unique id for a new conic. 
  int _get_new_conic_id ()
  {
    _conics_count++;
    return (_conics_count);
  }

  /*!
   * Protected constructor: Construct an arc which is a segment of the given
   * conic arc, with a new source and target points.
   * \param The copied arc.
   * \param source The new source point.
   * \param target The new target point.
   */
  Conic_arc_2 (const Self & arc,
	       const Point_2& source, const Point_2 & target) :
    _r(arc._r), _s(arc._s), _t(arc._t), _u(arc._u), _v(arc._v), _w(arc._w),
    _orient(arc._orient),
    _conic_id(arc._conic_id),
    _source(source),
    _target(target),
    _hyper_P(NULL)
  {
    _info = (arc._info & DEGREE_MASK) | 
      X_MON_UNDEFINED;

    // Copy the hyperbolic or circular data, if necessary.
    if (arc._hyper_P != NULL)
      _hyper_P = new Hyperbolic_arc_data (*(arc._hyper_P));
    else
      _hyper_P = NULL;

    // Mark whether this is a vertical segment.
    if ((arc._info & IS_VERTICAL_SEGMENT) != 0)
    {
      _info = _info | IS_VERTICAL_SEGMENT;
    }

    // Check whether the conic is x-monotone.
    if (is_x_monotone())
    {
      _info = (_info & ~X_MON_UNDEFINED) | X_MONOTONE;

      // Check whether the facing information is set for the orginating arc.
      // If it is, just copy it - otherwise calculate it if the degree is 2.
      Comparison_result facing_res = arc.facing();

      if (facing_res != EQUAL)
	_info = _info | (facing_res == LARGER ? FACING_UP : FACING_DOWN);
      else if (_orient != CGAL::COLLINEAR)
	_set_facing();
    }
    else
    {
      _info = (_info & ~X_MON_UNDEFINED);
    }

    // Mark the arc as valid.
    _info = (_info | IS_VALID);
  }

 public:

  /*!
   * Default constructor.
   */
  Conic_arc_2 () :
    _r(0), _s(0), _t(0), _u(0), _v(0), _w(0),
    _orient(CGAL::COLLINEAR),
    _conic_id(0),
    _info(X_MON_UNDEFINED),
    _hyper_P(NULL)
  {}

  /*!
   * Copy constructor.
   * \param arc The copied arc.
   */
  Conic_arc_2 (const Self & arc) :
    _r(arc._r), _s(arc._s), _t(arc._t), _u(arc._u), _v(arc._v), _w(arc._w),
    _orient(arc._orient),
    _conic_id(arc._conic_id),
    _source(arc._source),
    _target(arc._target),
    _info(arc._info),
    _hyper_P(NULL)
  {
    // Copy the hyperbolic or circular data, if necessary.
    if (arc._hyper_P != NULL)
      _hyper_P = new Hyperbolic_arc_data (*(arc._hyper_P));
    else
      _hyper_P = NULL;
  }

  /*! 
   * Construct a conic arc which lies on the conic:
   *   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
   * \param r,s,t,u,v,w The coefficients of the supporting conic curve.
   * \param orient The orientation of the arc.
   * \param source The source point.
   * \param target The target point.
   * \pre The source and the target must be on the conic boundary and must
   * not be the same.
   */
  Conic_arc_2 (const CfNT& r, const CfNT& s, const CfNT& t,
	       const CfNT& u, const CfNT& v, const CfNT& w,
	       const Orientation& orient,
	       const Point_2& source, const Point_2& target) :
    _r(r), _s(s), _t(t), _u(u), _v(v), _w(w),
    _source(source),
    _target(target),
    _info(X_MON_UNDEFINED),
    _hyper_P(NULL)
  {
    // Make sure the conic contains the two end-points on its boundary.
    const bool    source_on_boundary = _conic_has_on_boundary(source);
    CGAL_precondition(source_on_boundary);

    const bool    target_on_boundary = _conic_has_on_boundary(target);
    CGAL_precondition(target_on_boundary);

    // Make sure that the source and the target are not the same.
    const bool    source_not_equals_target = (source != target);
    CGAL_precondition(source_not_equals_target);      

    if (! (source_on_boundary && target_on_boundary && 
	   source_not_equals_target))
    {
      // In this case, the arc is invalid.
      return;
    }

    // Set the arc properties (no need to compute the orientation).
    _orient = orient;
    _set (false);
  }

  /*! 
   * Construct a conic arc which lies on the conic curve:
   *   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0 ,
   * and whose endpoints are given by the intersections of C with a given line,
   * such that the arc lies on the positive half-plane defined by l.
   * \param r,s,t,u,v,w The coefficients of the supporting conic curve.
   * \param l The intersecting line.
   * \pre The line must have two intersection points with the conic curve.
   *      The line cannot be vertical.
   */
  Conic_arc_2 (const CfNT& r, const CfNT& s, const CfNT& t,
	       const CfNT& u, const CfNT& v, const CfNT& w,
	       const Int_line_2& l) :
    _r(r), _s(s), _t(t), _u(u), _v(v), _w(w),
    _info(X_MON_UNDEFINED),
    _hyper_P(NULL)
  {
    // Make sure that the curve is of degree 2.
    const bool      has_degree_2 = (CGAL_NTS compare(r, 0) != EQUAL || 
				    CGAL_NTS compare(s, 0) != EQUAL ||
				    CGAL_NTS compare(t, 0) != EQUAL);

    CGAL_precondition (has_degree_2);

    if (! has_degree_2)
      // In this case, the arc is invalid.
      return;

    // Find the intersection points with the line l: a*x + b*y + c = 0.
    const CfNT   a = l.a();
    const CfNT   b = l.b();
    const CfNT   c = l.c();
    Point_2      pmid;

    if (CGAL_NTS compare(b, 0) != EQUAL)
    {
      // Find the x-coordinates of the intersection points of the conic curve
      // and the line y = -(a*x + c) / b:
      CoNT      xs[2];
      int       n_xs;

      n_xs = solve_quadratic_eq<CfNT,CoNT> (b*b*r + a*(a*s - b*t),
					    2*a*c*s - b*(c*t + a*v + b*u),
					    s*c*c + b*(b*w - c*v), 
					    xs);

      CGAL_precondition(n_xs == 2);
      
      _source = Point_2 (xs[0],
			 -(a*xs[0] + c) / CoNT(b));
      _target = Point_2 (xs[1],
			 -(a*xs[1] + c) / CoNT(b));

      // Get the conic points whose x-coordinate are in the middle of the
      // two endpoints.
      CoNT      x_mid = (xs[0] + xs[1]) / 2;
      CoNT      ys[2];
      int       n_ys;

      n_ys = _conic_get_y_coordinates (x_mid, ys);

      CGAL_precondition(n_ys > 0);

      if (CGAL_NTS compare (CoNT(a)*x_mid + CoNT(b)*ys[0] + CoNT(c), 
			    0) == LARGER)
      {
	pmid = Point_2 (x_mid, ys[0]);
      }
      else
      {
	CGAL_assertion(CGAL_NTS
		       compare (CoNT(a)*x_mid + CoNT(b)*ys[1] + CoNT(c),
				0) == LARGER);
	pmid = Point_2 (x_mid, ys[1]);
      }
    }
    else
    {
      CGAL_precondition(CGAL_NTS compare(a, 0) != EQUAL);

      // Find the intersection of the vertical line x = -c / a:
      CoNT       _x = CoNT(-c) / CoNT(a);
      CoNT       ys[2];
      int        n_ys;

      n_ys = _conic_get_y_coordinates (_x, ys);
      
      CGAL_precondition(n_ys == 2);

      _source = Point_2 (_x, ys[0]);
      _target = Point_2 (_x, ys[1]);

      // Get the conic points whose y-coordinate are in the middle of the
      // two endpoints.
      CoNT      y_mid = (ys[0] + ys[1]) / 2;
      CoNT      xs[2];
      int       n_xs;

      n_xs = _conic_get_x_coordinates (y_mid, xs);

      CGAL_precondition(n_xs > 0);

      if (CGAL_NTS compare(CoNT(a)*xs[0] + CoNT(b)*y_mid + CoNT(c), 
			   0) == LARGER)
      {
	pmid = Point_2 (xs[0], y_mid);
      }
      else
      {
	CGAL_assertion(CGAL_NTS
		       compare (CoNT(a)*xs[1] + CoNT(b)*y_mid + CoNT(c),
				0) == LARGER);
	pmid = Point_2 (xs[1], y_mid);
      }
    }

    // Determine the orientation: If the mid-point forms a left-turn with
    // the source and the target points, the orientation is positive (going
    // counterclockwise).
    // Otherwise, it is negative (going clockwise).
    static Alg_kernel                  ker;
    typename Alg_kernel::Orientation_2 orient_f = ker.orientation_2_object();
    
    if (orient_f(_source, pmid, _target) == LEFT_TURN)
      _orient = CGAL::COUNTERCLOCKWISE;
    else
      _orient = CGAL::CLOCKWISE;

    // Set the arc properties (no need to compute the orientation).
    _set (false);
  }

  /*! 
   * Construct a conic arc which is the full conic:
   *   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
   * \param r,s,t,u,v,w The coefficients of the conic curve.
   * \pre The conic C must be an ellipse (so 4rs - t^2 > 0).
   */
  Conic_arc_2 (const CfNT& r, const CfNT& s, const CfNT& t,
	       const CfNT& u, const CfNT& v, const CfNT& w) :
    _r(r), _s(s), _t(t), _u(u), _v(v), _w(w),
    _hyper_P(NULL)
  {
    // Make sure the given curve is an ellipse.
    bool    is_ellipse = (CGAL_NTS compare(4*r*s - t*t, CfNT(0)) == LARGER);
    CGAL_precondition(is_ellipse);
        
    if (! is_ellipse)
      // In this case, the arc is invalid.
      return;

    // Set the arc to be the full conic (and compute the orientation).
    _set_full (true);
  }

  /*!
   * Construct a conic arc from the given line segment.
   * \param seg The line segment.
   */
  Conic_arc_2 (const Int_segment_2& seg) :
    _hyper_P(NULL)
  {
    // Get the coordinates of the endpoints.
    const Int_point_2 src = seg.source();
    const Int_point_2 trg = seg.target();
    const CfNT        x1 = src.x(), y1 = src.y();
    const CfNT        x2 = trg.x(), y2 = trg.y();

    _source = Point_2 (CoNT(x1),CoNT(y1));
    _target = Point_2 (CoNT(x2),CoNT(y2));

    // Make sure that the source and the taget are not the same.
    CGAL_precondition(_source != _target);      

    // The supporting conic is r=s=t=0, and u*x + v*y + w = 0 should hold
    // for both the source (x1,y1) and the target (x2, y2).
    const CfNT _zero = 0;
    const CfNT _one = 1;

    if (CGAL_NTS compare (x1, x2) == EQUAL)
    {
      // The supporting conic is a vertical line, of the form x = CONST.
      _r = _zero;    _s = _zero;    _t =  _zero;
      _u = _one;
      _v = _zero;
      _w = -x1;
    }
    else
    {
      // The supporting line is A*x + B*y + C = 0, where:
      //
      //  A = y2 - y1,    B = x1 - x2,    C = x2*y1 - x1*y2 
      //
      _r = _zero;    _s = _zero;    _t =  _zero;
      _u = y2 - y1;
      _v = x1 - x2;
      _w = x2*y1 - x1*y2;
    }

    // The orientation is zero in case of a linear object.
    _orient = CGAL::COLLINEAR;

    // Set the arc properties (no need to compute the orientation).
    _set (false);
  }

  /*!
   * Set a circular arc that lies on a given circle.
   * \param circ The supporting circle.
   * \param orient The orientation of the circle.
   * \param source The source point.
   * \param target The target point.
   * \pre The source and the target must be on the conic boundary and must
   * not be the same.
   */  
  Conic_arc_2 (const Int_circle_2& circ,
	       const Orientation& orient,
	       const Point_2& source, const Point_2& target) :
    _source(source),
    _target(target),
    _info(X_MON_UNDEFINED),
    _hyper_P(NULL)
  {
    // Produce the correponding conic: if the circle centre is (x0,y0)
    // and it radius is R, that its equation is:
    //   x^2 + y^2 - 2*x0*x - 2*y0*y + (x0^2 + y0^2 - R^2) = 0
    // Since this equation describes a curve with a negative (clockwise) 
    // orientation, we multiply it by -1 if necessary to obtain a positive
    // (counterclockwise) orientation.
    const Int_point_2 cc = circ.center();
    const CfNT        x0 = cc.x(), y0 = cc.y();
    const CfNT        R_sq = circ.squared_radius();
    const CfNT        _zero = 0;

    if (orient == CGAL::COUNTERCLOCKWISE)
    {
      const CfNT _minus_one = -1;
      const CfNT _two       = 2;

      _r = _minus_one;
      _s = _minus_one;
      _t = _zero;
      _u = _two*x0;
      _v = _two*y0;
      _w = R_sq - x0*x0 - y0*y0;

      _orient = CGAL::COUNTERCLOCKWISE;
    }
    else
    {
      const CfNT _one       = 1;
      const CfNT _minus_two = -2;

      _r = _one;
      _s = _one;
      _t = _zero;
      _u = _minus_two*x0;
      _v = _minus_two*y0;
      _w = x0*x0 + y0*y0 - R_sq;

      _orient = CGAL::CLOCKWISE;
    }

    // Make sure the conic contains the two end-points on its boundary.
    const bool    source_on_boundary = _conic_has_on_boundary(source);
    CGAL_precondition(source_on_boundary);

    const bool    target_on_boundary = _conic_has_on_boundary(target);
    CGAL_precondition(target_on_boundary);

    // Make sure that the source and the target are not the same.
    const bool    source_not_equals_target = (source != target);
    CGAL_precondition(source_not_equals_target);      

    if (! (source_on_boundary && target_on_boundary && 
	   source_not_equals_target))
    {
      // In this case, the arc is invalid.
      return;
    }

    // Set the arc properties (no need to compute the orientation).
    _set (false);
  }

  /*!
   * Set a circular arc that corresponds a full circle:
   * \param circ The circle.
   */ 
  Conic_arc_2 (const Int_circle_2& circ) :
    _hyper_P(NULL)
  {
    // Produce the correponding conic: if the circle centre is (x0,y0)
    // and it radius is R, that its equation is:
    //   x^2 + y^2 - 2*x0*x - 2*y0*y + (x0^2 + y0^2 - R^2) = 0
    // Note that this equation describes a curve with a negative (clockwise) 
    // orientation.
    const Int_point_2 cc = circ.center();
    const CfNT        x0 = cc.x(), y0 = cc.y();
    const CfNT        R_sq = circ.squared_radius();
    const CfNT        _zero      = 0;
    const CfNT        _one       = 1;
    const CfNT        _minus_two = -2;

    _r = _one;
    _s = _one;
    _t = _zero;
    _u = _minus_two*x0;
    _v = _minus_two*y0;
    _w = x0*x0 + y0*y0 - R_sq;

    _orient = CGAL::CLOCKWISE;

    // Set the arc to be the full conic (no need to compute the orientation).
    _set_full (false);
  }

  /*!
   * Construct a circular arc from the given three points.
   * \param p1 The source point of the circular arc.
   * \param p2 A point between p1 and p3 on the arc.
   * \param p3 The target point of the circular arc.
   * \pre The three points must not be collinear.
   */
  Conic_arc_2 (const Int_point_2& p1,
	       const Int_point_2& p2,
	       const Int_point_2& p3) :
    _info(X_MON_UNDEFINED),
    _hyper_P(NULL)
  {
    // Make sure that the source and the target are not the same.
    const bool    source_not_equals_target = (p1 != p3);
    CGAL_precondition(source_not_equals_target);      

    if (!source_not_equals_target)
    {
      // In this case, the arc is invalid.
      return;
    }

    // Get the coordinates of the points.
    const CfNT        x1 = p1.x(), y1 = p1.y();
    const CfNT        x2 = p2.x(), y2 = p2.y();
    const CfNT        x3 = p3.x(), y3 = p3.y();
    const CfNT        _zero = 0;
    const CfNT        _two  = 2;
 
    // Compute the lines: A1*x + B1*y + C1 = 0,
    //               and: A2*x + B2*y + C2 = 0,
    // where:
    const CfNT A1 = _two*(x1 - x2);
    const CfNT B1 = _two*(y1 - y2);
    const CfNT C1 = y2*y2 - y1*y1 + x2*x2 - x1*x1;

    const CfNT A2 = _two*(x2 - x3);
    const CfNT B2 = _two*(y2 - y3);
    const CfNT C2 = y3*y3 - y2*y2 + x3*x3 - x2*x2;

    // Compute the coordinates of the intersection point between the
    // two lines, given by (Nx / D, Ny / D), where:
    const CfNT Nx = B1*C2 - B2*C1;
    const CfNT Ny = A2*C1 - A1*C2;
    const CfNT D = A1*B2 - A2*B1;

    // Make sure the three points are not collinear.
    const bool points_collinear = (CGAL_NTS compare(D, _zero) == EQUAL);
    
    CGAL_precondition(!points_collinear);      

    if (points_collinear)
    {
      // In this case, the arc is invalid.
      return;
    }

    // The equation of the underlying circle is given by:
    _r = D*D;
    _s = D*D;
    _t = _zero;
    _u = -_two*D*Nx;
    _v = -_two*D*Ny;
    _w = Nx*Nx + Ny*Ny - ((D*x2 - Nx)*(D*x2 - Nx) + (D*y2 - Ny)*(D*y2 - Ny));

    // Set the endpoints.
    _source = Point_2 (CoNT(x1),CoNT(y1));
    _target = Point_2 (CoNT(x3),CoNT(y3));

    // Determine the orientation: If the mid-point forms a left-turn with
    // the source and the target points, the orientation is positive (going
    // counterclockwise).
    // Otherwise, it is negative (going clockwise).
    static Alg_kernel                  ker;
    typename Alg_kernel::Orientation_2 orient_f = ker.orientation_2_object();
    Point_2                            pmid = Point_2(CoNT(x2), CoNT(y2));
    
    if (orient_f(_source, pmid, _target) == LEFT_TURN)
      _orient = CGAL::COUNTERCLOCKWISE;
    else
      _orient = CGAL::CLOCKWISE;

    // Set the arc properties (no need to compute the orientation).
    _set (false);
  }

  /*!
   * Construct a conic arc from the given five points, specified by the
   * points p1, p2, p3, p4 and p5.
   * \param p1 The source point of the given arc.
   * \param p2,p3,p4 Points lying on the conic arc, between p1 and p5.
   * \param p5 The target point of the given arc.
   * \pre No three points are collinear.
   */
  Conic_arc_2 (const Int_point_2& p1,
	       const Int_point_2& p2,
	       const Int_point_2& p3,
	       const Int_point_2& p4,
	       const Int_point_2& p5) :
    _info(X_MON_UNDEFINED),
    _hyper_P(NULL)
  {
    // Make sure that the source and the target are not the same.
    const bool    source_not_equals_target = (p1 != p5);
    CGAL_precondition(source_not_equals_target);      

    if (!source_not_equals_target)
    {
      // In this case, the arc is invalid.
      return;
    }

    // Get the coordinates of the points.
    const CfNT        x1 = p1.x(), y1 = p1.y();
    const CfNT        x2 = p2.x(), y2 = p2.y();
    const CfNT        x3 = p3.x(), y3 = p3.y();
    const CfNT        x4 = p4.x(), y4 = p4.y();
    const CfNT        x5 = p5.x(), y5 = p5.y();

    // Set a conic curve that passes through the five given point.
    typename Int_kernel::Conic_2   temp_conic;

    temp_conic.set (p1, p2, p3, p4, p5);

    // Get the conic coefficients.
    _r = temp_conic.r();
    _s = temp_conic.s();
    _t = temp_conic.t();
    _u = temp_conic.u();
    _v = temp_conic.v();
    _w = temp_conic.w();

    // Set the source and target points.
    _source = Point_2 (CoNT(x1), CoNT(y1));
    _target = Point_2 (CoNT(x5), CoNT(y5));

    // Determine the orientation: If one of the midpoints forms a left-turn
    // with the source and the target points, the orientation is positive
    // (going counterclockwise).
    // Otherwise, it is negative (going clockwise).
    static Int_kernel                  ker;
    typename Int_kernel::Orientation_2 orient_f = ker.orientation_2_object();
    const Orientation                  turn = orient_f(p1, p2, p5);
    bool                               same_turn;

    if (turn == LEFT_TURN)
    {
      _orient = CGAL::COUNTERCLOCKWISE;
      same_turn = (orient_f(p1, p3, p5) == LEFT_TURN) &&
	          (orient_f(p1, p4, p5) == LEFT_TURN);
    }
    else
    {
      _orient = CGAL::CLOCKWISE;
      same_turn = (orient_f(p1, p3, p5) != LEFT_TURN) &&
	          (orient_f(p1, p4, p5) != LEFT_TURN);
    }

    if (! same_turn)
      // In this case, the arc is invalid.
      return;

    // Set the arc properties (no need to compute the orientation).
    _set (false);
  }

  // Construct a conic arc which lies on the conic:
  //   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
  // The source and the target are specified by the intersection of the
  // conic with:
  //   C_1: r_1*x^2 + s_1*y^2 + t_1*xy + u_1*x + v_1*y + w_1 = 0
  //   C_2: r_2*x^2 + s_2*y^2 + t_2*xy + u_2*x + v_2*y + w_2 = 0
  // The user must also specify the source and the target with approximated
  // coordinates. The actual intersection points that best fits the source 
  // (or the target) will be selected.
  Conic_arc_2 (const CfNT& r, const CfNT& s, const CfNT& t,
	       const CfNT& u, const CfNT& v, const CfNT& w,
	       const Orientation& orient,
	       const Point_2& app_source,
	       const CfNT& r_1, const CfNT& s_1, const CfNT& t_1,
	       const CfNT& u_1, const CfNT& v_1, const CfNT& w_1,
	       const Point_2& app_target,
	       const CfNT& r_2, const CfNT& s_2, const CfNT& t_2,
	       const CfNT& u_2, const CfNT& v_2, const CfNT& w_2) :
    _r(r), _s(s), _t(t), _u(u), _v(v), _w(w),
    _info(X_MON_UNDEFINED),
    _hyper_P(NULL)
  {
    // Create the two auxiliary conic arcs.
    Conic_arc_2    aux_s_arc, aux_t_arc;

    aux_s_arc._r = r_1;       aux_t_arc._r = r_2;
    aux_s_arc._s = s_1;       aux_t_arc._s = s_2;
    aux_s_arc._t = t_1;       aux_t_arc._t = t_2;
    aux_s_arc._u = u_1;       aux_t_arc._u = u_2;
    aux_s_arc._v = v_1;       aux_t_arc._v = v_2;
    aux_s_arc._w = w_1;       aux_t_arc._w = w_2;

    // Compute the source and the target.
    const CfNT _zero = 0;
    int        my_info;
    CoNT       xs[4];        // The x coordinates of intersection points.
    int        n_xs;         // Number of x coordinates.
    CoNT       ys[4];        // The y coordinates of intersection points.
    int        n_ys;         // Number of y coordinates.
    Point_2    ipts[4];      // The intersection points.
    int        n_points;     // Their number.
    CoNT       dist, best_dist;
    int        ind, best_ind;
    int        k;

    if (CGAL_NTS compare(r, _zero) == EQUAL && 
	CGAL_NTS compare(s, _zero) == EQUAL &&
	CGAL_NTS compare(t, _zero) == EQUAL)
    {
      my_info = DEGREE_1;
    }
    else
    {
      my_info = DEGREE_2;
    }

    for (k = 0; k < 2; k++)
    {
      // Get the x- and y-coordinates of the intersection points with the
      // first (if k == 0) or the second (if k == 1) auxiliary curves.
      const Conic_arc_2& arc = (k == 0) ? aux_s_arc : aux_t_arc;
      int         arc_info;

      if (CGAL_NTS compare(arc._r, _zero) == EQUAL && 
	  CGAL_NTS compare(arc._s, _zero) == EQUAL &&
	  CGAL_NTS compare(arc._t, _zero) == EQUAL)
      {
	arc_info = DEGREE_1;
      }
      else
      {
	arc_info = DEGREE_2;
      }

      n_xs = _x_coordinates_of_intersection_points (_r, _s, _t, _u, _v, _w,
						    my_info,
						    arc._r, arc._s, arc._t, 
						    arc._u, arc._v, arc._w,
						    arc_info,
						    xs);
      
      n_ys = _y_coordinates_of_intersection_points (_r, _s, _t, _u, _v, _w,
						    my_info,
						    arc._r, arc._s, arc._t, 
						    arc._u, arc._v, arc._w,
						    arc_info,
						    ys);

      // Perform the pairing process of the x and y coordinates.
      n_points = _pair_intersection_points (arc,
					    xs, n_xs,
					    ys, n_ys,
					    ipts);

      
      if (n_points == 0)
	// In this case, the arc is invalid.
	return;

      // Choose the intersection point closest to the source
      // (or the target point).
      const Point_2&  p = (k == 0) ? app_source : app_target;

      best_ind = -1;
      for (ind = 0; ind < n_points; ind++)
      {
	dist = CGAL::squared_distance (p, ipts[ind]);

	if (best_ind < 0 || CGAL_NTS compare(best_dist, dist) == LARGER)
	{
	  best_dist = dist;
	  best_ind = ind;
	}
      }

      // Set the source (or the target).
      if (k == 0)
	_source = ipts[best_ind];
      else
	_target = ipts[best_ind];
    }

    // Make sure that the source and the taget are not the same.
    if (_source == _target)
      // In this case, the arc is invalid.
      return;

    // Set the prescribed orientation.
    _orient = orient;    

    // Set the arc properties (no need to compute the orientation).
    _set (false);
  }

  /*!
   * Destructor.
   */
  virtual ~Conic_arc_2 ()
  {
    if (_hyper_P)
      delete _hyper_P;
    _hyper_P = NULL;
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
    if (_hyper_P != NULL)
      delete _hyper_P;
    _hyper_P = NULL;

    // Copy the arc's attributes.
    _r = arc._r;
    _s = arc._s;
    _t = arc._t;
    _u = arc._u;
    _v = arc._v;
    _w = arc._w;

    _orient = arc._orient;
    _conic_id = arc._conic_id;
    _source = arc._source;
    _target = arc._target;
    _info = arc._info;

    // Duplicate the data for hyperbolic arcs.
    if (arc._hyper_P != NULL)
      _hyper_P = new Hyperbolic_arc_data (*(arc._hyper_P));

    return (*this);
  }

  /*! 
   * Get the coefficients of the underlying conic.
   */
  const CfNT& r () const {return (_r);}
  const CfNT& s () const {return (_s);}
  const CfNT& t () const {return (_t);}
  const CfNT& u () const {return (_u);}
  const CfNT& v () const {return (_v);}
  const CfNT& w () const {return (_w);}

  /*!
   * Check if the arc is valid.
   * \return Whether the arc is valid.
   */
  bool is_valid () const
  {
    return ((_info & IS_VALID) != 0);
  }

  /*!
   * Get the arc's source.
   * \return The source point.
   */
  const Point_2& source () const
  {
    return (_source);
  }

  /*!
   * Get the arc's target.
   * \return The target point.
   */
  const Point_2& target () const
  {
    return (_target);
  }

  /*!
   * Check whether the two arcs are the same (have the same graph).
   * \param arc The compared arc.
   * \return (true) if the two arcs are equal.
   */
  bool equals (const Self& arc) const
  {
    // Check if (*this) and arc are really the same object:
    if (this == &arc)
      return (true);

    // Check whether all arc features are the same.
    if (_orient == arc._orient)
    {
      // Same orientation: The base conics must be the same and the sources
      // and targets must be equal.
      return (this->has_same_base_conic(arc) &&
	      _source.equals(arc._source) &&
	      _target.equals(arc._target));
    }
    else
    {
      // Opposite orientation: The base conics must be the same and the sources
      // and targets must be flipped.
      return (this->has_same_base_conic(arc) &&
	      _source.equals(arc._target) &&
	      _target.equals(arc._source));
    }
  }

  /*!
   * Check whether the two arcs have the same base conic.
   * \param arc The compared arc.
   * \return (true) if the two base conics are the same (have the same graph).
   */
  bool has_same_base_conic (const Self& arc) const
  {
    // In case the two arcs originiat from the same base conic:
    if (_conic_id == arc._conic_id)
      return (true);

    // In case both arcs are collinear, check if they have the same
    // supporting lines.
    if (_orient == CGAL::COLLINEAR && arc._orient == CGAL::COLLINEAR)
    {
      typename Alg_kernel::Line_2  l1 (_source, _target);
      typename Alg_kernel::Line_2  l2 (arc._source, arc._target);

      return (l1 == l2);
    }

    // Check whether arc equals (*this) up to a constant factor.
    const CfNT _zero = 0;
    CfNT       factor1 = 1, factor2 = 1;

    if (CGAL_NTS compare(_r, _zero) != EQUAL)
      factor1 = _r;
    else if (CGAL_NTS compare(_s, _zero) != EQUAL)
      factor1 = _s;
    else if (CGAL_NTS compare(_t, _zero) != EQUAL)
      factor1 = _t;
    else if (CGAL_NTS compare(_u, _zero) != EQUAL)
      factor1 = _u;
    else if (CGAL_NTS compare(_v, _zero) != EQUAL)
      factor1 = _v;
    else if (CGAL_NTS compare(_w, _zero) != EQUAL)
      factor1 = _w;

    if (CGAL_NTS compare(arc._r, _zero) != EQUAL)
      factor2 = arc._r;
    else if (CGAL_NTS compare(arc._s, _zero) != EQUAL)
      factor2 = arc._s;
    else if (CGAL_NTS compare(arc._t, _zero) != EQUAL)
      factor2 = arc._t;
    else if (CGAL_NTS compare(arc._u, _zero) != EQUAL)
      factor2 = arc._u;
    else if (CGAL_NTS compare(arc._v, _zero) != EQUAL)
      factor2 = arc._v;
    else if (CGAL_NTS compare(arc._w, _zero) != EQUAL)
      factor2 = arc._w;

    return (CGAL_NTS compare (_r * factor2, arc._r * factor1) == EQUAL && 
	    CGAL_NTS compare (_s * factor2, arc._s * factor1) == EQUAL &&
	    CGAL_NTS compare (_t * factor2, arc._t * factor1) == EQUAL &&
	    CGAL_NTS compare (_u * factor2, arc._u * factor1) == EQUAL &&
	    CGAL_NTS compare (_v * factor2, arc._v * factor1) == EQUAL &&
	    CGAL_NTS compare (_w * factor2, arc._w * factor1) == EQUAL);
  }

  /*!
   * Check whether the arc is a full conic (i.e. a non-degenerate ellipse).
   * \return (true) if the arc represents a full conic.
   */
  bool is_full_conic () const
  {
    return ((_info & FULL_CONIC) != 0);
  }

  /*!
   * Check whether the arc is a circular arc.
   * \return (true) if the underlying conic curve is a circle.
   */
  bool is_circular() const
  {
    return (CGAL_NTS compare(_r, _s) == EQUAL && 
	    CGAL_NTS compare(_t, 0) == EQUAL);
  }
  
  /*!
   * Check whether the arc is a line segment.
   * \return (true) if the underlying conic curve is linear.
   */
  bool is_segment () const
  {
    // Notice that the orientation is COLLINEAR if the underlying curve has
    // a degree 1, or when it is a line-pair.
    return (_orient == CGAL::COLLINEAR);
  }

  /*!
   * Check whether the arc is a vertical segment.
   * \return (true) if the arc is a vertical segment.
   */
  bool is_vertical_segment () const
  {
    return ((_info & IS_VERTICAL_SEGMENT) != 0);
  }

  /*!
   * Get a bounding box for the conic arc.
   * \return The bounding box.
   */
  Bbox_2 bbox () const
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
    Point_2 tps[2];
    int        n_tps;
    int        i;

    n_tps = vertical_tangency_points (tps);
    for (i = 0; i < n_tps; i++)
    {
      if (CGAL::to_double(tps[i].x()) < x_min)
	x_min = CGAL::to_double(tps[i].x());
      if (CGAL::to_double(tps[i].x()) > x_max)
	x_max = CGAL::to_double(tps[i].x());
    }

    // Go over the horizontal tangency points and try to update the y-points.
    n_tps = horizontal_tangency_points (tps);
    for (i = 0; i < n_tps; i++)
    {
      if (CGAL::to_double(tps[i].y()) < y_min)
	y_min = CGAL::to_double(tps[i].y());
      if (CGAL::to_double(tps[i].y()) > y_max)
	y_max = CGAL::to_double(tps[i].y());
    }
    
    // Return the resulting bounding box.
    return (Bbox_2 (x_min, y_min, x_max, y_max));
  }

  /*!
   * Check whether the given point is on the conic arc.
   * \param q The query point.
   * \return (true) if the arc contains the point q.
   */
  bool contains_point (const Point_2& q) const
  { 
    // Check whether the conic contains the point (x,y).
    if (q.is_generating_conic_id(_conic_id) ||
	_conic_has_on_boundary(q))
    {
      // If the point is on the conic boundary, it is contained in the arc
      // either if the arc is a full conic, or if it is between the two
      // endpoints of the arc.      
      return (is_full_conic() || _is_between_endpoints(q));
    }

    // If the point is not on the conic boundary, it cannot be on the arc.
    return (false);
  }

  /*!
   * Calculate the vertical tangency points of the arc.
   * \param vpts The vertical tangency points -- should be allocated at the 
   * size of 2).
   * \return The number of vertical tangency points.
   */
  int vertical_tangency_points (Point_2* vpts) const
  {
    // No vertical tangency points for segments or for x-monotone curves:
    if ((_info & DEGREE_MASK) < 2 ||
	(_info & X_MON_UNDEFINED) == X_MONOTONE)
    {
      return (0);
    }

    // Calculate the vertical tangency points of the conic.
    Point_2 ps[2];
    int     n;

    n = _conic_vertical_tangency_points (ps);

    // Return only the points that are contained in the arc interior.
    int    m = 0;
    int    i;

    for (i = 0; i < n; i++)
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
   * \param hpts The horizontal tangency points -- should be allocated at the 
   * size of 2).
   * \return The number of horizontal tangency points.
   */
  int horizontal_tangency_points (Point_2* hpts) const
  {
    // No horizontal tangency points for segments:
    if ((_info & DEGREE_MASK) < 2)
      return (0);

    // Calculate the horizontal tangency points of the conic.
    Point_2    ps[2];
    int        n;

    n = _conic_horizontal_tangency_points (ps);

    // Return only the points that are contained in the arc interior.
    int    m = 0;
    int    i;

    for (i = 0; i < n; i++)
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
   * Check whether the arc is x-monotone.
   * \return (true) if the arc is x-monotone.
   */
  bool is_x_monotone() const 
  {
    // If the answer is pre-calculated (and stored in the _info field), just
    // return it:
    int    is_x_mon = _info & X_MON_UNDEFINED;

    if (is_x_mon == 0)
      return (false);
    else if (is_x_mon == X_MONOTONE)
      return (true);

    // Check the number of vertical tangency points.
    Point_2 vpts[2];

    return (vertical_tangency_points(vpts) == 0);
  }

  /*!
   * Find all points on the arc with a given x-coordinate.
   * \param p A placeholder for the x-coordinate.
   * \param ps The point on the arc at x(p) -- should be allocated at the 
   * size of 2.
   * \return The number of points found.
   */
  int get_points_at_x (const Point_2& p,
                       Point_2 *ps) const
  {
    // Get the y coordinates of the points on the conic.
    CoNT    ys[2];
    int     n;

    n = _conic_get_y_coordinates (p.x(), ys);
    
    // Find all the points that are contained in the arc.
    int   m = 0;
    int   i;

    for (i = 0; i < n; i++)
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
   * \param ps The point on the arc at y(p) -- should be allocated at the 
   * size of 2.
   * \return The number of points found.
   */
  int get_points_at_y (const Point_2& p,
                       Point_2 *ps) const
  {
    // Get the y coordinates of the points on the conic.
    CoNT    xs[2];
    int     n;

    n = _conic_get_x_coordinates (p.y(), xs);
    
    // Find all the points that are contained in the arc.
    int   m = 0;
    int   i;

    for (i = 0; i < n; i++)
    {
      ps[m] = Point_2 (xs[i], p.y(), _conic_id);

      if (is_full_conic() || _is_between_endpoints(ps[m]))
	m++;
    }

    // Return the number of points on the arc.
    return (m);
  }
  
  /*! 
   * Flip the conic arc: change its orientation and swap it source and target.
   * \return The flipped arc.
   */
  Self flip () const
  {

    // Create a base conic with opposite orientation:
    Self     opp_arc;

    opp_arc._r = -_r;
    opp_arc._s = -_s;
    opp_arc._t = -_t;
    opp_arc._u = -_u;
    opp_arc._v = -_v;
    opp_arc._w = -_w;
    
    if (_orient == CGAL::COUNTERCLOCKWISE)
      opp_arc._orient = CGAL::CLOCKWISE;
    else if (_orient == CGAL::CLOCKWISE)
      opp_arc._orient = CGAL::COUNTERCLOCKWISE;
    else
      opp_arc._orient = _orient;         // Linear arc (a segment).

    opp_arc._conic_id = _conic_id;

    // Exchange the source and the target.
    opp_arc._source = _target;
    opp_arc._target = _source;
    opp_arc._info = _info;         // These properties do not change.

    // Copy data for hyperbolic arcs.
    if (_hyper_P != NULL)
      opp_arc._hyper_P = new Hyperbolic_arc_data (*_hyper_P);
    else
      opp_arc._hyper_P = NULL;
 
    return (opp_arc);
  }

  /*! RWRW - Allow higher order derivatives.
   * Get the i'th order derivative by x of the conic at the point p=(x,y).
   * \param p The point where we derive.
   * \param i The order of the derivatives (either 1 or 2).
   * \param slope_numer The numerator of the slope.
   * \param slope_denom The denominator of the slope.
   * \pre i should be either 1 (first order) or 2 (second order).
   */
  void derive_by_x_at (const Point_2& p, const int& i,
		       CoNT& slope_numer, CoNT& slope_denom) const
  {
    // Make sure i is either 1 or 2.
    CGAL_precondition(i == 1 || i == 2);

    // Make sure p is contained in the arc.
    CGAL_precondition(contains_point(p));

    // The derivative by x of the conic 
    //   C: {r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0}
    // at the point p=(x,y) is given by:
    //
    //           2r*x + t*y + u       alpha 
    //   y' = - ---------------- = - -------
    //           2s*y + t*x + v       beta
    //
    const CoNT _two = 2;
    const CoNT sl_numer = _two*CoNT(_r)*p.x() + CoNT(_t)*p.y() + CoNT(_u);
    const CoNT sl_denom = _two*CoNT(_s)*p.y() + CoNT(_t)*p.x() + CoNT(_v);

    if (i == 1)
    {
      if (CGAL_NTS compare (sl_denom, 0) == LARGER)
      {
	slope_numer = -sl_numer;
	slope_denom = sl_denom;
      }
      else
      {
	slope_numer = sl_numer;
	slope_denom = -sl_denom;
      }

      return;
    }

    // The second derivative is given by:
    //
    //             s*alpha^2 - t*alpha*beta + r*beta^2
    //   y'' = -2 -------------------------------------
    //                           beta^3
    //
    const CoNT sl2_numer = CoNT(_s) * sl_numer * sl_numer -
                           CoNT(_t) * sl_numer * sl_denom +
                           CoNT(_r) * sl_denom * sl_denom;
    const CoNT sl2_denom = sl_denom * sl_denom * sl_denom;

    if (CGAL_NTS compare(sl_denom, 0) == LARGER) // So sl2_denom > 0 as well.
    {
      slope_numer = -_two *sl2_numer;
      slope_denom = sl2_denom;
    }
    else
    {
      slope_numer = _two *sl2_numer;
      slope_denom = -sl2_denom;
    }

    return;
  }

  /*! RWRW - Allow higher order derivatives.
   * Get the i'th order derivative by y of the conic at the point p=(x,y).
   * \param p The point where we derive.
   * \param i The order of the derivatives (either 1 or 2).
   * \param slope_numer The numerator of the slope.
   * \param slope_denom The denominator of the slope.
   * \pre i should be either 1 (first order) or 2 (second order).
   */
  void derive_by_y_at (const Point_2& p, const int& i,
		       CoNT& slope_numer, CoNT& slope_denom) const
  {
    // Make sure i is either 1 or 2.
    CGAL_precondition(i == 1 || i == 2);

    // Make sure p is contained in the arc.
    CGAL_precondition(contains_point(p));

    // The derivative by y of the conic 
    //   C: {r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0}
    // at the point p=(x,y) is given by:
    //
    //           2s*y + t*x + v     alpha 
    //   x' = - ---------------- = -------
    //           2r*x + t*y + u      beta
    //
    const CoNT _two = 2;
    const CoNT sl_numer = _two*CoNT(_s)*p.y() + CoNT(_t)*p.x() + CoNT(_v);
    const CoNT sl_denom = _two*CoNT(_r)*p.x() + CoNT(_t)*p.y() + CoNT(_u);

    if (i == 1)
    {
      if (CGAL_NTS compare(sl_denom, 0) == LARGER)
      {
	slope_numer = -sl_numer;
	slope_denom = sl_denom;
      }
      else
      {
	slope_numer = sl_numer;
	slope_denom = -sl_denom;
      }

      return;
    }

    // The second derivative is given by:
    //
    //             r*alpha^2 - t*alpha*beta + s*beta^2
    //   x'' = -2 -------------------------------------
    //                           beta^3
    //
    const CoNT sl2_numer = CoNT(_r) * sl_numer * sl_numer -
                           CoNT(_t) * sl_numer * sl_denom +
                           CoNT(_s) * sl_denom * sl_denom;
    const CoNT sl2_denom = sl_denom * sl_denom * sl_denom;

    if (CGAL_NTS compare(sl_denom, 0) == LARGER) // So sl2_denom > 0 as well.
    {
      slope_numer = -_two *sl2_numer;
      slope_denom = sl2_denom;
    }
    else
    {
      slope_numer = _two *sl2_numer;
      slope_denom = -sl2_denom;
    }

    return;
  }

  /*!
   * Calculate the intersection points with the given arc.
   * \param arc The arc to intersect.
   * \param ps The output intersection points. 
   *           This area must be allocated at the size of 4.
   * \param inter_list_P For caching purposes.
   * \pre The two arcs do not lie on the same conic.
   * \return The number of the actual intersection points.
   */
  int intersections_with (const Self& arc,
			  Point_2* ps
#ifdef CGAL_CONIC_ARC_USE_CACHING
			  ,std::list<Intersections> *inter_list_P = NULL
#endif
			  ) const
  {
    // The two conics must not be the same.
    CGAL_precondition (! has_same_base_conic(arc));

    // First make sure that (this->degree) is >= than (arc.degree).
    if ((arc._info & DEGREE_MASK) == DEGREE_2 && 
	(_info & DEGREE_MASK) == DEGREE_1)
    {
      return (arc.intersections_with (*this, ps
#ifdef CGAL_CONIC_ARC_USE_CACHING
				      ,inter_list_P
#endif
				      ));
    }

    // Deal with vertical segments.
    if (arc.is_vertical_segment())
    {
      if (is_vertical_segment())
      {
	// Two vertical segments intersect only if they overlap.
	return (0);
      }
      
      // Find all points on our arc that have the same x coordinate as
      // the other vertical segment.
      int         n_ys;
      Point_2     xps[2];
      int         j;
      int         n = 0;

      n_ys = get_points_at_x (arc._source, xps);
      
      for (j = 0; j < n_ys; j++)
      {
	// Store this point only if it is contained on the other arc.
	if (arc.contains_point(xps[j]))
	{
	  ps[n] = Point_2 (xps[j].x(), xps[j].y(),
			   _conic_id, arc._conic_id);
	  n++;
	}
      }
      
      return (n);
    }
    else if (is_vertical_segment())
    {
      // Find all points on the other arc that have the same x coordinate as
      // our vertical segment.
      int         n_ys;
      Point_2     xps[2];
      int         j;
      int         n = 0;

      n_ys = arc.get_points_at_x (_source, xps);
      
      for (j = 0; j < n_ys; j++)
      {
	// Store this point only if it is contained on the other arc.
	if (contains_point(xps[j]))
	{
	  ps[n] = Point_2 (xps[j].x(), xps[j].y(),
			   _conic_id, arc._conic_id);
	  n++;
	}
      }
      
      return (n);
    }

    // Find all intersection points between the two base conic curves.
    Point_2   ipts[4];             // The intersection points.
    int       n_points = 0;        // Their number.
    bool      calc_points = true;

#ifdef CGAL_CONIC_ARC_USE_CACHING
    Intersections inter;
    int           k;


    if (inter_list_P != NULL &&
	(_info & DEGREE_MASK) != DEGREE_1)
    {
      int           id1 = _conic_id;
      int           id2 = arc._conic_id;
    
      inter.id1 = id1 < id2 ? id1 : id2;
      inter.id2 = id1 > id2 ? id1 : id2;
    
      typename std::list<Intersections>::iterator iter;
      for (iter = inter_list_P->begin(); iter != inter_list_P->end(); iter++)
      {
	if ((*iter).id1 == inter.id1 && (*iter).id2 == inter.id2)
	{
	  n_points = (*iter).n_points;
	  for (k = 0; k < n_points; k++)
	    ipts[k] = (*iter).ps[k];

	  calc_points = false;
	}
      }
    }
#endif // (of ifdef CGAL_CONIC_ARC_USE_CACHING)

    if (calc_points)
    {
      // Find all potential x coordinates and y coordinates of the
      // intersection points.
      const CfNT _zero = 0;
      CoNT       xs[4];        // The x coordinates of intersection points.
      int        n_xs;         // Number of x coordinates.
      CoNT       ys[4];        // The y coordinates of intersection points.
      int        n_ys;         // Number of y coordinates.

      n_xs = _x_coordinates_of_intersection_points (_r, _s, _t, _u, _v, _w,
						    _info,
						    arc._r, arc._s, arc._t, 
						    arc._u, arc._v, arc._w,
						    arc._info,
						    xs);
     
      n_ys = _y_coordinates_of_intersection_points (_r, _s, _t, _u, _v, _w,
						    _info,
						    arc._r, arc._s, arc._t, 
						    arc._u, arc._v, arc._w,
						    arc._info,
						    ys);

      // Perform the pairing process of the x and y coordinates.
      n_points = _pair_intersection_points (arc,
					    xs, n_xs,
					    ys, n_ys,
					    ipts);

#ifdef CGAL_CONIC_ARC_USE_CACHING
      if (inter_list_P != NULL &&
	  (_info & DEGREE_MASK) != DEGREE_1)
      {
	inter.n_points = n_points;	
	for (k = 0; k < n_points; k++)
	  inter.ps[k] = ipts[k];

	inter_list_P->push_front(inter);
      }
    
#endif // (of ifdef CGAL_CONIC_ARC_USE_CACHING)
    }

    // Go over all intersection points between the two base conics and return
    // only those located on both arcs.
    int      n = 0;
    int      i;

    for (i = 0; i < n_points; i++)
    {
      // Check for an exact point.
      if (contains_point(ipts[i]) &&
	  arc.contains_point(ipts[i]))
      {
	ps[n] = ipts[i];
	n++;
      }
    }

    return (n);
  }

  /*!
   * Check whether the two arcs overlap, and if so - compute the overlapping
   * portions.
   * \param arc The other conic arc.
   * \param ovlp_arc The output overlapping sub-arc.
   *                 This area should be allocated to the size of 2.
   * \return The number of overlapping sub-arcs.
   */
  int overlaps (const Self& arc,
		Self* ovlp_arcs) const
  {
    // Two arcs can overlap only if their base conics are identical.
    if (! this->has_same_base_conic (arc))
      return (0);

    // If the two arcs are completely equal, return one of them as the
    // overlapping arc.
    int       orient1 = _orient;
    int       orient2 = arc._orient;
    bool      same_or = (orient1 == orient2);
    bool      identical = false;

    if (orient1 == 0)
    {
      // That mean both arcs are really segments, so they are identical
      // if their endpoints are the same.
      if ((_source.equals(arc._source) && _target.equals(arc._target)) ||
	  (_source.equals(arc._target) && _target.equals(arc._source)))
	identical = true;
    }
    else
    {
      // If those are really curves of degree 2, than the points curves
      // are identical only if their source and target are the same and the
      // orientation is the same, or vice-versa if the orientation is opposite.
      if ((same_or && 
	   _source.equals(arc._source) && _target.equals(arc._target)) ||
	  (!same_or && 
	   _source.equals(arc._target) && _target.equals(arc._source)))
	identical = true;
    }

    if (identical)
    {
      ovlp_arcs[0] = arc;
      return (1);
    }

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
    // and target (notice that in case of segments, when the orientation is 0,
    // we make sure the two segments have the same direction).
    const Point_2 *arc_sourceP;
    const Point_2 *arc_targetP;

    if (orient1 == 0)
      orient1 = (_source.compare_lex_xy(_target) 
		 == LARGER) ? 1 : -1;
    if (orient2 == 0)
      orient2 = (arc._source.compare_lex_xy(arc._target) 
		 == LARGER) ? 1 : -1;

    // Check the overlap cases:
    if (orient1 == orient2)
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
	  ovlp_arcs[0] = Self(*this,_source, *arc_targetP);
	  ovlp_arcs[1] = Self(*this, *arc_sourceP, _target);
	  return (2);
	}

	// Case 1 - *this:     +----------->     
        //            arc:       +=====>
	ovlp_arcs[0] = Self(*this, *arc_sourceP,*arc_targetP);
	return (1);
      }
      else
      {
	// Case 2 - *this:     +----------->     
        //            arc:               +=====>
	ovlp_arcs[0] = Self(*this, *arc_sourceP, _target);
	return (1);
      }
    }
    else if (_is_strictly_between_endpoints(*arc_targetP))
    {
      // Case 3 - *this:     +----------->     
      //            arc:   +=====>
      ovlp_arcs[0] = Self(*this, _source, *arc_targetP);
      return (1);
    }
    else if (arc._is_between_endpoints(_source) &&
             arc._is_between_endpoints(_target) &&
	     (arc._is_strictly_between_endpoints(_source) ||
              arc._is_strictly_between_endpoints(_target)))
    {
      // Case 4 - *this:     +----------->     
      //            arc:   +================>
      ovlp_arcs[0] = *this;
      return (1);
    }
    
    // If we reached here, there are no overlaps:
    return (0);
  }
	
  /*!
   * Check whether the arc is facing up or facing down.
   * \return LARGER if the arcs is facing up, or SMALLER if it is facing down.
   *         If the arc is a line segment, EQUAL is returned.
   */
  Comparison_result facing () const
  {
    if ((_info & FACING_MASK) == 0)
      return (EQUAL);
    else if ((_info & FACING_UP) != 0)
      return (LARGER);
    else
      return (SMALLER);
  }

 protected:

  /*!
   * Set the properties of a conic arc (for the usage of the constructors).
   * \param comp_orient Should we compute the orientation of the given curve.
   */
  void _set (const bool& comp_orient)
  {
    // Initialize the information bits.
    _info = X_MON_UNDEFINED;

    // Set the orientation of conic arc.
    typename Cartesian<CfNT>::Conic_2   temp_conic (_r, _s, _t, _u, _v, _w);

    if (comp_orient)
    {
      // Compute the orientation.
      _orient = temp_conic.orientation();
    }
    else if (_orient != temp_conic.orientation())
    {
      // If the computed orientation does not match the current value,
      // multiply all conic coefficients by -1 (negate the curve).
      _r = - _r;
      _s = - _s;
      _t = - _t;
      _u = - _u;
      _v = - _v;
      _w = - _w;
    }

    // Find the degree and make sure the conic is not invalid.
    const CfNT _zero = 0;
    int        deg;
 
    if (CGAL_NTS compare(_r, _zero) != EQUAL || 
	CGAL_NTS compare(_s, _zero) != EQUAL ||
	CGAL_NTS compare(_t, _zero) != EQUAL)
    {
      // In case one of the coefficients of x^2,y^2 or xy is not zero, the
      // degree is 2.
      deg = 2;
    }
    else if (CGAL_NTS compare(_u, _zero) != EQUAL || 
	     CGAL_NTS compare(_v, _zero) != EQUAL)
    {
      // In case of a line - the degree is 1.
      deg = 1;
      _orient = CGAL::COLLINEAR;
    }
    else
    {
      // Empty conic!
      deg = 0;
    }

    CGAL_precondition(deg > 0);

    // Store the degree information.
    _info = _info | deg;

    // In case the base conic is a hyperbola, build the hyperbolic data
    // (this happens when (4rs - t^2) < 0).
    const Comparison_result det = (CGAL_NTS compare(4*_r*_s, _t*_t));
    
    if (deg == 2 && det == SMALLER && _orient != CGAL::COLLINEAR)
      _build_hyperbolic_arc_data ();
    else
      _hyper_P = NULL;

    // In case of a non-degenerate parabola or a hyperbola, make sure 
    // the arc is not infinite.
    if (deg == 2 && det != LARGER && _orient != CGAL::COLLINEAR)
    {
      CGAL_precondition_code(
      const CoNT       _two = 2;
      const Point_2    p_mid ((_source.x() + _target.x()) / _two,
                              (_source.y() + _target.y()) / _two);
      Point_2          ps[2];

      bool  finite_at_x = (get_points_at_x(p_mid, ps) > 0);
      bool  finite_at_y = (get_points_at_y(p_mid, ps) > 0);
      );
      CGAL_precondition(finite_at_x && finite_at_y);
    }

    // If we reached here, the conic arc is legal: Get a new id for the conic.
    _conic_id = _get_new_conic_id();

    _source.set_generating_conics (_conic_id);
    _target.set_generating_conics (_conic_id);

    // Check whether the conic is x-monotone.
    if (is_x_monotone())
    {
      _info = (_info & ~X_MON_UNDEFINED) | X_MONOTONE;

      // In case the conic is od degree 2, determine where is it facing.
      if (_orient != CGAL::COLLINEAR)
	_set_facing();
    }
    else
    {
      _info = (_info & ~X_MON_UNDEFINED);
    }

    // Check if the arc is a vertical segment.
    if (_orient == CGAL::COLLINEAR)
    {
      // A vertical segment is contained in the degenerate conic: u*x + w = 0,
      // but it can also be a part of a line-pair -- in this case we just check
      // that the x-coordinates of the source and the target are the same.
      const CfNT _zero = 0;

      if (((_info & DEGREE_MASK) == 1 && 
	   CGAL_NTS compare(_v, _zero) == EQUAL) ||
	   _source.compare_x (_target) == EQUAL)
      {
	_info = _info | IS_VERTICAL_SEGMENT;
      }
    }

    // Mark that the arc is valid.
    _info = (_info | IS_VALID);

    return;
  }

  /*!
   * Set the properties of a conic arc that is really a full curve
   * (that is, an ellipse).
   * \param comp_orient Should we compute the orientation of the given curve.
   */
  void _set_full (const bool& comp_orient)
  {
    // Initialize the information bits.
    _info = 0;

    // Set the orientation of conic arc.
    typename Cartesian<CfNT>::Conic_2   temp_conic (_r, _s, _t, _u, _v, _w);

    if (comp_orient)
    {
      // Compute the orientation.
      _orient = temp_conic.orientation();
    }
    else if (_orient != temp_conic.orientation())
    {
      // If the computed orientation does not match the current value,
      // multiply all conic coefficients by -1 (negate the curve).
      _r = - _r;
      _s = - _s;
      _t = - _t;
      _u = - _u;
      _v = - _v;
      _w = - _w;
    }

    // Set the information: a full conic, which is obvoiusly not x-monotone.
    _info = DEGREE_2 | FULL_CONIC;
    _hyper_P = NULL;

    // Assign one of the vertical tangency points as both the source and
    // the target of the conic arc.
    Point_2    vpts[2];
    int        n_vpts;

    n_vpts = _conic_vertical_tangency_points (vpts);

    CGAL_assertion(n_vpts > 0);
    CGAL_assertion(_conic_has_on_boundary(vpts[0]));

    // If we reached here, the conic arc is legal: Get a new id for the conic.
    _conic_id = _get_new_conic_id();

    _source = _target = vpts[0];
    _source.set_generating_conics (_conic_id);
    _target.set_generating_conics (_conic_id);

    // Mark that the arc is valid.
    _info = (_info | IS_VALID);

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
    const int   or_fact = (_orient == CGAL::CLOCKWISE) ? -1 : 1;
    const CoNT  r = or_fact * CoNT(_r);
    const CoNT  s = or_fact * CoNT(_s);
    const CoNT  t = or_fact * CoNT(_t);
    const CoNT  cos_2phi = (r - s) / CGAL::sqrt((r-s)*(r-s) + t*t);
    const CoNT  _zero = 0;
    const CoNT  _one = 1;
    const CoNT  _two = 2;
    CoNT        sin_phi;
    CoNT        cos_phi;

    // Calculate sin(phi) and cos(phi) according to the half-angle formulae:
    // 
    //  sin(phi)^2 = 0.5 * (1 - cos(2*phi))
    //  cos(phi)^2 = 0.5 * (1 + cos(2*phi))
    if (CGAL_NTS compare (t, _zero) == EQUAL)
    {
      // sin(2*phi) == 0, so phi = 0 or phi = PI/2
      if (CGAL_NTS compare (cos_2phi, _zero) == LARGER)
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
    else if (CGAL_NTS compare (t, _zero) == LARGER)
    {
      // sin(2*phi) > 0 so 0 < phi < PI/2.
      sin_phi = CGAL::sqrt((_one + cos_2phi) / _two);
      cos_phi = CGAL::sqrt((_one - cos_2phi) / _two);
    }
    else
    {
      // sin(2*phi) < 0 so PI/2 < phi < PI.
      sin_phi = CGAL::sqrt((_one + cos_2phi) / _two);
      cos_phi = -CGAL::sqrt((_one - cos_2phi) / _two);
    }
    
    // Calculate the center (x0, y0) of the conic, given by the formulae:
    //
    //        t*v - 2*s*u                t*u - 2*r*v
    //  x0 = -------------   ,     y0 = -------------
    //        4*r*s - t^2                4*r*s - t^2
    //
    // The denominator (4*r*s - t^2) must be negative for hyperbolas.
    const CoNT  u = or_fact * CoNT(_u);
    const CoNT  v = or_fact * CoNT(_v);
    const CoNT  det = 4*r*s - t*t;
    CoNT        x0, y0;

    CGAL_assertion (CGAL_NTS compare(det, _zero) == SMALLER);
    
    x0 = (t*v - _two*s*u) / det;
    y0 = (t*u - _two*r*v) / det;
    
    // The axis separating the two branches of the hyperbola is now given by:
    // 
    //  cos(phi)*x + sin(phi)*y - (cos(phi)*x0 + sin(phi)*y0) = 0
    //
    _hyper_P = new Hyperbolic_arc_data;

    _hyper_P->a = cos_phi;
    _hyper_P->b = sin_phi;
    _hyper_P->c = - (cos_phi*x0 + sin_phi*y0);

    // Make sure that the two endpoints are located on the same branch
    // of the hyperbola.
    _hyper_P->side = _hyperbolic_arc_side(_source);

    CGAL_assertion (_hyper_P->side = _hyperbolic_arc_side(_target));

    return;
  }

  /*!
   * Find on which branch of the hyperbola is the given point located.
   * The point is assumed to be on the hyperbola.
   * \param p The query point.
   * \return The branch ID (either -1 or 1).
   */
  int _hyperbolic_arc_side (const Point_2& p) const
  {
    if (_hyper_P == NULL)
      return (0);

    CoNT       val;

    val = _hyper_P->a*p.x() + _hyper_P->b*p.y() + _hyper_P->c;
    return ((CGAL_NTS compare(val, 0) == LARGER) ? 1 : -1);
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
    if (p.equals(_source) || p.equals(_target))
      return (true);
    else
      return (_is_strictly_between_endpoints(p));
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
    if (_hyper_P != NULL)
    {
      if (_hyperbolic_arc_side(p) != _hyper_P->side)
	return (false);
    }

    // In case of a line-pair (curve of degree 2 with a collinear orientation)
    // first make sure that p is collinear with the source and target.
    if ((_info & DEGREE_MASK) == DEGREE_2 && _orient == CGAL::COLLINEAR)
    {
      static Alg_kernel                  ker;
      typename Alg_kernel::Orientation_2 orient_f = ker.orientation_2_object();

      if (orient_f(_source, p, _target) != COLLINEAR)
	return (false);
    }

    // Act according to the degree and orientation of the conic.
    if ((_info & DEGREE_MASK) == DEGREE_1 || _orient == CGAL::COLLINEAR)
    {
      // The underlying curve is a straight line (or a line-pair).
      if (is_vertical_segment())
      {
	// In case of a vertical segment - just check whether the y coordinate
	// of p is between those of the source's and of the target's.
	Comparison_result r1 = compare_y (p, _source);
	Comparison_result r2 = compare_y (p, _target);

	return ((r1 == SMALLER && r2 == LARGER) ||
		(r1 == LARGER && r2 == SMALLER));
      }
      else
      {
	// Otherwise, since the segment is x-monotone, just check whether the
	// x coordinate of p is between those of the source's and of the 
	// target's.
	Comparison_result r1 = compare_x (p, _source);
	Comparison_result r2 = compare_x (p, _target);

	return ((r1 == SMALLER && r2 == LARGER) ||
		(r1 == LARGER && r2 == SMALLER));
      }
    }
    else
    {
      // In case of a conic of degree 2, make a decision based on the conic's
      // orientation and whether (source,p,target) is a right or a left turn.
      static Alg_kernel                  ker;
      typename Alg_kernel::Orientation_2 orient_f = ker.orientation_2_object();

      if (_orient == CGAL::COUNTERCLOCKWISE)
	return (orient_f(_source, p, _target) == LEFT_TURN);
      else
	return (orient_f(_source, p, _target) == RIGHT_TURN);
    }
  }

  /*!
   * Check whether the underlying conic contains a point on its boundary.
   * \param q The query point.
   * \return (true) if the underlying conic contains the point on its boundary.
   */
  bool _conic_has_on_boundary (const Point_2& q) const
  {
    return (_conic_has_on_boundary (q.x(), q.y()));
  }

  /*!
   * Check whether the underlying conic contains (x,y) on its boundary.
   * \param x The x coordinate of the query point.
   * \param y The y coordinate of the query point.
   * \return (true) if the underlying conic contains the point on its boundary.
   */
  bool _conic_has_on_boundary (const CoNT& x, const CoNT& y) const
  {
    const CoNT _zero = 0;
    CoNT       val;

    // The point must satisfy: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0.
    val = (CoNT(_r)*x + CoNT(_t)*y + CoNT(_u))*x +
          (CoNT(_s)*y + CoNT(_v))*y + CoNT(_w);

    return (CGAL_NTS compare (val, _zero) == EQUAL);
  }

  /*!
   * Find the y coordinates of the underlying conic at a given x coordinate.
   * \param x The x coordinate.
   * \param ys The output y coordinates. 
   *           This area must be allocated at the size of 2.
   * \return The number of y coordinates computed (either 0, 1 or 2).
   */
  int _conic_get_y_coordinates (const CoNT& x,
                                CoNT *ys) const
  {
    // Solve the quadratic equation for a given x and find the y values:
    //  s*y^2 + (t*x + v)*y + (r*x^2 + u*x + w) = 0
    return (solve_quadratic_eq<CoNT,CoNT> 
	    (CoNT(_s),
	     CoNT(_t)*x + CoNT(_v),
	     (CoNT(_r)*x + CoNT(_u))*x + CoNT(_w),
	     ys));
  }

  /*!
   * Find the x coordinates of the underlying conic at a given y coordinate.
   * \param y The y coordinate.
   * \param xs The output x coordinates. 
   *           This area must be allocated at the size of 2.
   * \return The number of x coordinates computed (either 0, 1 or 2).
   */
  int _conic_get_x_coordinates (const CoNT& y,
                                CoNT *xs) const
  {
    // Solve the quadratic equation for a given y and find the x values:
    //  r*x^2 + (t*y + u)*x + (s*y^2 + v*y + w) = 0
    return (solve_quadratic_eq<CoNT,CoNT> 
	    (CoNT(_r),
	     CoNT(_t)*y + CoNT(_u),
	     (CoNT(_s)*y + CoNT(_v))*y + CoNT(_w),
	     xs));
  }
  
  /*!
   * Find the vertical tangency points of the undelying conic.
   * \param ps The output points of vertical tangency.
   *           This area must be allocated at the size of 2.
   * \return The number of vertical tangency points.
   */
  int _conic_vertical_tangency_points (Point_2* ps) const
  {
    const CfNT _zero = 0;

    // In case the base conic is of degree 1 (and not 2), the arc has no
    // vertical tangency points.
    if ((_info & DEGREE_MASK) == DEGREE_1 || 
	CGAL_NTS compare(_s,_zero) == EQUAL)
    {
      return (0);
    }

    // We are interested in the x coordinates where the quadratic equation:
    //  s*y^2 + (t*x + v)*y + (r*x^2 + u*x + w) = 0
    // has a single solution (obviously if s = 0, there are no such points).
    // We therefore demand that the discriminant of this equation is zero:
    //  (t*x + v)^2 - 4*s*(r*x^2 + u*x + w) = 0
    const CfNT _two = 2;
    const CfNT _four = 4;
    CoNT       xs[2];
    int        n_xs;

    n_xs = solve_quadratic_eq<CfNT,CoNT> (_t*_t - _four*_r*_s,
					  _two*_t*_v - _four*_s*_u,
					  _v*_v - _four*_s*_w,
					  xs);

    // Find the y-coordinates of the vertical tangency points.
    CoNT     ys[2];
    int      n_ys;

    if (CGAL_NTS compare(_t, _zero) == EQUAL)
    {
      // The two vertical tangency points have the same y coordinate:
      ys[0] = CoNT(-_v) / CoNT(_two*_s);
      n_ys = 1;
    }
    else
    {
      n_ys = solve_quadratic_eq<CfNT,CoNT> (_four*_r*_s*_s - _s*_t*_t,
					    _four*_r*_s*_v - _two*_s*_t*_u,
					    _r*_v*_v - _t*_u*_v + _t*_t*_w,
					    ys);
    }

    // Pair the x and y coordinates and obtain the vertical tangency points.
    int   n = 0;
    int   i, j;

    for (i = 0; i < n_xs; i++)
    {
      if (n_ys == 1)
      {
        ps[n] = Point_2 (xs[i], ys[0],
			 _conic_id);
	n++;
      }
      else
      {
	for (j = 0; j < n_ys; j++)
	{
	  if (CGAL_NTS compare (ys[j],
				-(CoNT(_t)*xs[i] + CoNT(_v)) / 
				CoNT(_two*_s)) == EQUAL)
	  {
	    ps[n] = Point_2 (xs[i], ys[j],
			     _conic_id);
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
    const CfNT _zero = 0;

    // In case the base conic is of degree 1 (and not 2), the arc has no
    // vertical tangency points.
    if ((_info & DEGREE_MASK) == DEGREE_1 || 
	CGAL_NTS compare(_r, _zero) == EQUAL)
    {
      return (0);
    }

    // We are interested in the y coordinates were the quadratic equation:
    //  r*x^2 + (t*y + u)*x + (s*y^2 + v*y + w) = 0
    // has a single solution (obviously if r = 0, there are no such points).
    // We therefore demand that the discriminant of this equation is zero:
    //  (t*y + u)^2 - 4*r*(s*y^2 + v*y + w) = 0
    const CfNT _two = 2;
    const CfNT _four = 4;
    int        n;
    CoNT       ys[2];

    n = solve_quadratic_eq<CfNT,CoNT> (_t*_t - _four*_r*_s,
				       _two*_t*_u - _four*_r*_v,
				       _u*_u - _four*_r*_w,
				       ys);

    // Compute the x coordinates and construct the horizontal tangency points.
    CoNT       x;
    int        i;

    for (i = 0; i < n; i++)
    {
      // Having computed y, x is the simgle solution to the quadratic equation
      // above, and since its discriminant is 0, x is simply given by:
      x = -(CoNT(_t)*ys[i] + CoNT(_u)) / CoNT(_two*_r);

      ps[i] = Point_2 (x, ys[i],
		       _conic_id);
    }
      
    return (n);
  }
  
  /*!
   * Set the facing information for the (x-monotone) arc: It is facing up if
   * it lies above the line segments that connect its two endpoints, and facing
   * down if it lies below it.
   */
  void _set_facing ()
  {
    // Check whether the arc (which is x-monotone of degree 2) lies above or 
    // below the segement that contects its two end-points (x1,y1) and (x2,y2).
    // To do that, we find the y coordinate of a point on the arc whose x
    // coordinate is (x1+x2)/2 and compare it to (y1+y2)/2.
    const CoNT   _two = 2;
    const CoNT   x_mid = (_source.x() + _target.x()) / _two;
    const CoNT   y_mid = (_source.y() + _target.y()) / _two;
    Point_2      p_mid (x_mid, y_mid);
    Point_2      ps[2];
    int          n_ps;

    n_ps = get_points_at_x (p_mid, ps);

    CGAL_assertion (n_ps == 1);

    Comparison_result res = ps[0].compare_y (p_mid);

    if (res == LARGER)
    {
      // The arc is above the connecting segment, so it is facing upwards.
      _info = _info | FACING_UP;
    }
    else if (res == SMALLER)
    {
      // The arc is below the connecting segment, so it is facing downwards.
      _info = _info | FACING_DOWN;
    }
    
    CGAL_assertion(res != EQUAL);
    return;
  }

  /*!
   * Compute the x-coordinates of the intersection points of the two conics:
   *   C1: r1*x^2 + s1*y^2 + t1*xy + u1*x + v1*y + w1 = 0,
   * and:
   *   C2: r2*x^2 + s2*y^2 + t2*xy + u2*x + v2*y + w2 = 0.
   * \param xs The output x-coordinates (must be allocated to the size of 4).
   * \return The number of distinct x-coordinates.
   */
  int _x_coordinates_of_intersection_points (const CfNT& r1, const CfNT& s1,
					     const CfNT& t1, const CfNT& u1,
					     const CfNT& v1, const CfNT& w1,
					     const int& info1,
					     const CfNT& r2, const CfNT& s2,
					     const CfNT& t2, const CfNT& u2,
					     const CfNT& v2, const CfNT& w2,
					     const int& info2,
					     CoNT* xs) const
  {
    if ((info1 & DEGREE_MASK) == DEGREE_2 &&
	(info2 & DEGREE_MASK) == DEGREE_1)
    {
      // If necessary, swap roles between the two curves, so that the first
      // curve always has the minimal degree.
      return (_x_coordinates_of_intersection_points (r2, s2, t2, u2, v2, w2, 
						     info2,
						     r1, s1, t1, u1, v1, w1, 
						     info1,
						     xs));
    }

    const CfNT _zero = 0;
    const CfNT _two = 2;
    CfNT       c[5];

    // Check the case that the first conic has no quadratic coefficients.
    if ((info1 & DEGREE_MASK) == DEGREE_1)
    {
      if (CGAL_NTS compare(v1, _zero) == EQUAL)
      {
	// The first line is u1*x + w1 = 0, therefore:
	xs[0] = CoNT(-w1) / CoNT(u1);
	return (1);
      }
      
      // We can write the first curve as: y = (u1*x + w1) / v1.
      if ((info2 & DEGREE_MASK) == DEGREE_1)
      {
	// The second curve is also a line. We therefore get the linear
	// equation c[1]*x + c[0] = 0:
	c[1] = v1*u2 - u1*v2;
	c[0] = v1*w2 - w1*v2;

	if (CGAL_NTS compare(c[1], _zero) == EQUAL)
	  return (0);

	xs[0] = CoNT(-c[0]) / CoNT(c[1]);
	return (1);
      }

      // We substitute this expression into the equation of the second
      // conic, and get the quadratic equation c[2]*x^2 + c[1]*x + c[0] = 0:
      c[2] = u1*u1*s2 - u1*v1*t2 + v1*v1*r2;
      c[1] = _two*u1*w1*s2 - u1*v1*v2 - v1*w1*t2 + v1*v1*u2;
      c[0] = w1*w1*s2 - v1*w1*v2 + v1*v1*w2;

      return (solve_quadratic_eq<CfNT, CoNT> (c[2], c[1], c[0], xs)); 
    }

    // At this stage, both curves have degree 2. We obtain a qaurtic polynomial
    // whose roots are the x-coordinates of the intersection points.
    if (CGAL_NTS compare(s1, _zero) == EQUAL && 
	CGAL_NTS compare(s2, _zero) == EQUAL)
    {
      // If both s1 and s2 are zero, we can write the two curves as:
      //   C1: (t1*x + v1)*y + (r1*x^2 + u1*x + w1) = 0
      //   C2: (t2*x + v2)*y + (r2*x^2 + u2*x + w2) = 0
      // By writing the resultant of these two polynomials we get:
      c[4] = _zero;
      c[3] = r2*t1 - r1*t2;
      c[2] = t1*u2 - t2*u1 + r2*v1 - r1*v2;
      c[1] = t1*w2 - t2*w1 + u2*v1 - u1*v2;
      c[0] = v1*w2 - v2*w1;
    }
    else
    {
      // We can write the two curves as:
      //   C1: (s1)*y^2 + (t1*x + v1)*y + (r1*x^2 + u1*x + w1) = 0
      //   C2: (s2)*y^2 + (t2*x + v2)*y + (r2*x^2 + u2*x + w2) = 0
      // By writing the resultant of these two polynomials we get:
      c[4] = -_two*s1*s2*r1*r2 + s1*t2*t2*r1 - s1*t2*t1*r2 +
	s1*s1*r2*r2 - s2*t1*r1*t2 + s2*t1*t1*r2 + s2*s2*r1*r1;

      c[3] = -t2*r1*v1*s2 - u2*t1*t2*s1 - v2*r1*t1*s2 -
	r2*t1*v2*s1 - _two*s1*s2*r1*u2 - t2*u1*t1*s2 + u2*t1*t1*s2 -
	r2*v1*t2*s1 + u1*t2*t2*s1 + _two*v2*r1*t2*s1 + _two*u2*r2*s1*s1 + 
	_two*r2*v1*t1*s2 + _two*u1*r1*s2*s2 - _two*s1*s2*u1*r2;

      c[2] = -r2*v1*v2*s1 + u2*u2*s1*s1 + _two*w2*r2*s1*s1 +
	_two*u2*v1*t1*s2 - u2*v1*t2*s1 + w2*t1*t1*s2 - _two*s1*s2*u1*u2 - 
	w2*t1*t2*s1 + v2*v2*r1*s1 + u1*u1*s2*s2 - v2*r1*v1*s2 +
	_two*w1*r1*s2*s2 - u2*t1*v2*s1 - t2*u1*v1*s2 - _two*s1*s2*r1*w2 -
	_two*s1*s2*w1*r2 + r2*v1*v1*s2 + w1*t2*t2*s1 - v2*u1*t1*s2 -
	t2*w1*t1*s2 + _two*v2*u1*t2*s1;

      c[1] = _two*w2*u2*s1*s1 + _two*w2*v1*t1*s2 - w2*v1*t2*s1 +
	_two*v2*w1*t2*s1 + _two*w1*u1*s2*s2 - v2*u1*v1*s2 - _two*s1*s2*u1*w2 -
	v2*w1*t1*s2 + u2*v1*v1*s2 - t2*w1*v1*s2 - w2*t1*v2*s1 + 
	v2*v2*u1*s1 - u2*v1*v2*s1 - _two*s1*s2*w1*u2;

      c[0] = s2*v1*v1*w2 - s1*v2*v1*w2 - s2*v1*w1*v2 + s2*s2*w1*w1 -
	_two*s1*s2*w1*w2 + s1*w1*v2*v2 + s1*s1*w2*w2;
    }

    // Now solve the quartic equation:
    //   c[4]*x^4 + c[3]*x^3 + c[2]*x^2 + c[1]*x + c[0] = 0.
    return (solve_quartic_eq<CfNT,CoNT> (c[4], c[3], c[2], c[1], c[0], xs));
  }

  /*!
   * Compute the y-coordinates of the intersection points of the two conics:
   *   C1: r1*x^2 + s1*y^2 + t1*xy + u1*x + v1*y + w1 = 0,
   * and:
   *   C2: r2*x^2 + s2*y^2 + t2*xy + u2*x + v2*y + w2 = 0.
   * \param ys The output y-coordinates (must be allocated to the size of 4).
   * \return The number of distinct y-coordinates.
   */
  int _y_coordinates_of_intersection_points (const CfNT& r1, const CfNT& s1,
					     const CfNT& t1, const CfNT& u1,
					     const CfNT& v1, const CfNT& w1,
					     const int& info1,
					     const CfNT& r2, const CfNT& s2,
					     const CfNT& t2, const CfNT& u2,
					     const CfNT& v2, const CfNT& w2,
					     const int& info2,
					     CoNT* ys) const
  {
    // Swap roles between x and y and compute the x-coordinates.
    return (_x_coordinates_of_intersection_points (s1, r1, t1, v1, u1, w1,
						   info1,
						   s2, r2, t2, v2, u2, w2,
						   info2,
						   ys));
  }

  /*! 
   * Pair the x coordinates and the y coordinates of the intersection point
   * of the underlying conics of (*this) and arc, and return a vector of 
   * intersection points.
   * \param arc The other arc.
   * \param xs The potential x coordinates.
   * \param n_xs Number of x coordinates.
   * \param ys The potential y coordinates.
   * \param n_ys Number of y coordinates.
   * \param ipts The points that lie on both conics.
   *             This area must be allocated to the size of 4.
   * \return The number of intersection points between the conics.
   */
  int _pair_intersection_points (const Self& arc,
				 const CoNT* xs, const int& n_xs,
				 const CoNT* ys, const int& n_ys,
				 Point_2* ipts) const
  {
    int        n_ipts = 0;
    int        i;
    int        j;

    for (i = 0; i < n_xs; i++)
    {
      for (j = 0; j < n_ys; j++)
      {
	// If the current pair of x and y coordinates lies on both underlying
	// conic curves, accept it.
	if (_conic_has_on_boundary (xs[i], ys[j]) &&
	    arc._conic_has_on_boundary (xs[i], ys[j]))
	{
	  CGAL_assertion(n_ipts < 4);

	  ipts[n_ipts] = Point_2 (xs[i], ys[j],
				  _conic_id, arc._conic_id);
     
	  n_ipts++;
	}
      }
    }
    
    return (n_ipts);
  }

};

#ifndef NO_OSTREAM_INSERT_CONIC_ARC_2
template <class CfNT, class Kernel>
std::ostream& operator<< (std::ostream& os, 
			  const Conic_arc_2<CfNT, Kernel> & arc)
{
  typedef typename Conic_arc_2<CfNT,Kernel>::Point_2 Point_2;

  const Point_2& source = arc.source();
  const Point_2& target = arc.target();

  os << "{" << arc.r() << "*x^2 + "
     << arc.s() << "*y^2 + "
     << arc.t() << "*xy + " 
     << arc.u() << "*x + "
     << arc.v() << "*y + "
     << arc.w() << "} :"
     << "(" << CGAL::to_double(source.x()) << "," 
     << CGAL::to_double(source.y()) << ") -> "
     << "(" << CGAL::to_double(target.x()) << "," 
     << CGAL::to_double(target.y()) << ")";

  return (os);
}
#endif // NO_OSTREAM_INSERT_CONIC_ARC_2

CGAL_END_NAMESPACE

#endif
