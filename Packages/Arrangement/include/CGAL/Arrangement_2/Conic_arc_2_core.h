// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related docmentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $$
// release_date  : $$
//
// file          : include/CGAL/Arrangement_2/Conic_arc_2.h
// package       : Arrangement (2.62)
// maintainer    : Ron Wein <wein@post.tau.ac.il>
// author(s)     : Ron Wein <wein@post.tau.ac.il>
// 
// coordinator   : Tel-Aviv University (Dan Halperin <danha@post.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_CONIC_ARC_2_CORE_H
#define CGAL_CONIC_ARC_2_CORE_H

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Conic_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Arrangement_2/Conic_arc_2_eq_core.h>

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
 * CfNT_ is a number type for representing the conic coefficients (integers).
 * Kernel_ is a geometric kernel, where Kernel_::FT is the number type for the
 * coordinates of points (which are algebraic numbers).
 */

static int _conics_count = 0;

template <class CfNT_, class Kernel_> class Arr_conic_traits_2;

template <class CfNT_, class Kernel_>
class Conic_arc_2
{
protected:
  typedef Conic_arc_2<CfNT_, Kernel_>  Self;
        
public:
  typedef CfNT_                        CfNT;
  typedef Kernel_                      Kernel;
  typedef typename Kernel::FT          CoNT;
    
  typedef Point_ex_2<Kernel>           Point_2;
  typedef Segment_2<Kernel>            Segment_2;
  typedef Circle_2<Kernel>             Circle_2;
  
 protected:

  friend class Arr_conic_traits_2<CfNT, Kernel>;

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
    IS_CIRCLE = 128,
    IS_HYPERBOLA = 256
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
                            // Bit 7     - Is the underlying curve a circle.
                            // Bit 8     - Is the underlying curve a hyperbola.

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

  // In case of a circle, is is convinient to store its center and radius.
  struct Circular_arc_data
  {
    CoNT       x0;
    CoNT       y0;
    CoNT       r;
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

  union
  {
    Hyperbolic_arc_data *hyper_P;
    Circular_arc_data   *circ_P;
  } _data;

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
    _target(target)
  {
    _info = (arc._info & DEGREE_MASK) | 
      X_MON_UNDEFINED;

    // Check whether the conic is x-monotone.
    if (is_x_monotone())
    {
      _info = (_info & ~X_MON_UNDEFINED) | X_MONOTONE;

      // Check whether the facing information is set for the orginating arc.
      // If it is, just copy it - otherwise calculate it if the degree is 2.
      Comparison_result facing_res = arc.facing();

      if (facing_res != EQUAL)
	_info = _info | (facing_res == LARGER ? FACING_UP : FACING_DOWN);
      else if ((_info & DEGREE_MASK) == DEGREE_2)
	_set_facing();
    }
    else
    {
      _info = (_info & ~X_MON_UNDEFINED);
    }

    // Copy the hyperbolic or circular data, if necessary.
    if ((arc._info & IS_HYPERBOLA) != 0)
    {
      _info = _info | IS_HYPERBOLA;
      _data.hyper_P = new Hyperbolic_arc_data (*(arc._data.hyper_P));
    }
    else if ((arc._info & IS_CIRCLE) != 0)
    {
      _info = _info | IS_CIRCLE;
      _data.circ_P = new Circular_arc_data (*(arc._data.circ_P));
    }
    else
    {
      _data.hyper_P = NULL;
    }
  }

 public:

  /*!
   * Default constructor.
   */
  Conic_arc_2 () :
    _r(0), _s(0), _t(0), _u(0), _v(0), _w(0),
    _orient(CGAL::COLLINEAR),
    _conic_id(0),
    _info(X_MON_UNDEFINED)
  {
    _data.hyper_P = NULL;
  }

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
    _info(arc._info)
  {
    // Copy the hyperbolic or circular data, if necessary.
    if ((arc._info & IS_HYPERBOLA) != 0)
    {
      _data.hyper_P = new Hyperbolic_arc_data (*(arc._data.hyper_P));
    }
    else if ((arc._info & IS_CIRCLE) != 0)
    {
      _data.circ_P = new Circular_arc_data (*(arc._data.circ_P));
    }
    else
    {
      _data.hyper_P = NULL;
    }
  }

  /*! 
   * Construct a conic arc which lies on the conic:
   *   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
   * \param orient The orientation of the arc (clockwise or couterclockwise).
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
    _info(X_MON_UNDEFINED)
  {
    // Make sure the conic contains the two end-points on its boundary.
    CGAL_precondition(_conic_has_on_boundary(source));
    CGAL_precondition(_conic_has_on_boundary(target));

    // Make sure that the source and the taget are not the same.
    CGAL_precondition(source != target);      

    // Set the arc properties (no need to compute the orientation).
    _orient = orient;
    _set (false);
  }

  /*! 
   * Construct a conic arc which is the full conic:
   *   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
   * \pre The conic C must be an ellipse (so 4rs - t^2 > 0).
   */
  Conic_arc_2 (const CfNT& r, const CfNT& s, const CfNT& t,
	       const CfNT& u, const CfNT& v, const CfNT& w) :
    _r(r), _s(s), _t(t), _u(u), _v(v), _w(w)
  {
    // Make sure the given curve is an ellipse.
    CGAL_precondition(CGAL_NTS compare(4*r*s - t*t, CfNT(0)) == LARGER);
        
    // Set the arc to be the full conic (and compute the orientation).
    _set_full (true);
  }

  /*!
   * Construct a conic arc from the given line segment, specified by its
   * end-points (x1,y1) and (x2, y2).
   */
  Conic_arc_2 (const CfNT& x1, const CfNT& y1,
	       const CfNT& x2, const CfNT& y2) :
    _source(CoNT(x1),CoNT(y1)),
    _target(CoNT(x2),CoNT(y2))
  {
    const CfNT _zero = 0;
    const CfNT _one = 1;

    // Make sure that the source and the taget are not the same.
    CGAL_precondition(_source != _target);      

    // The supporting conic is r=s=t=0, and u*x + v*y + w = 0 should hold
    // for both the source (x1,y1) and the target (x2, y2).
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
   * Set a circular arc that lies on the given circle:
   *   C: (x - x0)^2 + (y - y0)^2 = R^2
   * \param orient The orientation of the circle.
   * \param source The source point.
   * \param target The target point.
   * \pre The source and the target must be on the conic boundary and must
   * not be the same.
   */  
  Conic_arc_2 (const CfNT& x0, const CfNT& y0,
	       const CfNT& R,
	       const Orientation& orient,
	       const Point_2& source, const Point_2& target) :
    _source(source),
    _target(target),
    _info(X_MON_UNDEFINED)
  {
    // Produce the correponding conic: if the circle centre is (x0,y0)
    // and it radius is R, that its equation is:
    //   x^2 + y^2 - 2*x0*x - 2*y0*y + (x0^2 + y0^2 - R^2) = 0
    // Since this equation describes a curve with a negative (clockwise) 
    // orientation, we multiply it by -1 if necessary to obtain a positive
    // (counterclockwise) orientation.
    const CfNT _zero = 0;

    if (orient == CGAL::COUNTERCLOCKWISE)
    {
      const CfNT _minus_one = -1;
      const CfNT _two       = 2;

      _r = _minus_one;
      _s = _minus_one;
      _t = _zero;
      _u = _two*x0;
      _v = _two*y0;
      _w = R*R - x0*x0 - y0*y0;

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
      _w = x0*x0 + y0*y0 - R*R;

      _orient = CGAL::CLOCKWISE;
    }

    // Make sure the circle contains the two endpoints on its boundary.
    CGAL_precondition(_conic_has_on_boundary(source));
    CGAL_precondition(_conic_has_on_boundary(target));

    // Make sure that the source and the taget are not the same.
    CGAL_precondition(source != target);      

    // Prepare the auxiliary data structure.
    Circular_arc_data    *circ_data_P = new Circular_arc_data;

    circ_data_P->x0 = CoNT(x0);
    circ_data_P->y0 = CoNT(y0);
    circ_data_P->r = CoNT(R);

    // Set the arc properties (no need to compute the orientation).
    _set (false, circ_data_P);
  }

  /*!
   * Set a circular arc that corresponds to the full circle:
   *   C: (x - x0)^2 + (y - y0)^2 = R^2
   */ 
  Conic_arc_2 (const CfNT& x0, const CfNT& y0,
	       const CfNT& R)
  {
    // Produce the correponding conic: if the circle centre is (x0,y0)
    // and it radius is R, that its equation is:
    //   x^2 + y^2 - 2*x0*x - 2*y0*y + (x0^2 + y0^2 - R^2) = 0
    // Note that this equation describes a curve with a negative (clockwise) 
    // orientation.
    const CfNT _zero = 0;
    const CfNT _one       = 1;
    const CfNT _minus_two = -2;

    _r = _one;
    _s = _one;
    _t = _zero;
    _u = _minus_two*x0;
    _v = _minus_two*y0;
    _w = x0*x0 + y0*y0 - R*R;

    _orient = CGAL::CLOCKWISE;

    // Prepare the auxiliary data structure.
    Circular_arc_data    *circ_data_P = new Circular_arc_data;

    circ_data_P->x0 = CoNT(x0);
    circ_data_P->y0 = CoNT(y0);
    circ_data_P->r = CoNT(R);

    // Set the arc to be the full conic (no need to compute the orientation).
    _set_full (false, circ_data_P);
  }

  /*!
   * Construct a circular arc from the given three points, specified by the
   * coordinates (x1,y1), (x2,y2) and (x3, y3).
   * (x1,y1) is the arc source and (x3,y3) is the arc target.
   * \pre The three points must not be collinear.
   */
  Conic_arc_2 (const CfNT& x1, const CfNT& y1,
	       const CfNT& x2, const CfNT& y2,
	       const CfNT& x3, const CfNT& y3,
	       const bool& ) :
    _source(CoNT(x1),CoNT(y1)),
    _target(CoNT(x3),CoNT(y3)),
    _info(X_MON_UNDEFINED)
  {
    const CfNT _zero = 0;
    const CfNT _two  = 2;
    
    // Make sure that the source and the taget are not the same.
    CGAL_precondition(_source != _target);      

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
    CGAL_precondition_code(
    const bool points_collinear = (D == _zero);
    );
    CGAL_precondition(!points_collinear);

    // The equation of the underlying circle is given by:
    _r = D*D;
    _s = D*D;
    _t = _zero;
    _u = -_two*D*Nx;
    _v = -_two*D*Ny;
    _w = Nx*Nx + Ny*Ny - ((D*x2 - Nx)*(D*x2 - Nx) + (D*y2 - Ny)*(D*y2 - Ny));

    // Determine the orientation: If the mid-point forms a left-turn with
    // the source and the target points, the orientation is positive (going
    // counterclockwise).
    // Otherwise, it is negative (going clockwise).
    static Kernel                  ker;
    typename Kernel::Orientation_2 orient_f = ker.orientation_2_object();
    Point_2                        pmid = Point_2(CoNT(x2), CoNT(y2));
    
    if (orient_f(_source, pmid, _target) == LEFT_TURN)
      _orient = CGAL::COUNTERCLOCKWISE;
    else
      _orient = CGAL::CLOCKWISE;

    // Prepare the auxiliary data structure.
    Circular_arc_data    *circ_data_P = new Circular_arc_data;

    circ_data_P->x0 = CoNT(Nx) / CoNT(D);
    circ_data_P->y0 = CoNT(Ny) / CoNT(D);
    circ_data_P->r = 
      CGAL::sqrt (CoNT((D*x2 - Nx)*(D*x2 - Nx) + (D*y2 - Ny)*(D*y2 - Ny)) /
		  CoNT(D*D));

    // Set the arc properties (no need to compute the orientation).
    _set (false, circ_data_P);
  }

  /*!
   * Destructor.
   */
  virtual ~Conic_arc_2 ()
  {
    if ((_info & IS_HYPERBOLA) != 0)
      delete _data.hyper_P;
    else if ((_info & IS_CIRCLE) != 0)
      delete _data.circ_P;
    _data.hyper_P = NULL;
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
    if ((_info & IS_HYPERBOLA) != 0)
      delete _data.hyper_P;
    else if ((_info & IS_CIRCLE) != 0)
      delete _data.circ_P;
    _data.hyper_P = NULL;

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

    // Duplicate the data for hyperbolic or circular arcs.
    if ((arc._info & IS_HYPERBOLA) != 0)
    {
      _data.hyper_P = new Hyperbolic_arc_data (*(arc._data.hyper_P));
    }
    else if ((arc._info & IS_CIRCLE) != 0)
    {
      _data.circ_P = new Circular_arc_data (*(arc._data.circ_P));
    }

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
    return ((_info & IS_CIRCLE) != 0);
  }
  
  /*!
   * Check whether the arc is a line segment.
   * \return (true) if the underlying conic curve is of degree 1.
   */
  bool is_segment () const
  {
    return ((_info & DEGREE_MASK) == 1);
  }

  /*!
   * Check whether the arc is a vertical segment.
   * \return (true) if the arc is a vertical segment.
   */
  bool is_vertical_segment () const
  {
    // A vertical segment is contained in the degenerate conic: u*x + w = 0.
    const CfNT _zero = 0;

    return ((_info & DEGREE_MASK) == 1 && 
	    CGAL_NTS compare(_v, _zero) == EQUAL);
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
    
    for (int i = 0; i < n; i++)
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

    opp_arc._r = -r;
    opp_arc._s = -s;
    opp_arc._t = -t;
    opp_arc._u = -u;
    opp_arc._v = -v;
    opp_arc._w = -w;
    
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

    if ((_info & IS_HYPERBOLA) != 0)
    {
      opp_arc._data.hyper_P = new Hyperbolic_arc_data (*_data.hyper_P);
    }
    else if ((_info & IS_CIRCLE) != 0)
    {
      opp_arc._data.circ_P = new Circular_arc_data (*_data.circ_P);
    }
    else
    {
      opp_arc._data.hyper_P = NULL;
    }
 
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
      if (sl_denom > 0)
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

    if (sl_denom > 0) // so sl2_denom > 0 as well ...
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
      if (sl_denom > 0)
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

    if (sl_denom > 0) // so sl2_denom > 0 as well ...
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
   * \pre The tow arcs do not lie on the same conic.
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

      if (_s == _zero && arc._s != _zero)
      {
	n_xs = arc._x_coordinates_of_intersections_with (*this,
							 xs);
      }
      else
      {
	n_xs = _x_coordinates_of_intersections_with (arc,
						     xs);
      }

      if (_r == _zero && arc._r != _zero)
      {
	n_ys = arc._y_coordinates_of_intersections_with (*this,
							 ys);
      }
      else
      {
	n_ys = _y_coordinates_of_intersections_with (arc,
						     ys);
      }
    
      // Perform the pairing process od the x and y coordinates.
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
   * \param circ_data_P The center and radius of the base circle
   *                    (only if the arc lies on a circle).
   */
  void _set (const bool& comp_orient,
	     Circular_arc_data *circ_data_P = NULL)
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
 
    if (_r != _zero || _s != _zero || _t != _zero)
    {
      // In case one of the coefficients of x^2,y^2 or xy is not zero, the
      // degree is 2.
      deg = 2;
      CGAL_assertion (_orient != CGAL::COLLINEAR);
    }
    else if (_u != _zero || _v != _zero)
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
    if (deg == 2 && 4*_r*_s < _t*_t)
    {
      _info = _info | IS_HYPERBOLA;
      _build_hyperbolic_arc_data ();
    }
    // In case the base conic is a circle, set the circular data.
    else if (deg == 2 && circ_data_P != NULL)
    {
      _info = _info | IS_CIRCLE;
      _data.circ_P = circ_data_P;
    }
    else
    {
      _data.hyper_P = NULL;
    }

    // In case of a non-degenerate parabola or a hyperbola, make sure 
    // the arc is not infinite.
    if (deg == 2 && 4*_r*_s <= _t*_t)
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
      if ((_info & DEGREE_MASK) == DEGREE_2)
	_set_facing();
    }
    else
    {
      _info = (_info & ~X_MON_UNDEFINED);
    }

    return;
  }

  /*!
   * Set the properties of a conic arc that is really a full curve
   * (that is, an ellipse).
   * \param comp_orient Should we compute the orientation of the given curve.
   * \param circ_data_P The center and radius of the base circle
   *                    (only if the arc is a full circle).
   */
  void _set_full (const bool& comp_orient,
		  Circular_arc_data *circ_data_P = NULL)
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

    // Make sure the conic is a non-degenerate ellipse:
    // The coefficients should satisfy (4rs - t^2) > 0.
    CGAL_precondition(4*_r*_s > _t*_t);

    // Set the information: a full conic, which is obvoiusly not x-monotone.
    _info = DEGREE_2 | FULL_CONIC;

    // Check if the conic is really a circle.
    if (circ_data_P != NULL)
    {
      _info = _info | IS_CIRCLE;
      _data.circ_P = circ_data_P;
    }
    else
    {
      _data.hyper_P = NULL;
    }

    // Assign one of the vertical tangency points as both the source and
    // the target of the conic arc.
    Point_2    vpts[2];
    int        n_vpts;

    if (circ_data_P != NULL)
    {
      vpts[0] = Point_2 (CoNT(circ_data_P->x0 + circ_data_P->r), 
			 CoNT(circ_data_P->y0));
    }
    else
    {
      n_vpts = _conic_vertical_tangency_points (vpts);

      CGAL_assertion(n_vpts > 0);
      CGAL_assertion(_conic_has_on_boundary(vpts[0]));
    }

    // If we reached here, the conic arc is legal: Get a new id for the conic.
    _conic_id = _get_new_conic_id();

    _source.set_generating_conics (_conic_id);
    _target.set_generating_conics (_conic_id);

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
    if (t == _zero)
    {
      // sin(2*phi) == 0, so phi = 0 or phi = PI/2
      if (cos_2phi > _zero)
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
    else if (t > _zero)
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

    CGAL_assertion (det < _zero);
    
    x0 = (t*v - _two*s*u) / det;
    y0 = (t*u - _two*r*v) / det;
    
    // The axis separating the two branches of the hyperbola is now given by:
    // 
    //  cos(phi)*x + sin(phi)*y - (cos(phi)*x0 + sin(phi)*y0) = 0
    //
    _data.hyper_P = new Hyperbolic_arc_data;

    _data.hyper_P->a = cos_phi;
    _data.hyper_P->b = sin_phi;
    _data.hyper_P->c = - (cos_phi*x0 + sin_phi*y0);

    // Make sure that the two endpoints are located on the same branch
    // of the hyperbola.
    _data.hyper_P->side = _hyperbolic_arc_side(_source);

    CGAL_assertion (_data.hyper_P->side = _hyperbolic_arc_side(_target));

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
    if (_data.hyper_P == NULL)
      return (0);

    CoNT       val;

    val = _data.hyper_P->a*p.x() + _data.hyper_P->b*p.y() + _data.hyper_P->c;
    return ((val > 0) ? 1 : -1);
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
    if ((_info & IS_HYPERBOLA) != 0)
    {
      if (_hyperbolic_arc_side(p) != _data.hyper_P->side)
	return (false);
    }

    // Act according to the conic degree.
    if ((_info & DEGREE_MASK) == DEGREE_1)
    {
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
      static Kernel                  ker;
      typename Kernel::Orientation_2 orient_f = ker.orientation_2_object();

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
    val = CoNT(_r)*x*x + CoNT(_s)*y*y + CoNT(_t)*x*y +
          CoNT(_u)*x + CoNT(_v)*y + CoNT(_w);

    return (val == _zero);
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
    if ((_info & DEGREE_MASK) == DEGREE_1 || _s == _zero)
      return (0);

    // Special treatment for circles, where the vertical tangency points
    // are simply (x0-r,y0) and (x0+r,y0).
    if ((_info & IS_CIRCLE) != 0)
    {
      ps[0] = Point_2 (CoNT(_data.circ_P->x0 - _data.circ_P->r), 
		       CoNT(_data.circ_P->y0),
		       _conic_id);

      ps[1] = Point_2 (CoNT(_data.circ_P->x0 + _data.circ_P->r), 
                       CoNT(_data.circ_P->y0),
		       _conic_id);

      return (2);
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

    if (_t == _zero)
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
	  if (ys[j] == -(CoNT(_t)*xs[i] + CoNT(_v)) / CoNT(_two*_s))
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
    if ((_info & DEGREE_MASK) == DEGREE_1 || _r == _zero)
      return (0);

    // Special treatment for circles, where the horizontal tangency points
    // are simply (x0,y0-r) and (x0,y0+r).
    if ((_info & IS_CIRCLE) != 0)
    {
      ps[0] = Point_2 (CoNT(_data.circ_P->x0), 
                       CoNT(_data.circ_P->y0 - _data.circ_P->r),
		       _conic_id);

      ps[1] = Point_2 (CoNT(_data.circ_P->x0), 
                       CoNT(_data.circ_P->y0 + _data.circ_P->r),
		       _conic_id);

      return (2);
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
   * Calculate all x coordinates of intersection points between the two
   * base conics of (*this) and the given arc.
   * \param arc The arc whose underlying conic we intersect.
   * \param xs The output x coordinates.
   *           This area must be allocated to the size of 4.
   * \return The number of unique x coordinates.
   */
  int _x_coordinates_of_intersections_with (const Self& arc,
					    CoNT* xs) const
  {
    const CfNT _zero = 0;
    int        n_roots;         // The number of distinct x values.

    // Check whether both arcs are line segments.
    if ((_info & DEGREE_MASK) == DEGREE_1 && 
	(arc._info & DEGREE_MASK) == DEGREE_1)
    {
      // The two conics are: u*x + v*y + w = 0
      //                and: u'*x + v'*y + w' = 0
      // There's a single solution for x, which is:
      const CfNT denom = _v*arc._w - _w*arc._v;
      const CfNT numer = _u*arc._v - _v*arc._u;

      if (numer == _zero)
      {
	n_roots = 0;
      }
      else
      {
	xs[0] = CoNT(denom) / CoNT(numer);
	n_roots = 1;
      }
    }
    // Check whether the second arc is really a line segment.
    else if ((arc._info & DEGREE_MASK) == DEGREE_1)
    {
      // The two conics are: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
      //                and: a*x + b*y + c = 0
      // There are therefore 2 possible x-values, the solutions for:
      const CfNT  _two = 2;
      const CfNT& a = arc._u;
      const CfNT& b = arc._v;
      const CfNT& c = arc._w;

      n_roots = solve_quadratic_eq<CfNT, CoNT> 
	(a*a*_s + b*b*_r - a*b*_t,
	 _two*a*c*_s + b*b*_u - a*b*_v - b*c*_t,
	 c*c*_s + b*b*_w - b*c*_v,
	 xs);
    }
    // Check if the two arcs are circular arcs.
    else if ((_info & IS_CIRCLE) != 0 && (arc._info & IS_CIRCLE) != 0)
    {
      // Special treatment for two circles.
      // The two curves are: r*x^2 + r*y^2 + u*x + v*y + w = 0
      //                and: r'*x^2 + r'*y^2 + u'*x + v'*y + w' = 0
      //
      // Thus, r'*C1-r*C2 is a line whose equation is: a*x + b*y + = 0, where:
      const CfNT  _two = 2;
      const CfNT a = arc._r*_u - _r*arc._u;
      const CfNT b = arc._r*_v - _r*arc._v;
      const CfNT c = arc._r*_w - _r*arc._w;

      if (b == _zero)
      {
	// The line a*x + c = 0 connects both intersection points of the two
	// circles, so the both have an x-coordinate of -c/a.
	if (a == _zero)
	{
	  n_roots = 0;
	}
	else
	{
	  xs[0] = -CoNT(c) / CoNT(a);
	  n_roots = 1;
	}
      }
      else
      {
	// The intersection points of the two circles are the same as the
	// intersection points of one of the circles with a*x + b*y + c = 0.
	n_roots = solve_quadratic_eq<CfNT,CoNT> 
	  (a*a*_s + b*b*_r - a*b*_t,
	   _two*a*c*_s + b*b*_u - a*b*_v - b*c*_t,
	   c*c*_s + b*b*_w - b*c*_v,
	   xs);      
      }
    }
    // Check if the two arcs lie of canonic parabolas.
    else if (_s == _zero && _t == _zero &&
	     arc._s == _zero && arc._t == _zero)
    {
      // Special treatment for canonic parabolas whose axes are parallel
      // to the y axis.
      // The two curves are: r*x^2 + u*x + v*y + w = 0
      //                and: r'*x^2 + u'*x + v'*y + w' = 0
      // There are therefore 2 possible x-values, the solutions for:
      n_roots = solve_quadratic_eq<CfNT,CoNT> 
	(_r*arc._v - _v*arc._r, 
	 _u*arc._v - _v*arc._u, 
	 _w*arc._v - _v*arc._w,
	 xs);
    }
    else
    {
      // If the two conics are: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
      //                   and: r'*x^2 + s'*y^2 + t'*xy + u'*x + v'*y + w' = 0
      // Let us define:
      // A = s'*r - s*r'    if (s' != 0), otherwise   A = r'
      // B = s'*u - s*u'    if (s' != 0), otherwise   B = u'
      // C = s'*w - s*w'    if (s' != 0), otherwise   C = w'
      // D = s'*t - s*t'    if (s' != 0), otherwise   A = t'
      // E = s'*v - s*v'    if (s' != 0), otherwise   A = v'
      CfNT    A, B, C, D, E;

      if (arc._s != _zero)
      {
	const CfNT& s_tag = arc._s;

	A = s_tag*_r - _s*arc._r;
	B = s_tag*_u - _s*arc._u;
	C = s_tag*_w - _s*arc._w;
	D = s_tag*_t - _s*arc._t;
	E = s_tag*_v - _s*arc._v;
      }
      else
      {
	A = arc._r;
	B = arc._u;
	C = arc._w;
	D = arc._t;
	E = arc._v;
      }

      if (D == _zero && E == _zero)
      {
	// In this case: A*x^2 + B*x + C = 0, so:
	n_roots = solve_quadratic_eq<CfNT,CoNT> (A, B, C,
						 xs);      
      }
      else
      {
	// Now we have:
	//
	//         A*x^2 + B*x + C
	//  y = - -----------------
	//              D*x + E
	//
	// Applying this two our conic's equation yields a quartic equation:
	//
	//  c[4]*x^4 + c[3]*x^3 + c[2]*x^2 + c[1]*x + c[0] = 0
	const CfNT _two = 2;
	CfNT       cs[5];
	
	if (_t == _zero && arc._t == _zero)
	{
	  cs[4] = _s*A*A;
	  cs[3] = _two*_s*A*B;
	  cs[2] = _r*E*E + _two*_s*A*C + _s*B*B - _v*A*E;
	  cs[1] = _u*E*E + _two*_s*B*C - _v*B*E;
	  cs[0] = _w*E*E + _s*C*C - _v*C*E;
	}
	else
	{
	  const CfNT F = _t*E + _v*D;
	
	  cs[4] = _r*D*D + _s*A*A - _t*A*D;
	  cs[3] = _two*_r*D*E + _u*D*D + _two*_s*A*B - _t*B*D - F*A;
	  cs[2] = _r*E*E + _two*_u*D*E + _w*D*D + 
	          _two*_s*A*C + _s*B*B - _t*C*D - F*B - _v*A*E;
	  cs[1] = _u*E*E + _two*_w*D*E + _two*_s*B*C - F*C - _v*B*E;
	  cs[0] = _w*E*E + _s*C*C - _v*C*E;
	}

	// Solve the quartic equation.
	n_roots = solve_quartic_eq<CfNT,CoNT> 
	  (cs[4], cs[3], cs[2], cs[1], cs[0],
	   xs);
      }
   }

    return (n_roots);
  }

  /*!
   * Calculate all y coordinates of intersection points between the two
   * base conics of (*this) and the given arc.
   * \param arc The arc whose underlying conic we intersect.
   * \param ys The output y coordinates.
   *           This area must be allocated to the size of 4.
   * \return The number of unique y coordinates.
   */
  int _y_coordinates_of_intersections_with (const Self& arc,
					    CoNT* ys) const
  {
    const CfNT _zero = 0;
    int        n_roots;         // The number of distinct y values.

    // Check whether both arcs are line segments.
    if ((_info & DEGREE_MASK) == DEGREE_1 && 
	(arc._info & DEGREE_MASK) == DEGREE_1)
    {
      // The two conics are: u*x + v*y + w = 0
      //                and: u'*x + v'*y + w' = 0
      // There's a single solution for y, which is:
      const CfNT denom = _u*arc._w - _w*arc._u;
      const CfNT numer = _v*arc._u - _u*arc._v;

      if (numer == _zero)
      {
	n_roots = 0;
      }
      else
      {
	ys[0] = CoNT(denom) / CoNT(numer);
	n_roots = 1;
      }
    }
    // Check whether the second arc is really a line segment.
    else if ((arc._info & DEGREE_MASK) == DEGREE_1)
    {
      // The two conics are: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
      //                and: a*x + b*y + c = 0
      // There are therefore 2 possible y-values, the solutions for:
      const CfNT  _two = 2;
      const CfNT& a = arc._u;
      const CfNT& b = arc._v;
      const CfNT& c = arc._w;

      // Solve the equation.
      n_roots = solve_quadratic_eq<CfNT,CoNT> 
	(b*b*_r + a*a*_s - a*b*_t,
	 _two*b*c*_r + a*a*_v - a*b*_u - a*c*_t,
	 c*c*_r + a*a*_w - a*c*_u,
	 ys);
    }
    // Check if the two arcs are circular arcs.
    else if ((_info & IS_CIRCLE) != 0 && (arc._info & IS_CIRCLE) != 0)
    {
      // Special treatment for two circles.
      // The two curves are: r*x^2 + r*y^2 + u*x + v*y + w = 0
      //                and: r'*x^2 + r'*y^2 + u'*x + v'*y + w' = 0
      //
      // Thus, r'*C1-r*C2 is a line whose equation is: a*x + b*y + = 0, where:
      const CfNT _two = 2;
      const CfNT a = arc._r*_u - _r*arc._u;
      const CfNT b = arc._r*_v - _r*arc._v;
      const CfNT c = arc._r*_w - _r*arc._w;

      if (a == _zero)
      {
	// The line b*y + c = 0 connects both intersection points of the two
	// circles, so the both have a y-coordinate of -c/b.
	if (b == _zero)
	{
	  n_roots = 0;
	}
	else
	{
	  ys[0] = -CoNT(c) / CoNT(b);
	  n_roots = 1;
	}
      }
      else
      {
	// The intersection points of the two circles are the same as the
	// intersection points of one of the circles with a*x + b*y + c = 0.
	// Set the generating polynomial.
	n_roots = solve_quadratic_eq<CfNT,CoNT> 
	  (b*b*_r + a*a*_s - a*b*_t,
	   _two*b*c*_r + a*a*_v - a*b*_u - a*c*_t,
	   c*c*_r + a*a*_w - a*c*_u,
	   ys);
      }
    }
    // Check if the two arcs lie of canonic parabolas.
    else if (_r == _zero && _t == _zero &&
	     arc._r == _zero && arc._t == _zero)
    {
      // Special treatment for canonic parabolas whose axes are parallel
      // to the x axis.
      // The two curves are: s*y^2 + u*x + v*y + w = 0
      //                and: s'*y^2 + u'*x + v'*y + w' = 0
      // There are therefore 2 possible y-values, the solutions for:
      n_roots = solve_quadratic_eq<CfNT,CoNT> 
	(_s*arc._u - _u*arc._s, 
	 _v*arc._u - _u*arc._v, 
	 _w*arc._u - _u*arc._w,
	 ys);
    }
    else
    {
      // If the two conics are: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
      //                   and: r'*x^2 + s'*y^2 + t'*xy + u'*x + v'*y + w' = 0
      // Let us define:
      // A = r'*s - r*s'    if (r' != 0), otherwise   A = s'
      // B = r'*v - r*v'    if (r' != 0), otherwise   B = v'
      // C = r'*w - r*w'    if (r' != 0), otherwise   C = w'
      // D = r'*t - r*t'    if (r' != 0), otherwise   A = t'
      // E = r'*u - r*u'    if (r' != 0), otherwise   A = u'
      CfNT    A, B, C, D, E;

      if (arc._r != _zero)
      {
	const CfNT& r_tag = arc._r;

	A = r_tag*_s - _r*arc._s;
	B = r_tag*_v - _r*arc._v;
	C = r_tag*_w - _r*arc._w;
	D = r_tag*_t - _r*arc._t;
	E = r_tag*_u - _r*arc._u;
      }
      else
      {
	A = arc._s;
	B = arc._v;
	C = arc._w;
	D = arc._t;
	E = arc._u;
      }

      if (D == _zero && E == _zero)
      {
	// In this case: A*y^2 + B*y + C = 0, so:
	n_roots = solve_quadratic_eq<CfNT,CoNT> (A, B, C,
						 ys);
      }
      else
      {
	// Now we have:
	//
	//         A*y^2 + B*y + C
	//  x = - -----------------
	//              D*y + E
	//
	// Applying this two our conic's equation yields a quartic equation:
	//
	//  c[4]*y^4 + c[3]*y^3 + c[2]*y^2 + c[1]*y + c[0] = 0
	const CfNT _two = 2;
	CfNT       cs[5];

	if (_t == _zero && arc._t == _zero)
	{
	  cs[4] = _r*A*A;
	  cs[3] = _two*_r*A*B;
	  cs[2] = _s*E*E + _two*_r*A*C + _r*B*B - _u*A*E;
	  cs[1] = _v*E*E + _two*_r*B*C - _u*B*E;
	  cs[0] = _w*E*E + _r*C*C - _u*C*E;
	}
	else
	{
	  const CfNT F = _t*E + _u*D;
	
	  cs[4] = _s*D*D + _r*A*A - _t*A*D;
	  cs[3] = _two*_s*D*E + _v*D*D + _two*_r*A*B - _t*B*D - F*A;
	  cs[2] = _s*E*E + _two*_v*D*E + _w*D*D + 
	          _two*_r*A*C + _r*B*B - _t*C*D - F*B - _u*A*E;
	  cs[1] = _v*E*E + _two*_w*D*E + _two*_r*B*C - F*C - _u*B*E;
	  cs[0] = _w*E*E + _r*C*C - _u*C*E;
	}

	n_roots = solve_quartic_eq<CfNT,CoNT> 
	  (cs[4], cs[3], cs[2], cs[1], cs[0],
	   ys);
      }
    }

    return (n_roots);
  }

  /*! 
   * Pair the x coordinates and the y coordinates of the intersection point
   * of the underlying conics of (*this) and arc, and return a vector of 
   * intersection points.
   * \param arc The other arc.
   * \param xs The potentail x coordinates.
   * \param n_xs Number of x coordinates.
   * \param ys The potentail y coordinates.
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
	// If the current pair of x and y coordinates lies on both underlying,
	// accept it.
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

public:

  /*!
   * Get a segment (if the arc is indeed a line segment).
   * \return A segment equivalent ot the arc.
   * \pre The conic arc is indeed a line segment (it is of degree 1).
   */
  Segment_2 segment() const
  {
    CGAL_precondition(is_segment());

    return (Segment_2 (_source, _target));
  }

  /*!
   * Get the underlying circle (if the arc is indeed a circular arc).
   * \return The underlying circle.
   * \pre The conic arc is indeed a circular arc.
   */
  Circle_2 circle() const
  {
    CGAL_precondition(is_circular());

    // Create the appropriate circle.
    const CoNT _zero = 0;
    const CoNT _two = 2;
    CoNT       x0, y0, r2;

    if (_r > CfNT(0))
    {
      // Positive orientation. The conic has the form:
      //  x^2 + y^2 - (2*x0)*x - (2*y0)*y + (x0^2 + y0^2 - r^2) = 0 
      x0 = -CoNT(_u) / _two;
      y0 = -CoNT(_v) / _two;
      r2 = x0*x0 + y0*y0 - CoNT(_w);
    }
    else
    {
      // Negative orientation:
      //  - x^2 - y^2 + (2*x0)*x + (2*y0)*y + (r^2 - x0^2 - y0^2) = 0 
      x0 = CoNT(_u) / _two;
      y0 = CoNT(_v) / _two;
      r2 = CoNT(_w) - (x0*x0 + y0*y0);
    }

    return (Circle_2 (Point_2(x0, y0), r2));
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

  os << "{" << CGAL::to_double(arc.r()) << "*x^2 + "
     << CGAL::to_double(arc.s()) << "*y^2 + "
     << CGAL::to_double(arc.t()) << "*xy + " 
     << CGAL::to_double(arc.u()) << "*x + "
     << CGAL::to_double(arc.v()) << "*y + "
     << CGAL::to_double(arc.w()) << "} :"
     << "(" << CGAL::to_double(source.x()) << "," 
     << CGAL::to_double(source.y()) << ") -> "
     << "(" << CGAL::to_double(target.x()) << "," 
     << CGAL::to_double(target.y()) << ")";

  return (os);
}
#endif // NO_OSTREAM_INSERT_CONIC_ARC_2

CGAL_END_NAMESPACE

#endif
