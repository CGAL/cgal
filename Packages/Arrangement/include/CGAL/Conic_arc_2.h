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
// file          : include/CGAL/Conic_arc_2.h
// package       : Arrangement (2.19)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Ron Wein <wein@post.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_CONIC_ARC_2_H
#define CGAL_CONIC_ARC_2_H

#include <CGAL/basic.h>
#include <list>

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Conic_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Bbox_2.h>
#include <fstream>

#include <CGAL/Conic_arc_2_eq.h>

#define CGAL_CONIC_ARC_USE_BOUNDING_BOX
#define CGAL_CONIC_ARC_USE_CACHING

CGAL_BEGIN_NAMESPACE

enum
{
  REFLECT_IN_X = 1,
  REFLECT_IN_Y = 2,
  REFLECTION_FACTOR = 4
};

template <class R>
class Point_2_ex : public R::Point_2
{
public:
  typedef typename R::Point_2 Base;
  typedef typename R::RT      RT;

  typedef enum 
  {
    User_defined,
    Tangency,
    Intersection_exact,
    Intersection_approx,
    Ray_shooting_exact,
    Ray_shooting_approx
  } Type;

private:
  Type       _type;
  int        conic_id1;
  int        conic_id2;

public:

  // Constructors.
  Point_2_ex () :
    Base(),
    _type(User_defined),
    conic_id1(0),
    conic_id2(0)
  {}

  Point_2_ex (const RT & hx, const RT & hy, const RT& hz, const Type & type) :
    Base(hx,hy,hz),
    _type(type),
    conic_id1(0),
    conic_id2(0)
  {}

  Point_2_ex (const RT & hx, const RT & hy, const Type & type,
              const int & id1 = 0, const int & id2 = 0) :
    Base(hx,hy),
    _type(type),
    conic_id1(id1),
    conic_id2(id2)
  {}

  Point_2_ex (const RT& hx, const RT& hy) :
    Base(hx,hy),
    _type(User_defined),
    conic_id1(0),
    conic_id2(0)
  {}

  // Check if the point is approximate.
  bool is_approximate () const
  {
    return (_type == Intersection_approx || _type == Ray_shooting_approx);
  }

  // Check if the given conic generates the points.
  bool is_generating_conic_id (const int& id) const
  {
    return (id == conic_id1 || id == conic_id2);
  }

  // Compare the co-ordinates of two given points.
  Comparison_result compare_x (const Point_2_ex<R>& p) const
  {
    if (is_approximate() || p.is_approximate())
      return (eps_compare<APNT>(TO_APNT(x()), TO_APNT(p.x())));
    
    return (CGAL::compare (x(), p.x()));
  }

  Comparison_result compare_y (const Point_2_ex<R>& p) const
  {
    if (is_approximate() || p.is_approximate())
      return (eps_compare<APNT>(TO_APNT(y()), TO_APNT(p.y())));
    
    return (CGAL::compare (y(), p.y()));
  }

  bool equals (const Point_2_ex<R>& p) const
  {
    if (is_approximate() || p.is_approximate())
      return (eps_compare<APNT>(TO_APNT(x()), TO_APNT(p.x())) == EQUAL &&
              eps_compare<APNT>(TO_APNT(y()), TO_APNT(p.y())) == EQUAL);
    
    return (CGAL::compare (x(), p.x()) == EQUAL &&
            CGAL::compare (y(), p.y()) == EQUAL);
  }

  Comparison_result compare_lex_xy (const Point_2_ex<R>& p) const
  {
    Comparison_result   res = this->compare_x (p);

    if (res != EQUAL)
      return (res);
    
    return (this->compare_y (p));
  }

  // Reflect a point.
  Point_2_ex<R> reflect_in_y () const
  {
    Point_2_ex<R> ref_point (-hx(), hy(), hw(), _type);

    if (conic_id1 != 0)
      ref_point.conic_id1 = conic_id1 ^ REFLECT_IN_Y;
    if (conic_id2 != 0)
      ref_point.conic_id2 = conic_id2 ^ REFLECT_IN_Y;

    return (ref_point);
  }

  Point_2_ex<R> reflect_in_x_and_y () const
  {
    Point_2_ex<R> ref_point (-hx(), -hy(), hw(), _type);

    if (conic_id1 != 0)
      ref_point.conic_id1 = conic_id1 ^ (REFLECT_IN_X | REFLECT_IN_Y);
    if (conic_id2 != 0)
      ref_point.conic_id2 = conic_id2 ^ (REFLECT_IN_X | REFLECT_IN_Y);

    return (ref_point);
  }

};

// ----------------------------------------------------------------------------
// Representation of a conic curve.
//

static int _conics_count = 0;
template <class _NT> class Arr_conic_traits_2;

template <class NT>
class Conic_arc_2
{
  friend class Arr_conic_traits_2<NT>;

 public:

  typedef Cartesian<NT>        R;
  typedef Point_2_ex<R>        Point_2;
  typedef typename R::Conic_2  Conic_2;

  // Obsolete, for backward compatibility
  typedef Point_2              Point;
  typedef Conic_2              Conic;

 private:

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
    FACING_MASK = 32 + 64
  };

  Conic    _conic;              // The conic that contains the arc.
  int      _conic_id;           // The id of the conic.
  Point    _source;             // The source of the arc. 
  Point    _target;             // The target of the arc.
  int      _info;               // A bit array with extra information:
                                // Bit 0 & 1 - The degree of the conic 
                                //             (either 1 or 2).
                                // Bit 2     - Whether the arc is a full conic.
                                // Bit 3 & 4 - Whether the arc is x-monotone
                                //             (00, 01 or 11 - undefined).
                                // Bit 5 & 6 - Indicate whether the arc is
                                //             facing upwards or downwards (for
                                //             x-monotone curves of degree 2).
#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
  Bbox_2   _bbox;               // A bounding box for the arc.  
#endif

  // For arcs whose base is a hyperbola we store the axis (a*x + b*y + c = 0)
  // which separates the two bracnes of the hyperbola. We also store the side
  // (-1 or 1) that the arc occupies.
  struct Hyperbolic_arc_data
  {
    NT       a;
    NT       b;
    NT       c;
    int      side;
  };

#ifdef CGAL_CONIC_ARC_USE_CACHING
  struct Intersections
  {
    int   id1;
    int   id2;
    int   n_points;
    Point ps[4];
  };
#endif

  Hyperbolic_arc_data *_hyper_data_P;

  // Produce a unique id for a new conic. 
  int _get_new_conic_id ()
  {
    _conics_count++;
    return (REFLECTION_FACTOR * _conics_count);
  }

  // Private constructor.
  Conic_arc_2 (const Conic_arc_2 & arc,
               const Point & source, const Point & target,
               const bool & is_full) :
    _conic(arc._conic),
    _conic_id(arc._conic_id),
    _source(source),
    _target(target),
    _hyper_data_P(NULL)
  {
    CGAL_precondition(is_full || ! source.equals(target) );

    _info = (arc._info & DEGREE_MASK) | 
      (is_full ? FULL_CONIC : 0) | 
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

    // Copy the hyperbolic data, if necessary.
    if (arc._hyper_data_P != NULL)
    {
      _hyper_data_P = new Hyperbolic_arc_data (*(arc._hyper_data_P));
    }

#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
    _bbox = bounding_box();         // Compute the bounding box.
#endif
  }

 public:

  // Default constructor.
  Conic_arc_2 () :
    _conic_id(0),
    _info(X_MON_UNDEFINED),
    _hyper_data_P(NULL)
  {}

  // Copy constructor.
  Conic_arc_2 (const Conic_arc_2<NT>& arc) :
    _conic(arc._conic),
    _conic_id(arc._conic_id),
    _source(arc._source),
    _target(arc._target),
    _info(arc._info),
    _hyper_data_P(NULL)
  {
    if (arc._hyper_data_P != NULL)
    {
      _hyper_data_P = new Hyperbolic_arc_data (*(arc._hyper_data_P));
    }

#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
    _bbox = arc._bbox;              // Copy the bounding box.  
#endif
  }

  // Construct a conic arc.
  // The source and the target must be on the conic boundary and must
  // not be the same.
  // If the error epsilon is provided and > 0, the endpoints of the conic may
  // be slightly modified (by epsilon) so they are located on the boundary.
  Conic_arc_2 (const Conic & conic,
               const Point & source, const Point & target,
               const double & error_eps = 0) :
    _conic(conic),
    _conic_id(0),
    _source(source),
    _target(target),
    _info(X_MON_UNDEFINED),
    _hyper_data_P(NULL)
  {
    static const NT _zero = 0;
    int             deg;
 
    // Find the degree and make sure the conic is not invalid.
    if (_conic.r() != _zero || _conic.s() != _zero || _conic.t() != _zero)
    {
      // In case one of the coefficients of x^2,y^2 or xy is not zero, the
      // degree is 2.
      deg = 2;
    }
    else if (_conic.u() != _zero || _conic.v() != _zero)
    {
      // In case of a line - the degree is 1.
      deg = 1;
    }
    else
    {
      // Empty conic!
      deg = 0;
    }

    CGAL_precondition(deg > 0);

    _info = _info | deg;

    // If it is allowed to slightly move the endpoints by epsilon, check the
    // endpoints before the precondition.
    if (error_eps > 0)
    {
      Point *p_ps[2];
      int   i;

      p_ps[0] = &_source;
      p_ps[1] = &_target;

      for (i = 0; i < 2; i++)
      {
        if (! conic.has_on_boundary(*(p_ps[i])))
        {
          // If the current endpoint is not located on the conic boundary,
          // try to find a y co-ordinate of a point on the conic boundary
          // that has the same x co-ordinate which is epsilon away for the 
          // original y.
          const NT x = (p_ps[i])->x();
          NT       ys[2];
          int      n_ys;
          int      j;

          n_ys = _conic_get_y_coordinates (x, ys);

          for (j = 0; j < n_ys; j++)
          {
            if (fabs(CGAL::to_double((p_ps[i])->y() - ys[j])) < error_eps)
            {
              *(p_ps[i]) = Point (x, ys[j], Point::User_defined);
              break;
            }
          }
        }
      }
    }

    // Make sure the conic contains the two endpoints on its boundary.
    CGAL_precondition(_conic.has_on_boundary(_source));
    CGAL_precondition(_conic.has_on_boundary(_target));

    // Make sure that the source and the target are not the same.
    CGAL_precondition(! _source.equals(_target));      

    // In case the base conic is a hyperbola, build the hyperbolic data.
    if (deg == 2 && _conic.is_hyperbola())
    {
      _build_hyperbolic_arc_data ();
    }

    // In case of a non-degenerate parabola or a hyperbola, make sure 
    // the arc is not infinite.
    if (deg == 2 && ! _conic.is_ellipse())
    {   
      const NT    _two = 2;
      const Point p_mid ((source.x() + target.x()) / _two,
                         (source.y() + target.y()) / _two);
      Point       ps[2];

      bool  finite_at_x = (this->get_points_at_x(p_mid, ps) > 0);
      bool  finite_at_y = (this->get_points_at_y(p_mid, ps) > 0);

      CGAL_precondition(finite_at_x && finite_at_y);
    }

    // If we reached here, the conic arc is legal: Get a new id for the conic.
    _conic_id = _get_new_conic_id();

    _source = Point (_source.x(), _source.y(), 
                     Point::User_defined,
                     _conic_id);

    _target = Point (_target.x(), _target.y(), 
                     Point::User_defined,
                     _conic_id);

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

#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
    _bbox = bounding_box();         // Compute the bounding box. 
#endif
  }

  // Construct a conic arc which is basically a full conic (an ellipse).
  Conic_arc_2 (const Conic& conic) :
    _conic(conic),
    _info(0),
    _hyper_data_P(NULL)
  {
    static const NT _zero = 0;
    int             deg;
 
    // Find the degree and make sure the conic is not invalid.
    if (_conic.r() != _zero || _conic.s() != _zero || _conic.t() != _zero)
    {
      // In case one of the coefficients of x^2,y^2 or xy is not zero, the
      // degree is 2.
      deg = 2;
    }
    else if (_conic.u() != _zero || _conic.v() != _zero)
    {
      // In case of a line - the degree is 1.
      deg = 1;
    }
    else
    {
      // Empty conic!
      deg = 0;
    }

    // Make sure the conic is a non-degenerate ellipse.
    CGAL_precondition(deg == 2 &&_conic.is_ellipse());

    // Set the information: a full conic, which is obvoiusly not x-monotone.
    _info = DEGREE_2 | FULL_CONIC;

    // Assign one of the vertical tangency points as both the source and
    // the target of the conic arc.
    Point     vpts[2];
    int       n_vpts;

    n_vpts = _conic_vertical_tangency_points (vpts);

    CGAL_assertion(n_vpts > 0);
    CGAL_assertion(_conic.has_on_boundary(vpts[0]));

    _source = vpts[0];
    _target = vpts[0];

    // If we reached here, the conic arc is legal: Get a new id for the conic.
    _conic_id = _get_new_conic_id();

#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
    _bbox = bounding_box();       // Compute the bounding box.
#endif
  }
      
  // Destructor.
  virtual ~Conic_arc_2 ()
  {
    if (_hyper_data_P != NULL)
      delete _hyper_data_P;
    _hyper_data_P = NULL;
  }

  // Assignment operator.
  const Conic_arc_2<NT>& operator= (const Conic_arc_2<NT>& arc)
  {
    if (this == &arc)
      return (*this);

    _conic = arc._conic;
    _conic_id = arc._conic_id;
    _source = arc._source;
    _target = arc._target;
    _info = arc._info;

    // Duplicate the data for hyperbolic arcs.
    if (_hyper_data_P != NULL)
      delete _hyper_data_P;
    _hyper_data_P = NULL;

    if (arc._hyper_data_P != NULL)
      _hyper_data_P = new Hyperbolic_arc_data (*(arc._hyper_data_P));

#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
    _bbox = arc._bbox;            // Copy the bounding box.
#endif

    return (*this);
  }

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

  // Check whether the two arcs have the same base conic.
  bool has_same_base_conic (const Conic_arc_2<NT>& arc) const
  {
    if (_conic_id == arc._conic_id)
      return (true);
    else
      return (_conic == arc._conic);
  }

  // Check whether the arc is a full conic (i.e. a non-degenerate ellipse).
  bool is_full_conic () const
  {
    return ((_info & FULL_CONIC) != 0);
  }

  // Check whether the curve is a sgement.
  bool is_segment () const
  {
    return ((_info & DEGREE_MASK) == 1);
  }

  // Check whether the curve is a vertical segment.
  bool is_vertical_segment () const
  {
    // A vertical segment is contained in the degenerate conic: u*x + w = 0.
    static const NT _zero = 0;

    return ((_info & DEGREE_MASK) == 1 && _conic.v() == _zero);
  }

  // Check whether the curve is a horizontal segment.
  bool is_horizontal_segment () const
  {
    // A vertical segment is contained in the degenerate conic: v*y + w = 0.
    static const NT _zero = 0;

    return ((_info & DEGREE_MASK) == 1 && _conic.u() == _zero);
  }

  // Get a bounding box for the conic arc.
  Bbox_2 bounding_box () const
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
    Point  tps[2];
    int    n_tps;
    int    i;

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

  // Check whether the given point is on the arc.
  bool contains_point (const Point& p) const
  { 
    // Check whether the conic contains the point (x,y).
    if (p.is_generating_conic_id(_conic_id) ||
      (p.is_approximate() && _conic_has_approx_point_on_boundary (p)) ||
        _conic.has_on_boundary(p))
    {
      // If the point is on the conic boundary, it is contained in the arc
      // either if the arc is a full conic, or if it is between the two
      // endpoints of the arc.      
      return (is_full_conic() || _is_between_endpoints(p));
    }
    else
    {
      // If the point is not on the conic boundary, it cannot be on the arc.
      return (false);
    }
  }

  // Calculate the vertical tangency points of the arc (ps should be allocated
  // at the size of 2).
  // The function return the number of vertical tangency points.
  int vertical_tangency_points (Point* vpts) const
  {
    // No vertical tangency points for segments or for x-monotone curves:
    if ((_info & DEGREE_MASK) < 2 || (_info & X_MON_UNDEFINED) == X_MONOTONE)
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
    if ((_info & DEGREE_MASK) < 2)
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
    // If the answer is pre-calculated (and stored in the _info field), just
    // return it:
    int    is_x_mon = _info & X_MON_UNDEFINED;

    if (is_x_mon == 0)
      return (false);
    else if (is_x_mon == X_MONOTONE)
      return (true);

    // Check the number of vertical tangency points.
    Point   vpts[2];

    return (vertical_tangency_points(vpts) == 0);
  }

  // Find all points on the arc with a given x-coordinate: ps should be
  // allocated to the size of 2.
  // The function return the number of points found.
  int get_points_at_x (const Point& p,
                       Point *ps) const
  {
    // Make sure the conic is not a vertical segment.
    CGAL_precondition(!is_vertical_segment());

    // Get the y co-ordinates of the points on the conic.
    NT    ys[2];
    int   n;

    n = _conic_get_y_coordinates (p.x(), ys);
    
    // Find all the points that are contained in the arc.
    int   m = 0;
    
    for (int i = 0; i < n; i++)
    {
      ps[m] = Point (p.x(), ys[i], 
                     p.is_approximate() ? Point::Ray_shooting_approx :
                     Point::Ray_shooting_exact,
                     _conic_id);

      if (is_full_conic() || _is_between_endpoints(ps[m]))
        m++;
    }

    // Return the number of points on the arc.
    return (m);
  }
  
  // Find all points on the arc with a given y-coordinate: ps should be
  // allocated to the size of 2.
  // The function return the number of points found.
  int get_points_at_y (const Point& p,
                       Point *ps) const
  {
    // Make sure the conic is not a horizontal segment.
    CGAL_precondition(!is_horizontal_segment());

    // Get the x co-ordinates of the points on the conic.
    NT    xs[2];
    int   n;

    n = _conic_get_x_coordinates (p.y(), xs);
    
    // Find all the points that are contained in the arc.
    int   m = 0;
    
    for (int i = 0; i < n; i++)
    {
      ps[m] = Point (xs[i], p.y(),
                     p.is_approximate() ? Point::Ray_shooting_approx :
                     Point::Ray_shooting_exact,
                     _conic_id);

      if (is_full_conic() || _is_between_endpoints(ps[m]))
        m++;
    }

    // Return the number of points on the arc.
    return (m);
  }
  
  // Reflect the curve in the y axis.
  Conic_arc_2 reflect_in_y () const
  {
    // Reflect the base conic in y:
    Conic ref_conic (_conic.r(), _conic.s(), - _conic.t(),
                     - _conic.u(), _conic.v(), _conic.w());

    // Create the reflected curve:
    Conic_arc_2 ref_arc;

    ref_arc._conic = ref_conic;    
    ref_arc._conic_id = _conic_id ^ REFLECT_IN_Y; 
    ref_arc._source = _source.reflect_in_y();
    ref_arc._target = _target.reflect_in_y();
    ref_arc._info = _info;         // These properties do not change.

    if (_hyper_data_P != NULL)
    {
      ref_arc._hyper_data_P = new Hyperbolic_arc_data (*_hyper_data_P);
      ref_arc._hyper_data_P->a = - _hyper_data_P->a;
    }

#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
    ref_arc._bbox = ref_arc.bounding_box();  // Compute the bounding box.
#endif

    return (ref_arc);
  }

  // Reflect the curve in the x and y axes.
  Conic_arc_2 reflect_in_x_and_y () const
  {
    // Reflect the base conic in x and y:
    Conic ref_conic (_conic.r(), _conic.s(), _conic.t(),
                     - _conic.u(), - _conic.v(), _conic.w());

    // Create the reflected curve:
    Conic_arc_2 ref_arc;

    ref_arc._conic = ref_conic;    
    ref_arc._conic_id = _conic_id ^ (REFLECT_IN_X | REFLECT_IN_Y);
    ref_arc._source = _source.reflect_in_x_and_y();
    ref_arc._target = _target.reflect_in_x_and_y();
    ref_arc._info = _info;         // These properties do not change.

    if (_hyper_data_P != NULL)
    {
      ref_arc._hyper_data_P = new Hyperbolic_arc_data (*_hyper_data_P);
      ref_arc._hyper_data_P->a = - _hyper_data_P->a;
      ref_arc._hyper_data_P->b = - _hyper_data_P->b;
    }

#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
    ref_arc._bbox = ref_arc.bounding_box();  // Compute the bounding box.
#endif

    return (ref_arc);
  }

  // Get the i'th order derivative by x of the conic at the point p=(x,y).
  // Note that i should be either 1 (first order) or 2 (second order).
  void derive_by_x_at (const Point& p, const int& i,
                       NT& slope_numer, NT& slope_denom) const
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
    const NT _two = 2;
    const NT sl_numer = _two*_conic.r()*p.x() + _conic.t()*p.y() + _conic.u();
    const NT sl_denom = _two*_conic.s()*p.y() + _conic.t()*p.x() + _conic.v();

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
    const NT sl2_numer = _conic.s() * sl_numer * sl_numer -
                         _conic.t() * sl_numer * sl_denom +
                         _conic.r() * sl_denom * sl_denom;
    const NT sl2_denom = sl_denom * sl_denom * sl_denom;

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

  // Get the i'th order derivative by y of the conic at the point p=(x,y).
  // Note that i should be either 1 (first order) or 2 (second order).
  void derive_by_y_at (const Point& p, const int& i,
                       NT& slope_numer, NT& slope_denom) const
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
    const NT _two = 2;
    const NT sl_numer = _two*_conic.s()*p.y() + _conic.t()*p.x() + _conic.v();
    const NT sl_denom = _two*_conic.r()*p.x() + _conic.t()*p.y() + _conic.u();

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
    const NT sl2_numer = _conic.r() * sl_numer * sl_numer -
                         _conic.t() * sl_numer * sl_denom +
                         _conic.s() * sl_denom * sl_denom;
    const NT sl2_denom = sl_denom * sl_denom * sl_denom;

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

  // Calculate the intersection points between the arc and the given arc.
  // ps must be allocated at the size of 4.
  // The function returns the number of the actual intersection points.
  int intersections_with (const Conic_arc_2<NT> & arc, Point * ps
#ifdef CGAL_CONIC_ARC_USE_CACHING
                          ,std::list<Intersections> *inter_list_P = NULL
#endif
                          ) const
  {
#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
    // Perform quick rejections, if possible.
    if (! do_overlap(_bbox, arc._bbox))
      return (false);
#endif

    // First make sure that (this->degree) is >= than (arc.degree).
    if ((arc._info & DEGREE_MASK) == DEGREE_2 &&
        (_info & DEGREE_MASK) == DEGREE_1)
    {
      return (arc.intersections_with (*this, ps));
    }

    static const NT _zero = 0;

    // The two conics must not be the same.
    CGAL_precondition(_conic != arc._conic);

    // Deal with vertical segments.
    if (arc.is_vertical_segment())
    {
      if (is_vertical_segment())
      {
        // Two vertical segments intersect only if they overlap.
        return (0);
      }
      
      // Find all points on our arc that have the same x co-ordinate as
      // the other vertical segment.
      int         n_ys;
      Point       xps[2];
      int         j;
      int         n = 0;

      n_ys = get_points_at_x (arc._source, xps);
      
      for (j = 0; j < n_ys; j++)
      {
        // Store this point only if it is contained on the other arc.
        if (arc.contains_point(xps[j]))
        {
          // Return an exact point:
          ps[n] = Point (xps[j].x(), xps[j].y(),
                         Point::Intersection_exact,
                         _conic_id, arc._conic_id);
          n++;
        }
      }
      
      return (n);
    }
    else if (is_vertical_segment())
    {
      // Find all points on the other arc that have the same x co-ordinate as
      // our vertical segment.
      int         n_ys;
      Point       xps[2];
      int         j;
      int         n = 0;

      n_ys = arc.get_points_at_x (_source, xps);
      
      for (j = 0; j < n_ys; j++)
      {
        // Store this point only if it is contained on the other arc.
        if (contains_point(xps[j]))
        {
          // Return an exact point:
          ps[n] = Point (xps[j].x(), xps[j].y(),
                         Point::Intersection_exact,
                         _conic_id, arc._conic_id);
          n++;
        }
      }
      
      return (n);
    }

    // Find all intersection points between the two base conic curves.
    Point     ipts[4];             // The intersection points.
    int       n_points = 0;        // Their number.
    bool      calc_points = true;

#ifdef CGAL_CONIC_ARC_USE_CACHING
    Intersections inter;
    int           k;

    if (inter_list_P != NULL && (_info & DEGREE_MASK) != DEGREE_1)
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
      //bool b_print = true; // cout
      //cout << (*this) << endl << arc << endl;

      // Find all potential x co-ordinates and y co-ordinates of the
      // intersection points.
      NT     xs[4];
      int    n_xs;             // Total number of x co-ordinates.
      int    n_approx_xs;      // Number of approximate x co-ordinates.
      NT     ys[4];
      int    n_ys;             // Total number of y co-ordinates.
      int    n_approx_ys;      // Number of approximate y co-ordinates.
      bool   x_approx, y_approx;
      int    i, j;

      if (_conic.s() == _zero && arc._conic.s() != _zero)
      {
        n_xs = arc._x_coordinates_of_intersections_with (*this,
                                                         xs,
                                                         n_approx_xs);
      }
      else
      {
        n_xs = _x_coordinates_of_intersections_with (arc,
                                                     xs,
                                                     n_approx_xs);
      }

      if (_conic.r() == _zero && arc._conic.r() != _zero)
      {
        n_ys = arc._y_coordinates_of_intersections_with (*this,
                                                         ys, 
                                                         n_approx_ys);
      }
      else
      {
        n_ys = _y_coordinates_of_intersections_with (arc,
                                                     ys, 
                                                     n_approx_ys);
      }
    
      for (i = 0; i < n_xs && n_points < 4; i++)
      {
        x_approx = i >= (n_xs - n_approx_xs);

        for (j = 0; j < n_ys && n_points < 4; j++)
        {
          y_approx = j >= (n_ys - n_approx_ys);

          //if (b_print)
          //  cout << "(" << CGAL::to_double(xs[i]) << ","
          //         << CGAL::to_double(ys[j]) << ") : ";

          if (x_approx || y_approx)
          {
            // At least one co-ordinate is approximate:
            // Create an approximate intersection point and check whether
            // it really lies on both base conics.
            ipts[n_points] = Point (xs[i], ys[j],
                                    Point::Intersection_approx,
                                    _conic_id, arc._conic_id);

            if (_conic_has_approx_point_on_boundary(ipts[n_points]) &&
                arc._conic_has_approx_point_on_boundary(ipts[n_points]))
            {
              //if (b_print) cout << "Approximated, OK !!!" << endl;
              n_points++;
            }
            //else
            //  if (b_print) cout << "Approximated, no." << endl;
          }
          else
          {
            // Both co-ordinates are exact:
            // Create an exact intersection point and check whether
            // it really lies on both base conics.
            ipts[n_points] = Point (xs[i], ys[j],
                                    Point::Intersection_exact,
                                    _conic_id, arc._conic_id);

            if (_conic.has_on_boundary(ipts[n_points]) &&
                arc._conic.has_on_boundary(ipts[n_points]))
            {
              //if (b_print) cout << "Exact, OK !!!" << endl;
              n_points++;
            }
            //else
            //  if (b_print) cout << "Exact, no." << endl;
          }
        }
      }

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
      if (ipts[i].is_approximate())
      {
        // Check for an approximate point.
        if ((is_full_conic() ||
             _is_approximate_between_endpoints(ipts[i])) &&
            (arc.is_full_conic() ||
             arc._is_approximate_between_endpoints(ipts[i])))
        {
          ps[n] = ipts[i];
          n++;
        }
      }
      else
      {
        // Check for an exact point.
        if (contains_point(ipts[i]) &&
            arc.contains_point(ipts[i]))
        {
          ps[n] = ipts[i];
          n++;
        }
      }
    }

    return (n);
  }

  // Check whether the two arcs overlap.
  // The function computes the number of overlapping arcs (2 at most), and
  // return their number (0 means there is no overlap).
  int overlaps (const Conic_arc_2<NT>& arc,
                Conic_arc_2<NT>* ovlp_arcs) const
  {
    // Two arcs can overlap only if their base conics are identical.
    if (_conic != arc._conic)
      return (0);

    // If the two arcs are completely equal, return one of them as the
    // overlapping arc.
    int       orient1 = _conic.orientation();
    int       orient2 = arc._conic.orientation();
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
    const Point *arc_sourceP;
    const Point *arc_targetP;

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
          ovlp_arcs[0] = Conic_arc_2<NT>(*this,_source, *arc_targetP, false);
          ovlp_arcs[1] = Conic_arc_2<NT>(*this, *arc_sourceP, _target, false);
          //ovlp_arcs[0] = Conic_arc_2<NT>(_conic, _source, *arc_targetP);
          //ovlp_arcs[1] = Conic_arc_2<NT>(_conic, *arc_sourceP, _target);
          return (2);
        }

        // Case 1 - *this:     +----------->     
        //            arc:       +=====>
        ovlp_arcs[0] = 
          Conic_arc_2<NT>(*this, *arc_sourceP,*arc_targetP, false);
        //ovlp_arcs[0] = Conic_arc_2<NT>(_conic, *arc_sourceP,*arc_targetP);
        return (1);
      }
      else
      {
        // Case 2 - *this:     +----------->     
        //            arc:               +=====>
        ovlp_arcs[0] = Conic_arc_2<NT>(*this, *arc_sourceP, _target, false);
        //ovlp_arcs[0] = Conic_arc_2<NT>(_conic, *arc_sourceP, _target);
        return (1);
      }
    }
    else if (_is_strictly_between_endpoints(*arc_targetP))
    {
      // Case 3 - *this:     +----------->     
      //            arc:   +=====>
      ovlp_arcs[0] = Conic_arc_2<NT>(*this, _source, *arc_targetP, false);
      //ovlp_arcs[0] = Conic_arc_2<NT>(_conic, _source, *arc_targetP);
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
        
  // Check whether the arc is facing up (LARGER is then returned),
  // facing down (SMALLER is then returned).
  // At any other case the function returns EQUAL.
  Comparison_result facing () const
  {
    if ((_info & FACING_MASK) == 0)
      return (EQUAL);
    else if ((_info & FACING_UP) != 0)
      return (LARGER);
    else
      return (SMALLER);
  }

 private:

  // Build the data for hyperbolic arc, contaning the characterization of the
  // hyperbolic branch the arc is placed on.
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
    const NT        r = _conic.orientation() * _conic.r();
    const NT        s = _conic.orientation() * _conic.s();
    const NT        t = _conic.orientation() * _conic.t();
    const NT        cos_2phi = (r - s) / CGAL::sqrt((r-s)*(r-s) + t*t);
    static const NT _zero(0);
    static const NT _half(0.5);
    static const NT _one(1);
    static const NT _two(2);
    NT              sin_phi;
    NT              cos_phi;

    // Calculate sin(phi) and cos(phi) according to the half-edge formulae:
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
    
    // Calculate the centre (x0, y0) of the conic, given by the formulae:
    //
    //        t*v - 2*s*u                t*u - 2*r*v
    //  x0 = -------------   ,     y0 = -------------
    //        4*r*s - t^2                4*r*s - t^2
    //
    // The denominator (4*r*s - t^2) must be negative for hyperbolas.
    const NT    u = _conic.orientation() * _conic.u();
    const NT    v = _conic.orientation() * _conic.v();
    const NT    det = NT(4)*r*s - t*t;
    NT          x0, y0;

    CGAL_assertion (det < _zero);
    
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

  // Find on which branch of the hyperbola is the given point located.
  // The point is assumed to be on the hyperbola.
  int _hyperbolic_arc_side (const Point& p) const
  {
    if (_hyper_data_P == NULL)
      return (0);

    NT       val;

    val = _hyper_data_P->a*p.x() + _hyper_data_P->b*p.y() + _hyper_data_P->c;
    return ((val > 0) ? 1 : -1);
  }
 
  // Check whether the given point is between the source and the target.
  // The point is assumed to be on the conic's boundary.
  bool _is_between_endpoints (const Point& p) const
  {
    if (p.equals(_source) || p.equals(_target))
      return (true);
    else
      return (_is_strictly_between_endpoints(p));
  }

  bool _is_approximate_between_endpoints (const Point& p) const
  {
    if ((eps_compare<APNT> (TO_APNT(_source.x()), TO_APNT(p.x())) == EQUAL &&
         eps_compare<APNT> (TO_APNT(_source.y()), TO_APNT(p.y())) == EQUAL) ||
        (eps_compare<APNT> (TO_APNT(_target.x()), TO_APNT(p.x())) == EQUAL &&
         eps_compare<APNT> (TO_APNT(_target.y()), TO_APNT(p.y())) == EQUAL))
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

    // In case of a hyperbolic arc, make sure the point is located on the
    // same branch as the arc.
    if (_hyper_data_P != NULL)
    {
      if (_hyperbolic_arc_side(p) != _hyper_data_P->side)
        return (false);
    }

    // Act according to the conic degree.
    if ((_info & DEGREE_MASK) == DEGREE_1)
    {
      if (is_vertical_segment())
      {
        // In case of a vertical segment - just check whether the y co-ordinate
        // of p is between those of the source's and of the target's.
        Comparison_result r1 = compare_y (p, _source);
        Comparison_result r2 = compare_y (p, _target);

        return ((r1 == SMALLER && r2 == LARGER) ||
                (r1 == LARGER && r2 == SMALLER));
      }
      else
      {
        // Otherwise, since the segment is x-monotone, just check whether the
        // x co-ordinate of p is between those of the source's and of the 
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
      if (_conic.orientation() == 1)
        return (left_turn<R>(_source, p, _target));
      else
        return (right_turn<R>(_source, p, _target));
    }
  }

  // Find the y-coordinates of the conic at a given x-coordinate.
  int _conic_get_y_coordinates (const NT& x,
                                NT *ys) const
  {
    // Solve the quadratic equation for a given x and find the y values:
    //  s*y^2 + (t*x + v)*y + (r*x^2 + u*x + w) = 0
    return (solve_quadratic_eq (_conic.s(),
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
    return (solve_quadratic_eq (_conic.r(),
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
    static const NT _four = 4;

    if ((_info & DEGREE_MASK) == DEGREE_1 || _conic.s() == _zero)
      return (0);

    // We are interested in the x co-ordinates were the quadratic equation:
    //  s*y^2 + (t*x + v)*y + (r*x^2 + u*x + w) = 0
    // has a single solution (obviously if s = 0, there are no such points).
    // We therefore demand that the discriminant of this equation is zero:
    //  (t*x + v)^2 - 4*s*(r*x^2 + u*x + w) = 0
    const NT r = _conic.r();
    const NT s = _conic.s();
    const NT t = _conic.t();
    const NT u = _conic.u();
    const NT v = _conic.v();
    const NT w = _conic.w();
    int      n;
    NT       xs[2];
    NT       y;

    n = solve_quadratic_eq (t*t - _four*r*s,
                            _two*t*v - _four*s*u,
                            v*v - _four*s*w,
                            xs);

    for (int i = 0; i < n; i++)
    {
      // Having computed x, y is the simgle solution to the quadratic equation
      // above, and since its discriminant is 0, y is simply given by:
      y = -(t*xs[i] + v) / (_two*s);

      ps[i] = Point (xs[i], y,
                     Point::Tangency,
                     _conic_id);
    }
      
    return (n);
  }

  // Find the horizontal tangency points of the conic.
  int _conic_horizontal_tangency_points (Point* ps) const
  {
    // In case the base conic is of degree 1 (and not 2), the arc has no
    // vertical tangency points.
    static const NT _zero = 0;
    static const NT _two = 2;
    static const NT _four = 4;

    if ((_info & DEGREE_MASK) == DEGREE_1 || _conic.r() == _zero)
      return (0);

    // We are interested in the y co-ordinates were the quadratic equation:
    //  r*x^2 + (t*y + u)*x + (s*y^2 + v*y + w) = 0
    // has a single solution (obviously if r = 0, there are no such points).
    // We therefore demand that the discriminant of this equation is zero:
    //  (t*y + u)^2 - 4*r*(s*y^2 + v*y + w) = 0
    const NT r = _conic.r();
    const NT s = _conic.s();
    const NT t = _conic.t();
    const NT u = _conic.u();
    const NT v = _conic.v();
    const NT w = _conic.w();
    int      n;
    NT       ys[2];
    NT       x;

    n = solve_quadratic_eq (t*t - _four*r*s,
                            _two*t*u - _four*r*v,
                            u*u - _four*r*w,
                            ys);

    for (int i = 0; i < n; i++)
    {
      // Having computed y, x is the simgle solution to the quadratic equation
      // above, and since its discriminant is 0, x is simply given by:
      x = -(t*ys[i] + u) / (_two*r);

      ps[i] = Point (x, ys[i],
                     Point::Tangency,
                     _conic_id);
    }
      
    return (n);
  }
  
  // Check whether the base conic contain the given approximate point on its
  // boundary.
  bool _conic_has_approx_point_on_boundary (const Point& p) const
  {
    APNT r = TO_APNT(_conic.r());
    APNT s = TO_APNT(_conic.s());
    APNT t = TO_APNT(_conic.t());
    APNT u = TO_APNT(_conic.u());
    APNT v = TO_APNT(_conic.v());
    APNT w = TO_APNT(_conic.w());
    APNT x = TO_APNT(p.x());
    APNT y = TO_APNT(p.y());

    APNT value = (r*x + u)*x + (s*y + t*x + v)*y + w;

    return (APNT_CZERO(value));
  }

  // Set the facing information for the arc.
  void _set_facing ()
  {
    // Check whether the arc (which is x-monotone of degree 2) lies above or 
    // below the segement that contects its two end-points (x1,y1) and (x2,y2).
    // To do that, we find the y co-ordinate of a point on the arc whose x
    // co-ordinate is (x1+x2)/2 and compare it to (y1+y2)/2.
    static const NT _two = 2;
    const NT        x_mid = (_source.x() + _target.x()) / _two;
    const NT        y_mid = (_source.y() + _target.y()) / _two;
    Point           p_mid (x_mid, y_mid);
    Point           ps[2];
    int             n_ps;

    n_ps = get_points_at_x (p_mid, ps);
    CGAL_assertion (n_ps == 1);

    Comparison_result res = CGAL::compare (ps[0].y(), y_mid);

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

  // Calculate all x co-ordinates of intersection points between the two
  // base curves of (*this) and the given arc.
  int _x_coordinates_of_intersections_with (const Conic_arc_2<NT>& arc,
                                            NT* xs, int& n_approx) const
  {
    static const NT _zero = 0;
    int             n_roots;         // The number of distinct x values.

    if ((_info & DEGREE_MASK) == DEGREE_1 && 
        (arc._info & DEGREE_MASK) == DEGREE_1)
    {
      // The two conics are: u*x + v*y + w = 0
      //                and: u'*x + v'*y + w' = 0
      // There's a single solution for x, which is:
      const NT denom = _conic.v()*arc._conic.w() - _conic.w()*arc._conic.v();
      const NT numer = _conic.u()*arc._conic.v() - _conic.v()*arc._conic.u();

      if (numer == _zero)
      {
        n_roots = 0;
      }
      else
      {
        xs[0] = denom / numer;
        n_roots = 1;
        n_approx = 0;
      }
    }
    else if ((arc._info & DEGREE_MASK) == DEGREE_1)
    {
      // The two conics are: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
      //                and: a*x + b*y + c = 0
      // There are therefore 2 possible x-values, the solutions for:
      const NT a = arc._conic.u();
      const NT b = arc._conic.v();
      const NT c = arc._conic.w();
      const NT r = _conic.r();
      const NT s = _conic.s();
      const NT t = _conic.t();
      const NT u = _conic.u();
      const NT v = _conic.v();
      const NT w = _conic.w();

      n_roots = solve_quadratic_eq (a*a*s + b*b*r - a*b*t,
                                    2*a*c*s + b*b*u - a*b*v - b*c*t,
                                    c*c*s + b*b*w - b*c*v,
                                    xs);      
      n_approx = 0;
    }
    else if (_conic.s() == _zero && _conic.t() == _zero &&
             arc._conic.s() == _zero && arc._conic.t() == _zero)
    {
      // Special treatment for canonic parabolas whose axes are parallel
      // to the y axis.
      // The two curves are: r*x^2 + u*x + v*y + w = 0
      //                and: r'*x^2 + u'*x + v'*y + w' = 0
      // There are therefore 2 possible x-values, the solutions for:
      NT       c[3];

      c[2] = _conic.r()*arc._conic.v() - _conic.v()*arc._conic.r();
      c[1] = _conic.u()*arc._conic.v() - _conic.v()*arc._conic.u();
      c[0] = _conic.w()*arc._conic.v() - _conic.v()*arc._conic.w();

      n_roots = solve_quadratic_eq (c[2], c[1], c[0],
                                    xs);      
      n_approx = 0;
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
      const NT r = _conic.r();
      const NT s = _conic.s();
      const NT t = _conic.t();
      const NT u = _conic.u();
      const NT v = _conic.v();
      const NT w = _conic.w();
      NT    A, B, C, D, E;

      if (arc._conic.s() != _zero)
      {
        const NT s_tag = arc._conic.s();

        A = s_tag*r - s*arc._conic.r();
        B = s_tag*u - s*arc._conic.u();
        C = s_tag*w - s*arc._conic.w();
        D = s_tag*t - s*arc._conic.t();
        E = s_tag*v - s*arc._conic.v();
      }
      else
      {
        A = arc._conic.r();
        B = arc._conic.u();
        C = arc._conic.w();
        D = arc._conic.t();
        E = arc._conic.v();
      }

      // Now we have:
      //
      //         A*x^2 + B*x + C
      //  y = - -----------------
      //              D*x + E
      //
      // Applying this two our conic's equation yields a quartic equation:
      //
      //  c[4]*x^4 + c[3]*x^3 + c[2]*x^2 + c[1]*x + c[0] = 0
      const NT _two = 2;
      NT       c[5];

      if (t == _zero && arc._conic.t() == _zero)
      {
        c[4] = s*A*A;
        c[3] = _two*s*A*B;
        c[2] = r*E*E + _two*s*A*C + s*B*B - v*A*E;
        c[1] = u*E*E + _two*s*B*C - v*B*E;
        c[0] = w*E*E + s*C*C - v*C*E;
      }
      else
      {
        const NT F = t*E + v*D;
        
        c[4] = r*D*D + s*A*A - t*A*D;
        c[3] = _two*r*D*E + u*D*D + _two*s*A*B - t*B*D - F*A;
        c[2] = r*E*E + _two*u*D*E + w*D*D + 
               _two*s*A*C + s*B*B - t*C*D - F*B - v*A*E;
        c[1] = u*E*E + _two*w*D*E + _two*s*B*C - F*C - v*B*E;
        c[0] = w*E*E + s*C*C - v*C*E;
      }

      n_roots = solve_quartic_eq (c[4], c[3], c[2], c[1], c[0],
                                  xs, n_approx);
    }

    return (n_roots);
  }

  // Calculate all y co-ordinates of intersection points between the two
  // base curves of (*this) and the given arc.
  int _y_coordinates_of_intersections_with (const Conic_arc_2<NT>& arc,
                                            NT* ys, int& n_approx) const
  {
    static const NT _zero = 0;
    int             n_roots;         // The number of distinct y values.

    if ((_info & DEGREE_MASK) == DEGREE_1 && 
        (arc._info & DEGREE_MASK) == DEGREE_1)
    {
      // The two conics are: u*x + v*y + w = 0
      //                and: u'*x + v'*y + w' = 0
      // There's a single solution for y, which is:
      const NT denom = _conic.u()*arc._conic.w() - _conic.w()*arc._conic.u();
      const NT numer = _conic.v()*arc._conic.u() - _conic.u()*arc._conic.v();

      if (numer == _zero)
      {
        n_roots = 0;
      }
      else
      {
        ys[0] = denom / numer;
        n_roots = 1;
        n_approx = 0;
      }
    }
    else if ((arc._info & DEGREE_MASK) == DEGREE_1)
    {
      // The two conics are: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
      //                and: a*x + b*y + c = 0
      // There are therefore 2 possible y-values, the solutions for:
      const NT a = arc._conic.u();
      const NT b = arc._conic.v();
      const NT c = arc._conic.w();
      const NT r = _conic.r();
      const NT s = _conic.s();
      const NT t = _conic.t();
      const NT u = _conic.u();
      const NT v = _conic.v();
      const NT w = _conic.w();

      n_roots = solve_quadratic_eq (b*b*r + a*a*s - a*b*t,
                                    2*b*c*r + a*a*v - a*b*u - a*c*t,
                                    c*c*r + a*a*w - a*c*u,
                                    ys);      
      n_approx = 0;
    }
    else if (_conic.r() == _zero && _conic.t() == _zero &&
             arc._conic.r() == _zero && arc._conic.t() == _zero)
    {
      // Special treatment for canonic parabolas whose axes are parallel
      // to the x axis.
      // The two curves are: s*y^2 + u*x + v*y + w = 0
      //                and: s'*y^2 + u'*x + v'*y + w' = 0
      // There are therefore 2 possible y-values, the solutions for:
      NT       c[3];

      c[2] = _conic.s()*arc._conic.u() - _conic.u()*arc._conic.s();
      c[1] = _conic.v()*arc._conic.u() - _conic.u()*arc._conic.v();
      c[0] = _conic.w()*arc._conic.u() - _conic.u()*arc._conic.w();

      n_roots = solve_quadratic_eq (c[2], c[1], c[0],
                                    ys);      
      n_approx = 0;
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
      const NT r = _conic.r();
      const NT s = _conic.s();
      const NT t = _conic.t();
      const NT u = _conic.u();
      const NT v = _conic.v();
      const NT w = _conic.w();
      NT    A, B, C, D, E;

      if (arc._conic.r() != _zero)
      {
        const NT r_tag = arc._conic.r();

        A = r_tag*s - r*arc._conic.s();
        B = r_tag*v - r*arc._conic.v();
        C = r_tag*w - r*arc._conic.w();
        D = r_tag*t - r*arc._conic.t();
        E = r_tag*u - r*arc._conic.u();
      }
      else
      {
        A = arc._conic.s();
        B = arc._conic.v();
        C = arc._conic.w();
        D = arc._conic.t();
        E = arc._conic.u();
      }

      // Now we have:
      //
      //         A*y^2 + B*y + C
      //  x = - -----------------
      //              D*y + E
      //
      // Applying this two our conic's equation yields a quartic equation:
      //
      //  c[4]*y^4 + c[3]*y^3 + c[2]*y^2 + c[1]*y + c[0] = 0
      const NT _two = 2;
      NT       c[5];

      if (t == _zero && arc._conic.t() == _zero)
      {
        c[4] = r*A*A;
        c[3] = _two*r*A*B;
        c[2] = s*E*E + _two*r*A*C + r*B*B - u*A*E;
        c[1] = v*E*E + _two*r*B*C - u*B*E;
        c[0] = w*E*E + r*C*C - u*C*E;
      }
      else
      {
        const NT F = t*E + u*D;
        
        c[4] = s*D*D + r*A*A - t*A*D;
        c[3] = _two*s*D*E + v*D*D + _two*r*A*B - t*B*D - F*A;
        c[2] = s*E*E + _two*v*D*E + w*D*D + 
               _two*r*A*C + r*B*B - t*C*D - F*B - u*A*E;
        c[1] = v*E*E + _two*w*D*E + _two*r*B*C - F*C - u*B*E;
        c[0] = w*E*E + r*C*C - u*C*E;
      }

      n_roots = solve_quartic_eq (c[4], c[3], c[2], c[1], c[0],
                                  ys, n_approx);
    }

    return (n_roots);
  }

};

#ifndef NO_OSTREAM_INSERT_CONIC_ARC_2
template <class NT>
std::ostream& operator<< (std::ostream& os, const Conic_arc_2<NT>& arc)
{
  Conic_arc_2<NT>::Conic conic = arc.conic();
  Conic_arc_2<NT>::Point source = arc.source(), target = arc.target();

  os << "{" << CGAL::to_double(conic.r()) << "*x^2 + "
     << CGAL::to_double(conic.s()) << "*y^2 + "
     << CGAL::to_double(conic.t()) << "*xy + " 
     << CGAL::to_double(conic.u()) << "*x + "
     << CGAL::to_double(conic.v()) << "*y + "
     << CGAL::to_double(conic.w()) << "} :"
     << "(" << CGAL::to_double(source.x()) << "," 
     << CGAL::to_double(source.y()) << ") -> "
     << "(" << CGAL::to_double(target.x()) << "," 
     << CGAL::to_double(target.y()) << ")";

  return (os);
}
#endif // NO_OSTREAM_INSERT_CONIC_ARC_2

CGAL_END_NAMESPACE

#endif
