#ifndef CGAL_CONIC_ARC_2_H
#define CGAL_CONIC_ARC_2_H

#include <CGAL/basic.h>
#include <list>
#include <vector>

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Conic_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Arrangement_2/Conic_arc_2_eq.h>
#include <CGAL/Arrangement_2/Conic_arc_2_point.h>

#include <fstream>

#define CGAL_CONIC_ARC_USE_BOUNDING_BOX
#define CGAL_CONIC_ARC_USE_CACHING
#define CGAL_CONIC_ARC_USE_FILTER

CGAL_BEGIN_NAMESPACE

enum
{
  REFLECT_IN_X = 1,
  REFLECT_IN_Y = 2,
  REFLECTION_FACTOR = 4
};

// ----------------------------------------------------------------------------
// Representation of a conic curve.
//

static int _conics_count = 0;
template <class _NT> class Arr_conic_traits;

template <class NT>
class Conic_arc_2
{
  friend class Arr_conic_traits<NT>;

 public:

  typedef Cartesian<NT>        Kernel;
  typedef Point_2_ex<Kernel>   Point_2;
  typedef Conic_2<Kernel>      Conic_2;
  typedef Circle_2<Kernel>     Circle_2;
  typedef Segment_2<Kernel>    Segment_2;

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
    FACING_MASK = 32 + 64,
    IS_CIRCLE = 128,
    IS_HYPERBOLA = 256
  };

  Conic_2  _conic;          // The conic that contains the arc.
  int      _conic_id;       // The id of the conic.
  Point_2  _source;         // The source of the arc. 
  Point_2  _target;         // The target of the arc.
  int      _info;           // A bit array with extra information:
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
#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
  Bbox_2   _bbox;           // A bounding box for the arc.  
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

  // In case of a circle, is is convinient to store its center and radius.
  struct Circular_arc_data
  {
    NT       x0;
    NT       y0;
    NT       r;
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
    return (REFLECTION_FACTOR * _conics_count);
  }

  // Private constructor.
  Conic_arc_2 (const Conic_arc_2& arc,
	       const Point_2& source, const Point_2& target,
	       const bool& is_full) :
    _conic(arc._conic),
    _conic_id(arc._conic_id),
    _source(source),
    _target(target)
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

#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
    _bbox = bounding_box();         // Compute the bounding box.
#endif
  }

 public:

  // Default constructor.
  Conic_arc_2 () :
    _conic_id(0),
    _info(X_MON_UNDEFINED)
  {
    _data.hyper_P = NULL;
  }

  // Copy constructor.
  Conic_arc_2 (const Conic_arc_2<NT>& arc) :
    _conic(arc._conic),
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

#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
    _bbox = arc._bbox;              // Copy the bounding box.  
#endif
  }

  // Construct a conic arc which lies on the conic:
  //   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
  // The source and the target must be on the conic boundary and must
  // not be the same.
  Conic_arc_2 (const NT& r, const NT& s, const NT& t,
	       const NT& u, const NT& v, const NT& w,
	       const Point_2& source, const Point_2& target) :
    _conic_id(0),
    _source(source),
    _target(target),
    _info(X_MON_UNDEFINED)
  {
    // Create a conic from the coefficients.
    Conic_2   conic (r, s, t, u, v, w);

    // Make sure the conic contains the two end-points on its boundary.
    CGAL_precondition(conic.has_on_boundary(source));
    CGAL_precondition(conic.has_on_boundary(target));

    // Make sure that the source and the taget are not the same.
    CGAL_precondition(source != target);      

    // Set the arc properties.
    _set (conic, source, target);
  }

  // Construct a conic arc which is the full conic:
  //   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
  // The conic C must be an ellipse (so 4rs - t^2 > 0).
  Conic_arc_2 (const NT& r, const NT& s, const NT& t,
	       const NT& u, const NT& v, const NT& w)
  {
    // Make sure the given curve is an ellipse.
    CGAL_precondition(4*r*s - t*t > 0);
    
    // Create a conic from the coefficients.
    Conic_2   conic (r, s, t, u, v, w);
    
    // Set the arc to be the full conic.
    _set_full (conic);
  }

  // Construct a conic arc from the given line segment.
  Conic_arc_2 (const Segment_2& segment)
  {
    // Use the source and target of the given segment.
    Point_2      source = Point_2(segment.source().x(), segment.source().y());
    Point_2      target = Point_2(segment.target().x(), segment.target().y());
    Conic_2      conic;
    static const NT _zero = 0;
    static const NT _one = 1;

    // Make sure that the source and the taget are not the same.
    CGAL_precondition(source != target);      

    // The supporting conic is r=s=t=0, and u*x + v*y + w = 0 should hold
    // for both the source (x1,y1) and the target (x2, y2).
    if (source.x() == target.x())
    {
      // The supporting conic is a vertical line, of the form x = CONST.
      conic.set (_zero, _zero, _zero,    // r = s = t = 0
		 _one,                   // u = 1
		 _zero,                  // v = 0
		 -source.x());           // w = -CONST
    }
    else
    {
      // The supporting line is A*x + B*y + C = 0, where:
      //
      //  A = y2 - y1,    B = x1 - x2,    C = x2*y1 - x1*y2 
      //
      const NT    A = (target.y() - source.y());
      const NT    B = (source.x() - target.x());
      const NT    C = (target.x()*source.y() - source.x()*target.y());

      // Now we can set:
      conic.set (_zero, _zero, _zero,    // r = s = t = 0
		 A,                      // u = A
		 B,                      // v = B
		 C);                     // w = C
    }

    // Set the arc properties.
    _set (conic, source, target);
  }

  // Set a circular arc that lies on the given circle with the given
  // end-point.
  // Note that the orientation of the input circle is preserved.
  Conic_arc_2 (const Circle_2& circle,
	       const Point_2& source, const Point_2& target)
  {
    // Make sure the circle contains the two endpoints on its boundary.
    CGAL_precondition(circle.has_on_boundary(source));
    CGAL_precondition(circle.has_on_boundary(target));

    // Make sure that the source and the taget are not the same.
    CGAL_precondition(source != target);      

    // Produce the correponding conic: if the circle centre is (x0,y0)
    // and it radius is r, that its equation is:
    //   x^2 + y^2 - 2*x0*x - 2*y0*y + (x0^2 + y0^2 - r^2) = 0
    // Since this equation describes a curve with a negative orientation,
    // we multiply it by -1 if necessary to preserve the original orientation
    // of the input circle.
    static const NT _zero = 0;
    static const NT _one = 1;
    static const NT _minus_one = -1;
    static const NT _two = 2;
    static const NT _minus_two = -2;
    const NT    x0 = circle.center().x();
    const NT    y0 = circle.center().y();
    const NT    r_squared = circle.squared_radius();
    Conic_2     conic;

    if (circle.orientation() == CGAL::COUNTERCLOCKWISE)
    {
      conic.set (_minus_one, _minus_one,      // r = s = -1
		 _zero,                       // t = 0
		 _two*x0,
		 _two*y0,
		 r_squared - x0*x0 - y0*y0);
    }
    else
    {
      conic.set (_one, _one,                  // r = s = 1
		 _zero,                       // t = 0
		 _minus_two*x0,
		 _minus_two*y0,
		 x0*x0 + y0*y0 - r_squared);
    }

    // Prepare the auxiliary data structure.
    Circular_arc_data    *circ_data_P = new Circular_arc_data;

    circ_data_P->x0 = x0;
    circ_data_P->y0 = y0;
    circ_data_P->r = CGAL::sqrt(circle.squared_radius());

    // Set the arc properties.
    _set (conic, source, target, circ_data_P);
  }

  // Construct an arc from the full circle.
  Conic_arc_2 (const Circle_2& circle)
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
    Conic_2     conic (_one, _one,                  // r = s = 1
		       _zero,                       // t = 0
		       _minus_two*x0,
		       _minus_two*y0,
		       x0*x0 + y0*y0 - r_squared);

    // Prepare the auxiliary data structure.
    Circular_arc_data    *circ_data_P = new Circular_arc_data;

    circ_data_P->x0 = x0;
    circ_data_P->y0 = y0;
    circ_data_P->r = CGAL::sqrt(circle.squared_radius());

    // Set the arc to be the full conic.
    _set_full (conic, circ_data_P);
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
  Conic_arc_2 (const NT& r, const NT& s, const NT& t,
	       const NT& u, const NT& v, const NT& w,
	       const Point_2& app_source,
	       const NT& r_1, const NT& s_1, const NT& t_1,
	       const NT& u_1, const NT& v_1, const NT& w_1,
	       const Point_2& app_target,
	       const NT& r_2, const NT& s_2, const NT& t_2,
	       const NT& u_2, const NT& v_2, const NT& w_2) :
    _conic_id(0),
    _info(X_MON_UNDEFINED)
  {
    // Create a conic from the given coefficients.
    Conic_2   conic (r, s, t, u, v, w);
    
    _conic = conic;
    if (r == 0 && s == 0 && t == 0)
      _info = _info | DEGREE_1;
    else
      _info = _info | DEGREE_2;

    // Create fictitious arcs for the source and target computation.
    Conic_2     conic_s (r_1, s_1, t_1, u_1, v_1, w_1);
    Conic_2     conic_t (r_2, s_2, t_2, u_2, v_2, w_2);
    Conic_arc_2 arc_s;
    Conic_arc_2 arc_t;

    arc_s._conic = conic_s;
    if (r_1 == 0 && s_1 == 0 && t_1 == 0)
      arc_s._info = arc_s._info | DEGREE_1;
    else
      arc_s._info = arc_s._info | DEGREE_2;

    arc_t._conic = conic_t;
    if (r_2 == 0 && s_2 == 0 && t_2 == 0)
      arc_t._info = arc_t._info | DEGREE_1;
    else
      arc_t._info = arc_t._info | DEGREE_2;
    
    // Compute the source and the target.
    Point_2     source;
    Point_2     target;
    Conic_arc_2   *arc_P;
    const Point_2 *app_P;
    Point_2       *end_P;
    int         i, j;
    Point_2     ipts[4];         // The intersection points.
    int         n_points = 0;    // Their number.
    NT          xs[4];
    int         x_mults[4];
    int         n_xs;            // Total number of x co-ordinates.
    int         n_approx_xs;     // Number of approximate x co-ordinates.
    NT          ys[4];
    int         y_mults[4];
    int         n_ys;            // Total number of y co-ordinates.
    int         n_approx_ys;     // Number of approximate y co-ordinates.
    int             x_deg[2];    // The generating polynomial for the x values.
    std::vector<NT> x_coeffs[2];
    int             y_deg[2];    // The generating polynomial for the y values.
    std::vector<NT> y_coeffs[2];
    APNT        dist, best_dist = 0;
    int         index;

    for (i = 0; i < 2; i++)
    {
      if (i == 0)
      {
	arc_P = &arc_s;
	app_P = &app_source;
	end_P = &source;
      }
      else
      {
	arc_P = &arc_t;
	app_P = &app_target;
	end_P = &target;
      }

      n_xs = _x_coordinates_of_intersections_with (*arc_P,
						   xs, x_mults,
						   n_approx_xs,
						   x_deg[i],
						   x_coeffs[i]);

      n_ys = _y_coordinates_of_intersections_with (*arc_P,
						   ys, y_mults,
						   n_approx_ys,
						   y_deg[i], 
						   y_coeffs[i]);
    
      n_points = _pair_intersection_points (*arc_P,
					    n_xs,
					    xs, x_mults,
					    n_approx_xs,
					    n_ys,
					    ys, y_mults,
					    n_approx_ys,
					    ipts);

      // Match the best point.
      index = -1;
      for (j = 0; j < n_points; j++)
      {
	dist = TO_APNT((ipts[j].x() - app_P->x())*(ipts[j].x() - app_P->x()) +
		       (ipts[j].y() - app_P->y())*(ipts[j].y() - app_P->y()));
	
	if (index == -1 || dist < best_dist)
	{
	  index = j;
	  best_dist = dist;
	}
      }

      CGAL_assertion(index != -1);
      *end_P = ipts[index];
    }
    
    // Set the arc.
    _set (conic, source, target);
    
    // Make sure the end-point carry all information about their background
    // polynomials.
    _source.attach_polynomials (x_deg[0], x_coeffs[0],
				y_deg[0], y_coeffs[0]);

    _target.attach_polynomials (x_deg[1], x_coeffs[1],
				y_deg[1], y_coeffs[1]);
  }

  // Destructor.
  virtual ~Conic_arc_2 ()
  {
    if ((_info & IS_HYPERBOLA) != 0)
      delete _data.hyper_P;
    else if ((_info & IS_CIRCLE) != 0)
      delete _data.circ_P;
    _data.hyper_P = NULL;
  }

  // Assignment operator.
  const Conic_arc_2<NT>& operator= (const Conic_arc_2<NT>& arc)
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
    _conic = arc._conic;
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

#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
    _bbox = arc._bbox;            // Copy the bounding box.
#endif

    return (*this);
  }

  // Get the arc's base conic.
  const Conic_2& conic () const
  {
    return (_conic);
  }

  // Get the arc's source.
  const Point_2& source () const
  {
    return (_source);
  }

  // Get the arc's target.
  const Point_2& target () const
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
    Point_2  tps[2];
    int      n_tps;
    int      i;

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
  bool contains_point (const Point_2& p) const
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
  int vertical_tangency_points (Point_2* vpts) const
  {
    // No vertical tangency points for segments or for x-monotone curves:
    if ((_info & DEGREE_MASK) < 2 ||
	(_info & X_MON_UNDEFINED) == X_MONOTONE)
      return (0);

    // Calculate the vertical tangency points of the conic.
    Point_2  ps[2];
    int      n;

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
  int horizontal_tangency_points (Point_2* hpts) const
  {
    // No horizontal tangency points for segments:
    if ((_info & DEGREE_MASK) < 2)
      return (0);

    // Calculate the horizontal tangency points of the conic.
    Point_2  ps[2];
    int      n;

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
    Point_2   vpts[2];

    return (vertical_tangency_points(vpts) == 0);
  }

  // Find all points on the arc with a given x-coordinate: ps should be
  // allocated to the size of 2.
  // The function return the number of points found.
  int get_points_at_x (const Point_2& p,
                       Point_2 *ps) const
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
      ps[m] = Point_2 (p.x(), ys[i], 
		       p.is_approximate() ? Point_2::Ray_shooting_approx :
		                            Point_2::Ray_shooting_exact,
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
  int get_points_at_y (const Point_2& p,
                       Point_2 *ps) const
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
      ps[m] = Point_2 (xs[i], p.y(),
		       p.is_approximate() ? Point_2::Ray_shooting_approx :
		                            Point_2::Ray_shooting_exact,
		       _conic_id);

      if (is_full_conic() || _is_between_endpoints(ps[m]))
	m++;
    }

    // Return the number of points on the arc.
    return (m);
  }
  
  // Return a flipped conic arc.
  Conic_arc_2 flip () const
  {
    // Create an identical conic with an opposite orientation.
    Conic_2         opp_conic (-_conic.r(), -_conic.s(), -_conic.t(),
			       -_conic.u(), -_conic.v(), -_conic.w());


    // Create the reflected curve (exchange the source and the target):
    Conic_arc_2     opp_arc;

    opp_arc._conic = opp_conic;    
    opp_arc._conic_id = _conic_id; 
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
 
#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
    opp_arc._bbox = _bbox;         // The bounding box is the same.
#endif

    return (opp_arc);
  }

  // Reflect the curve in the y axis.
  Conic_arc_2 reflect_in_y () const
  {
    // Reflect the base conic in y:
    Conic_2 ref_conic (  _conic.r(),
		         _conic.s(),
		       - _conic.t(),
		       - _conic.u(),
		         _conic.v(),
		         _conic.w());

    // Create the reflected curve:
    Conic_arc_2 ref_arc;

    ref_arc._conic = ref_conic;    
    ref_arc._conic_id = _conic_id ^ REFLECT_IN_Y; 
    ref_arc._source = _source.reflect_in_y();
    ref_arc._target = _target.reflect_in_y();
    ref_arc._info = _info;         // These properties do not change.

    if ((_info & IS_HYPERBOLA) != 0)
    {
      ref_arc._data.hyper_P = new Hyperbolic_arc_data (*_data.hyper_P);
      ref_arc._data.hyper_P->a = - _data.hyper_P->a;
    }
    else if ((_info & IS_CIRCLE) != 0)
    {
      ref_arc._data.circ_P = new Circular_arc_data (*_data.circ_P);
      ref_arc._data.circ_P->x0 = - _data.circ_P->x0;
    }
    else
    {
      ref_arc._data.hyper_P = NULL;
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
    Conic_2 ref_conic (  _conic.r(),
		         _conic.s(),
		         _conic.t(),
		       - _conic.u(),
		       - _conic.v(),
		         _conic.w());

    // Create the reflected curve:
    Conic_arc_2 ref_arc;

    ref_arc._conic = ref_conic;    
    ref_arc._conic_id = _conic_id ^ (REFLECT_IN_X | REFLECT_IN_Y);
    ref_arc._source = _source.reflect_in_x_and_y();
    ref_arc._target = _target.reflect_in_x_and_y();
    ref_arc._info = _info;         // These properties do not change.

    if ((_info & IS_HYPERBOLA) != 0)
    {
      ref_arc._data.hyper_P = new Hyperbolic_arc_data (*_data.hyper_P);
      ref_arc._data.hyper_P->a = - _data.hyper_P->a;
      ref_arc._data.hyper_P->b = - _data.hyper_P->b;
    }
    else if ((_info & IS_CIRCLE) != 0)
    {
      ref_arc._data.circ_P = new Circular_arc_data (*_data.circ_P);
      ref_arc._data.circ_P->x0 = - _data.circ_P->x0;
      ref_arc._data.circ_P->y0 = - _data.circ_P->y0;
    }
    else
    {
      ref_arc._data.hyper_P = NULL;
    }


#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
    ref_arc._bbox = ref_arc.bounding_box();  // Compute the bounding box.
#endif

    return (ref_arc);
  }

  // Get the i'th order derivative by x of the conic at the point p=(x,y).
  // Note that i should be either 1 (first order) or 2 (second order).
  void derive_by_x_at (const Point_2& p, const int& i,
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
  void derive_by_y_at (const Point_2& p, const int& i,
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
  int intersections_with (const Conic_arc_2<NT>& arc,
			  Point_2* ps
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
      return (arc.intersections_with (*this, ps
#ifdef CGAL_CONIC_ARC_USE_CACHING
				      ,inter_list_P
#endif
				      ));
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
      Point_2     xps[2];
      int         j;
      int         n = 0;

      n_ys = get_points_at_x (arc._source, xps);
      
      for (j = 0; j < n_ys; j++)
      {
	// Store this point only if it is contained on the other arc.
	if (arc.contains_point(xps[j]))
	{
	  // Return an exact point:
	  ps[n] = Point_2 (xps[j].x(), xps[j].y(),
			   Point_2::Intersection_exact,
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
      Point_2     xps[2];
      int         j;
      int         n = 0;

      n_ys = arc.get_points_at_x (_source, xps);
      
      for (j = 0; j < n_ys; j++)
      {
	// Store this point only if it is contained on the other arc.
	if (contains_point(xps[j]))
	{
	  // Return an exact point:
	  ps[n] = Point_2 (xps[j].x(), xps[j].y(),
			   Point_2::Intersection_exact,
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
    // REF_TRICK
    int       reflect_this = (_conic_id & (REFLECT_IN_X | REFLECT_IN_Y));
    int       k;

    CGAL_assertion (reflect_this == 
                    (arc._conic_id & (REFLECT_IN_X | REFLECT_IN_Y)));

#ifdef CGAL_CONIC_ARC_USE_CACHING
    Intersections inter;

    if (inter_list_P != NULL &&
	(_info & DEGREE_MASK) != DEGREE_1)
    {
      // REF-TRICK
      int           id1 = _conic_id / REFLECTION_FACTOR;
      int           id2 = arc._conic_id / REFLECTION_FACTOR;
      //int           id1 = _conic_id;
      //int           id2 = arc._conic_id;
    
      inter.id1 = id1 < id2 ? id1 : id2;
      inter.id2 = id1 > id2 ? id1 : id2;
    
      typename std::list<Intersections>::iterator iter;
      for (iter = inter_list_P->begin(); iter != inter_list_P->end(); iter++)
      {
	if ((*iter).id1 == inter.id1 && (*iter).id2 == inter.id2)
	{
	  n_points = (*iter).n_points;
  
	  for (k = 0; k < n_points; k++)
	  {
	    // REF-TRICK
	    if (reflect_this == (REFLECT_IN_X | REFLECT_IN_Y))
	      ipts[k] = (*iter).ps[k].reflect_in_x_and_y();
	    else if (reflect_this == (REFLECT_IN_Y))
	      ipts[k] = (*iter).ps[k].reflect_in_y();
	    else
	      ipts[k] = (*iter).ps[k];
	    //ipts[k] = (*iter).ps[k];
	  }
	  calc_points = false;
	}
      }
    }
#endif // (of ifdef CGAL_CONIC_ARC_USE_CACHING)

    if (calc_points)
    {
      // Find all potential x co-ordinates and y co-ordinates of the
      // intersection points.
      NT     xs[4];
      int    x_mults[4];
      int    n_xs;             // Total number of x co-ordinates.
      int    n_approx_xs;      // Number of approximate x co-ordinates.
      NT     ys[4];
      int    y_mults[4];
      int    n_ys;             // Total number of y co-ordinates.
      int    n_approx_ys;      // Number of approximate y co-ordinates.
      int             x_deg;   // The generating polynomial for the x values.
      std::vector<NT> x_coeffs;
      int             y_deg;   // The generating polynomial for the y values.
      std::vector<NT> y_coeffs;

      if (_conic.s() == _zero && arc._conic.s() != _zero)
      {
	n_xs = arc._x_coordinates_of_intersections_with (*this,
							 xs, x_mults,
							 n_approx_xs,
							 x_deg,
							 x_coeffs);
      }
      else
      {
	n_xs = _x_coordinates_of_intersections_with (arc,
						     xs, x_mults,
						     n_approx_xs,
						     x_deg,
						     x_coeffs);
      }

      if (_conic.r() == _zero && arc._conic.r() != _zero)
      {
	n_ys = arc._y_coordinates_of_intersections_with (*this,
							 ys, y_mults,
							 n_approx_ys,
							 y_deg, y_coeffs);
      }
      else
      {
	n_ys = _y_coordinates_of_intersections_with (arc,
						     ys, y_mults,
						     n_approx_ys,
						     y_deg, y_coeffs);
      }
    
      // Perform the pairing process od the x and y coordinates.
      n_points = _pair_intersection_points (arc,
					    n_xs,
					    xs, x_mults,
					    n_approx_xs,
					    n_ys,
					    ys, y_mults,
					    n_approx_ys,
					    ipts);

      // Attach the generating polynomials for the intersection points.
      for (k = 0; k < n_points; k++)
      {
	ipts[k].attach_polynomials (x_deg, x_coeffs,
				    y_deg, y_coeffs);
      }

#ifdef CGAL_CONIC_ARC_USE_CACHING
      if (inter_list_P != NULL &&
	  (_info & DEGREE_MASK) != DEGREE_1)
      {
	inter.n_points = n_points;
	
	for (k = 0; k < n_points; k++)
	{
	    // REF-TRICK
	    if (reflect_this == (REFLECT_IN_X | REFLECT_IN_Y))
	      inter.ps[k] = ipts[k].reflect_in_x_and_y();
	    else if (reflect_this == (REFLECT_IN_Y))
	      inter.ps[k] = ipts[k].reflect_in_y();
	    else
	      inter.ps[k] = ipts[k];
	    //inter.ps[k] = ipts[k];
	}

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
	  ovlp_arcs[0] = Conic_arc_2<NT>(*this,_source, *arc_targetP, false);
	  ovlp_arcs[1] = Conic_arc_2<NT>(*this, *arc_sourceP, _target, false);
	  //ovlp_arcs[0] = Conic_arc_2<NT>(_conic, _source, *arc_targetP);
	  //ovlp_arcs[1] = Conic_arc_2<NT>(_conic, *arc_sourceP, _target);
	  return (2);
	}

	// Case 1 - *this:     +----------->     
        //            arc:       +=====>
	ovlp_arcs[0] = Conic_arc_2<NT>(*this, *arc_sourceP,*arc_targetP, false);
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

  // Set the properties of a conic arc (for the usage of the constructors).
  // The source and the target are assumed be on the conic boundary.
  void _set (const Conic_2& conic,
	     const Point_2& source, const Point_2& target,
	     Circular_arc_data *circ_data_P = NULL)
  {
    // Set the data members.
    _conic = conic;
    _conic_id = 0;
    _source = source;
    _target = target;
    _info = X_MON_UNDEFINED;
  
    // Find the degree and make sure the conic is not invalid.
    static const NT _zero = 0;
    int             deg;
 
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

    // In case the base conic is a hyperbola, build the hyperbolic data.
    if (deg == 2 && _conic.is_hyperbola())
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
    if (deg == 2 && ! _conic.is_ellipse())
    {
      CGAL_precondition_code(
      const NT      _two = 2;
      const Point_2 p_mid ((source.x() + target.x()) / _two,
			   (source.y() + target.y()) / _two);
      Point_2       ps[2];

      bool  finite_at_x = (this->get_points_at_x(p_mid, ps) > 0);
      bool  finite_at_y = (this->get_points_at_y(p_mid, ps) > 0);
      );
      CGAL_precondition(finite_at_x && finite_at_y);
    }

    // If we reached here, the conic arc is legal: Get a new id for the conic.
    _conic_id = _get_new_conic_id();

    _source = Point_2 (_source.x(), _source.y(), 
		       Point_2::User_defined,
		       _conic_id);

    _target = Point_2 (_target.x(), _target.y(), 
		       Point_2::User_defined,
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

    return;
  }

  // Set an arc which is basically a full conic (an ellipse).
  void _set_full (const Conic_2& conic,
		  Circular_arc_data *circ_data_P = NULL)
  {
    // Set the data members.
    _conic = conic;
    _conic_id = 0;
    _info = 0;

    static const NT _zero = 0;
 
    // Make sure the conic is a non-degenerate ellipse.
    CGAL_precondition(_conic.is_ellipse());

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
    Point_2   vpts[2];
    int       n_vpts;

    if (circ_data_P != NULL)
    {
      vpts[0] = Point_2 (circ_data_P->x0 + circ_data_P->r, 
			 circ_data_P->y0);
    }
    else
    {
      n_vpts = _conic_vertical_tangency_points (vpts);

      CGAL_assertion(n_vpts > 0);
      CGAL_assertion(_conic.has_on_boundary(vpts[0]));
    }

    // If we reached here, the conic arc is legal: Get a new id for the conic.
    _conic_id = _get_new_conic_id();

    _source = Point_2 (vpts[0].x(), vpts[0].y(),
		       Point_2::User_defined,
		       _conic_id);

    _target = Point_2 (vpts[0].x(), vpts[0].y(),
		       Point_2::User_defined,
		       _conic_id);

#ifdef CGAL_CONIC_ARC_USE_BOUNDING_BOX
    _bbox = bounding_box();       // Compute the bounding box.
#endif

    return;
  }

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

  // Find on which branch of the hyperbola is the given point located.
  // The point is assumed to be on the hyperbola.
  int _hyperbolic_arc_side (const Point_2& p) const
  {
    if (_data.hyper_P == NULL)
      return (0);

    NT       val;

    val = _data.hyper_P->a*p.x() + _data.hyper_P->b*p.y() + _data.hyper_P->c;
    return ((val > 0) ? 1 : -1);
  }
 
  // Check whether the given point is between the source and the target.
  // The point is assumed to be on the conic's boundary.
  bool _is_between_endpoints (const Point_2& p) const
  {
    if (p.equals(_source) || p.equals(_target))
      return (true);
    else
      return (_is_strictly_between_endpoints(p));
  }

  // Check whether the given point is strictly between the source and the
  // target (but not any of them).
  // The point is assumed to be on the conic's boundary.
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
      static Kernel                ker;
      typename Kernel::Orientation_2 orient_f = ker.orientation_2_object();
      
      if (_conic.orientation() == 1)
	return (orient_f(_source, p, _target) == LEFTTURN);
      else
	return (orient_f(_source, p, _target) == RIGHTTURN);
    }
  }

  // Find the y-coordinates of the conic at a given x-coordinate.
  int _conic_get_y_coordinates (const NT& x,
                                NT *ys) const
  {
    int    mults[2];

    // Solve the quadratic equation for a given x and find the y values:
    //  s*y^2 + (t*x + v)*y + (r*x^2 + u*x + w) = 0
    return (solve_quadratic_eq (_conic.s(),
				x*_conic.t() + _conic.v(),
				x*(x*_conic.r() + _conic.u()) + _conic.w(),
				ys, mults));
  }

  // Find the x-coordinates of the conic at a given y-coordinate.
  int _conic_get_x_coordinates (const NT& y,
                                NT *xs) const
  {
    int    mults[2];

    // Solve the quadratic equation for a given y and find the x values:
    //  r*x^2 + (t*y + u)*x + (s*y^2 + v*y + w) = 0
    return (solve_quadratic_eq (_conic.r(),
				y*_conic.t() + _conic.u(),
				y*(y*_conic.s() + _conic.v()) + _conic.w(),
				xs, mults));
  }
  
  // Find the vertical tangency points of the conic.
  int _conic_vertical_tangency_points (Point_2* ps) const
  {
    // In case the base conic is of degree 1 (and not 2), the arc has no
    // vertical tangency points.
    static const NT _zero = 0;
    static const NT _two = 2;
    static const NT _four = 4;

    if ((_info & DEGREE_MASK) == DEGREE_1 || _conic.s() == _zero)
      return (0);

    // Special treatment for circles, where the vertical tangency points
    // are simply (x0-r,y0) and (x0+r,y0).
    if ((_info & IS_CIRCLE) != 0)
    {
      ps[0] = Point_2 (_data.circ_P->x0 - _data.circ_P->r, 
		       _data.circ_P->y0,
		       Point_2::Tangency,
		       _conic_id);
      ps[1] = Point_2 (_data.circ_P->x0 + _data.circ_P->r, 
		       _data.circ_P->y0,
		       Point_2::Tangency,
		       _conic_id);

      return (2);
    }

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
    int      n_xs;
    int      x_mults[2];
    NT       ys[2];
    int      n_ys;
    int      y_mults[2];

    n_xs = solve_quadratic_eq (t*t - _four*r*s,
			       _two*t*v - _four*s*u,
			       v*v - _four*s*w,
			       xs, x_mults);

    // Store the generating polynomial for the x-coordinates.
    int             x_deg = 2;
    std::vector<NT> x_coeffs(3);
    int             y_deg;
    std::vector<NT> y_coeffs;

    x_coeffs[2] = t*t - _four*r*s;
    x_coeffs[1] = _two*t*v - _four*s*u;
    x_coeffs[0] = v*v - _four*s*w; 

    // Find the y-coordinates of the vertical tangency points.
    if (t == _zero)
    {
      n_ys = 1;
      ys[0] = -v / (_two*s);
      y_mults[0] = 2;

      y_deg = 0;
   }
    else
    {
      n_ys = solve_quadratic_eq (_four*r*s*s - s*t*t,
				 _four*r*s*v - _two*s*t*u,
				 r*v*v - t*u*v + t*t*w,
				 ys, y_mults);

      // Store the generating polynomial for the y-coordinates.
      y_deg = 2;
      y_coeffs.resize(3);
      y_coeffs[2] = _four*r*s*s - s*t*t;
      y_coeffs[1] = _four*r*s*v - _two*s*t*u;
      y_coeffs[0] = r*v*v - t*u*v + t*t*w; 
    }

    // Pair the x and y coordinates and obtain the vertical tangency points.
    n = 0;
    for (int i = 0; i < n_xs; i++)
    {
      if (n_ys == 1)
      {
	ps[n] = Point_2 (xs[i], ys[0],
			 Point_2::Tangency,
			 _conic_id);
	n++;
      }
      else
      {
	for (int j = 0; j < n_ys; j++)
	{
	  if (ys[j] == -(t*xs[i] + v) / (_two*s))
	  {
	    ps[n] = Point_2 (xs[i], ys[j],
			     Point_2::Tangency,
			     _conic_id);
	    n++;
	    break;
	  }
	}
      }
    }

    // Attach the generating polynomials.
    for (int i = 0; i < n; i++)
    {
      ps[i].attach_polynomials (x_deg, x_coeffs,
				y_deg, y_coeffs);
    }

    return (n);
  }

  // Find the horizontal tangency points of the conic.
  int _conic_horizontal_tangency_points (Point_2* ps) const
  {
    // In case the base conic is of degree 1 (and not 2), the arc has no
    // vertical tangency points.
    static const NT _zero = 0;
    static const NT _two = 2;
    static const NT _four = 4;

    if ((_info & DEGREE_MASK) == DEGREE_1 || _conic.r() == _zero)
      return (0);

    // Special treatment for circles, where the horizontal tangency points
    // are simply (x0,y0-r) and (x0,y0+r).
    if ((_info & IS_CIRCLE) != 0)
    {
      ps[0] = Point_2 (_data.circ_P->x0, 
		       _data.circ_P->y0 - _data.circ_P->r,
		       Point_2::Tangency,
		       _conic_id);
      ps[1] = Point_2 (_data.circ_P->x0, 
		       _data.circ_P->y0 + _data.circ_P->r,
		       Point_2::Tangency,
		       _conic_id);

      return (2);
    }

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
    int      y_mults[2];
    NT       x;

    n = solve_quadratic_eq (t*t - _four*r*s,
			    _two*t*u - _four*r*v,
			    u*u - _four*r*w,
			    ys, y_mults);

    for (int i = 0; i < n; i++)
    {
      // Having computed y, x is the simgle solution to the quadratic equation
      // above, and since its discriminant is 0, x is simply given by:
      x = -(t*ys[i] + u) / (_two*r);

      ps[i] = Point_2 (x, ys[i],
		       Point_2::Tangency,
		       _conic_id);
    }
      
    return (n);
  }
  
  // Check whether the base conic contains the given approximate point on its
  // boundary.
  bool _conic_has_approx_point_on_boundary (const Point_2& p) const
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
    
    // RWRW: Correct this !!!!
    return (APNT_ABS(value) < 0.001);
    //return (eps_compare<APNT>(value, 0) == EQUAL);
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
    Point_2         p_mid (x_mid, y_mid);
    Point_2         ps[2];
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
					    NT* xs, int* x_mults, 
					    int& n_approx,
					    int& x_deg,
					    std::vector<NT>& x_coeffs) const
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
	x_mults[0] = 1;
	n_roots = 1;
	n_approx = 0;
      }

      // Do not set any coefficients' information in this case.
      x_deg = 0;
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

      // Set the generating polynomial.
      x_deg = 2;
      x_coeffs.resize(3);
      x_coeffs[2] = a*a*s + b*b*r - a*b*t;
      x_coeffs[1] = 2*a*c*s + b*b*u - a*b*v - b*c*t;
      x_coeffs[0] = c*c*s + b*b*w - b*c*v;

      // Find the roots.
      n_roots = solve_quadratic_eq (x_coeffs[2],
				    x_coeffs[1],
				    x_coeffs[0],
				    xs, x_mults);      
      n_approx = 0;
    }
    else if ((_info & IS_CIRCLE) != 0 && (arc._info & IS_CIRCLE) != 0)
    {
      // Special treatment for two circles.
      // The two curves are: r*x^2 + r*y^2 + u*x + v*y + w = 0
      //                and: r'*x^2 + r'*y^2 + u'*x + v'*y + w' = 0
      //
      // Thus, r'*C1-r*C2 is a line whose equation is: a*x + b*y + = 0, where:
      const NT r = _conic.r();
      const NT s = _conic.s();
      const NT t = _conic.t();
      const NT u = _conic.u();
      const NT v = _conic.v();
      const NT w = _conic.w();
      const NT a = arc._conic.r()*u - r*arc._conic.u();
      const NT b = arc._conic.r()*v - r*arc._conic.v();
      const NT c = arc._conic.r()*w - r*arc._conic.w();

      if (b == _zero)
      {
	// The line a*x + c = 0 connects both intersection points of the two
	// circles, so the both have an x-coordinate of -c/a.
	if (a == _zero)
	{
	  n_roots = 0;
	  n_approx = 0;
	}
	else
	{
	  xs[0] = -c / a;
	  x_mults[0] = 2;
	  n_roots = 1;
	  n_approx = 0;
	}
	// Do not set any coefficients' information in this case.
	x_deg = 0;
      }
      else
      {
	// The intersection points of the two circles are the same as the
	// intersection points of one of the circles with a*x + b*y + c = 0.
	// Set the generating polynomial.
	x_deg = 2;
	x_coeffs.resize(3);
	x_coeffs[2] = a*a*s + b*b*r - a*b*t;
	x_coeffs[1] = 2*a*c*s + b*b*u - a*b*v - b*c*t;
	x_coeffs[0] = c*c*s + b*b*w - b*c*v;

	// Find the roots.
	n_roots = solve_quadratic_eq (x_coeffs[2],
				      x_coeffs[1],
				      x_coeffs[0],
				      xs, x_mults);      
	n_approx = 0;
      }
    }
    else if (_conic.s() == _zero && _conic.t() == _zero &&
	     arc._conic.s() == _zero && arc._conic.t() == _zero)
    {
      // Special treatment for canonic parabolas whose axes are parallel
      // to the y axis.
      // The two curves are: r*x^2 + u*x + v*y + w = 0
      //                and: r'*x^2 + u'*x + v'*y + w' = 0
      // There are therefore 2 possible x-values, the solutions for:

      // Set the generating polynomial.
      x_deg = 2;
      x_coeffs.resize(3);
      x_coeffs[2] = _conic.r()*arc._conic.v() - _conic.v()*arc._conic.r();
      x_coeffs[1] = _conic.u()*arc._conic.v() - _conic.v()*arc._conic.u();
      x_coeffs[0] = _conic.w()*arc._conic.v() - _conic.v()*arc._conic.w();

      // Solve the equation.
      n_roots = solve_quadratic_eq (x_coeffs[2], 
				    x_coeffs[1], 
				    x_coeffs[0],
				    xs, x_mults);      
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

      if (D == _zero && E == _zero)
      {
	// In this case: A*x^2 + B*x + C = 0, so:
	x_deg = 2;
	x_coeffs.resize(3);

	x_coeffs[2] = A;
	x_coeffs[1] = B;
	x_coeffs[0] = C;

	// Solve the equation.
	n_roots = solve_quadratic_eq (x_coeffs[2], 
				      x_coeffs[1], 
				      x_coeffs[0],
				      xs, x_mults);      
	n_approx = 0;
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
	const NT _two = 2;
      
	// Set the generating polynomial.
	x_deg = 4;
	x_coeffs.resize(5);
	
	if (t == _zero && arc._conic.t() == _zero)
	{
	  x_coeffs[4] = s*A*A;
	  x_coeffs[3] = _two*s*A*B;
	  x_coeffs[2] = r*E*E + _two*s*A*C + s*B*B - v*A*E;
	  x_coeffs[1] = u*E*E + _two*s*B*C - v*B*E;
	  x_coeffs[0] = w*E*E + s*C*C - v*C*E;
	}
	else
	{
	  const NT F = t*E + v*D;
	
	  x_coeffs[4] = r*D*D + s*A*A - t*A*D;
	  x_coeffs[3] = _two*r*D*E + u*D*D + _two*s*A*B - t*B*D - F*A;
	  x_coeffs[2] = r*E*E + _two*u*D*E + w*D*D + 
	                _two*s*A*C + s*B*B - t*C*D - F*B - v*A*E;
	  x_coeffs[1] = u*E*E + _two*w*D*E + _two*s*B*C - F*C - v*B*E;
	  x_coeffs[0] = w*E*E + s*C*C - v*C*E;
	}

	// Try to factor out intersection points that occur at vertical
	// tangency points of one of the curves.
	Point_2   vpts[4];
	int       n_vpts;
	NT        red_coeffs[5], temp_coeffs[5];
	int       red_degree = x_deg, temp_degree;
	int       k, l;
	bool      reduced = false;
	bool      first_time;

	n_vpts = this->_conic_vertical_tangency_points (vpts);
	n_vpts += arc._conic_vertical_tangency_points (vpts + n_vpts);

	for (l = 0; l <= x_deg; l++)
	  red_coeffs[l] = x_coeffs[l];

	n_roots = 0;
	for (k = 0; k < n_vpts; k++)
	{
	  first_time = true;
	  while (factor_root<NT> (red_coeffs, red_degree,
				  vpts[k].x(),
				  temp_coeffs, temp_degree))
	  {
	    if (first_time)
	    {
	      xs[n_roots] = vpts[k].x();
	      x_mults[n_roots] = 1;
	      n_roots++;
	      first_time = false;
	    }
	    else
	    {
	      x_mults[n_roots-1]++;
	    }
	    
	    red_coeffs[red_degree] = _zero;
	    for (l = 0; l <= temp_degree; l++)
	      red_coeffs[l] = temp_coeffs[l];
	    red_degree = temp_degree;
	    reduced = true;
	  }
	}

	if (reduced)
	{
	  n_roots += solve_quartic_eq (red_coeffs[4], 
				       red_coeffs[3], 
				       red_coeffs[2], 
				       red_coeffs[1], 
				       red_coeffs[0],
				       xs + n_roots, x_mults + n_roots,
				       n_approx);
	}
	else
	{
	  // Solve the quartic equation.
	  n_roots = solve_quartic_eq (x_coeffs[4], 
				      x_coeffs[3], 
				      x_coeffs[2], 
				      x_coeffs[1], 
				      x_coeffs[0],
				      xs, x_mults,
				      n_approx);
	}
      }
   }

    return (n_roots);
  }

  // Calculate all y co-ordinates of intersection points between the two
  // base curves of (*this) and the given arc.
  int _y_coordinates_of_intersections_with (const Conic_arc_2<NT>& arc,
					    NT* ys, int* y_mults,
					    int& n_approx,
					    int& y_deg,
					    std::vector<NT>& y_coeffs) const
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
	y_mults[0] = 1;
	n_roots = 1;
	n_approx = 0;
      }

      // Do not store coefficients information in this case.
      y_deg = 0;
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

      // Store the generating polynomial information in this case.
      y_deg = 2;
      y_coeffs.resize(3);
      y_coeffs[2] = b*b*r + a*a*s - a*b*t;
      y_coeffs[1] = 2*b*c*r + a*a*v - a*b*u - a*c*t;
      y_coeffs[0] = c*c*r + a*a*w - a*c*u;

      // Solve the equation.
      n_roots = solve_quadratic_eq (y_coeffs[2],
				    y_coeffs[1],
				    y_coeffs[0],
				    ys, y_mults);      
      n_approx = 0;
    }
    else if ((_info & IS_CIRCLE) != 0 && (arc._info & IS_CIRCLE) != 0)
    {
      // Special treatment for two circles.
      // The two curves are: r*x^2 + r*y^2 + u*x + v*y + w = 0
      //                and: r'*x^2 + r'*y^2 + u'*x + v'*y + w' = 0
      //
      // Thus, r'*C1-r*C2 is a line whose equation is: a*x + b*y + = 0, where:
      const NT r = _conic.r();
      const NT s = _conic.s();
      const NT t = _conic.t();
      const NT u = _conic.u();
      const NT v = _conic.v();
      const NT w = _conic.w();
      const NT a = arc._conic.r()*u - r*arc._conic.u();
      const NT b = arc._conic.r()*v - r*arc._conic.v();
      const NT c = arc._conic.r()*w - r*arc._conic.w();

      if (a == _zero)
      {
	// The line b*y + c = 0 connects both intersection points of the two
	// circles, so the both have a y-coordinate of -c/b.
	if (b == _zero)
	{
	  n_roots = 0;
	  n_approx = 0;
	}
	else
	{
	  ys[0] = -c / b;
	  y_mults[0] = 2;
	  n_roots = 1;
	  n_approx = 0;
	}

	// Do not set any coefficients' information in this case.
	y_deg = 0;
      }
      else
      {
	// The intersection points of the two circles are the same as the
	// intersection points of one of the circles with a*x + b*y + c = 0.
	// Set the generating polynomial.
	y_deg = 2;
	y_coeffs.resize(3);
	y_coeffs[2] = b*b*r + a*a*s - a*b*t;
	y_coeffs[1] = 2*b*c*r + a*a*v - a*b*u - a*c*t;
	y_coeffs[0] = c*c*r + a*a*w - a*c*u;

	// Find the roots.
	n_roots = solve_quadratic_eq (y_coeffs[2],
				      y_coeffs[1],
				      y_coeffs[0],
				      ys, y_mults);      
	n_approx = 0;
      }
    }
    else if (_conic.r() == _zero && _conic.t() == _zero &&
	     arc._conic.r() == _zero && arc._conic.t() == _zero)
    {
      // Special treatment for canonic parabolas whose axes are parallel
      // to the x axis.
      // The two curves are: s*y^2 + u*x + v*y + w = 0
      //                and: s'*y^2 + u'*x + v'*y + w' = 0
      // There are therefore 2 possible y-values, the solutions for:

      // Store the generating polynomial information in this case.
      y_deg = 2;
      y_coeffs.resize(3);
      y_coeffs[2] = _conic.s()*arc._conic.u() - _conic.u()*arc._conic.s();
      y_coeffs[1] = _conic.v()*arc._conic.u() - _conic.u()*arc._conic.v();
      y_coeffs[0] = _conic.w()*arc._conic.u() - _conic.u()*arc._conic.w();

      // Solve the equation:
      n_roots = solve_quadratic_eq (y_coeffs[2], 
				    y_coeffs[1], 
				    y_coeffs[0],
				    ys, y_mults);      
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

      if (D == _zero && E == _zero)
      {
	// In this case: A*y^2 + B*y + C = 0, so:
	y_deg = 2;
	y_coeffs.resize(3);

	y_coeffs[2] = A;
	y_coeffs[1] = B;
	y_coeffs[0] = C;

	// Solve the equation.
	n_roots = solve_quadratic_eq (y_coeffs[2], 
				      y_coeffs[1], 
				      y_coeffs[0],
				      ys, y_mults);      
	n_approx = 0;
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
	const NT _two = 2;

	// Store the generating polynomial information in this case.
	y_deg = 4;
	y_coeffs.resize(5);

	if (t == _zero && arc._conic.t() == _zero)
	{
	  y_coeffs[4] = r*A*A;
	  y_coeffs[3] = _two*r*A*B;
	  y_coeffs[2] = s*E*E + _two*r*A*C + r*B*B - u*A*E;
	  y_coeffs[1] = v*E*E + _two*r*B*C - u*B*E;
	  y_coeffs[0] = w*E*E + r*C*C - u*C*E;
	}
	else
	{
	  const NT F = t*E + u*D;
	
	  y_coeffs[4] = s*D*D + r*A*A - t*A*D;
	  y_coeffs[3] = _two*s*D*E + v*D*D + _two*r*A*B - t*B*D - F*A;
	  y_coeffs[2] = s*E*E + _two*v*D*E + w*D*D + 
	                _two*r*A*C + r*B*B - t*C*D - F*B - u*A*E;
	  y_coeffs[1] = v*E*E + _two*w*D*E + _two*r*B*C - F*C - u*B*E;
	  y_coeffs[0] = w*E*E + r*C*C - u*C*E;
	}

	// Try to factor out intersection points that occur at vertical
	// tangency points of one of the curves.
	Point_2   vpts[4];
	int       n_vpts;
	NT        red_coeffs[5], temp_coeffs[5];
	int       red_degree = y_deg, temp_degree;
	int       k, l;
	bool      reduced = false;
	bool      first_time;

	n_vpts = this->_conic_vertical_tangency_points (vpts);
	n_vpts += arc._conic_vertical_tangency_points (vpts + n_vpts);

	for (l = 0; l <= y_deg; l++)
	  red_coeffs[l] = y_coeffs[l];

	n_roots = 0;
	for (k = 0; k < n_vpts; k++)
	{
	  first_time = true;
	  while (factor_root<NT> (red_coeffs, red_degree,
				  vpts[k].y(),
				  temp_coeffs, temp_degree))
	  {
	    if (first_time)
	    {
	      ys[n_roots] = vpts[k].y();
	      y_mults[n_roots] = 1;
	      n_roots++;
	      first_time = false;
	    }
	    else
	    {
	      y_mults[n_roots-1]++;
	    }
	    
	    red_coeffs[red_degree] = _zero;
	    for (l = 0; l <= temp_degree; l++)
	      red_coeffs[l] = temp_coeffs[l];
	    red_degree = temp_degree;
	    reduced = true;
	  }
	}

	if (reduced)
	{
	  n_roots += solve_quartic_eq (red_coeffs[4], 
				       red_coeffs[3], 
				       red_coeffs[2], 
				       red_coeffs[1], 
				       red_coeffs[0],
				       ys + n_roots, y_mults + n_roots,
				       n_approx);
	}
	else
	{
	  // Solve the equation:
	  n_roots = solve_quartic_eq (y_coeffs[4], 
				      y_coeffs[3], 
				      y_coeffs[2],
				      y_coeffs[1],
				      y_coeffs[0],
				      ys, y_mults,
				      n_approx);
        }
      }
    }

    return (n_roots);
  }

  // Pair the x coordinates and the y coordinates of the intersection point
  // of (*this) and arc and return a vector of intersection points. 
  int _pair_intersection_points (const Conic_arc_2<NT>& arc,
				 const int& n_xs,
				 const NT* xs, const int* x_mults,
				 const int& n_approx_xs,
				 const int& n_ys,
				 const NT* ys, const int* y_mults,
				 const int& n_approx_ys,
				 Point_2* ipts) const
  {
    // Calculate the minimal number of intersection points.
    int          x_total = 0, y_total = 0;
    int          min_points;
    int          i, j;

    for (i = 0; i < n_xs; i++)
      x_total += x_mults[i];

    for (j = 0; j < n_ys; j++)
      y_total += y_mults[j];

    min_points = (x_total < y_total) ? x_total : y_total;

    // Go over all x coordinates.
    const APNT   r1 = TO_APNT(_conic.r());
    const APNT   s1 = TO_APNT(_conic.s());
    const APNT   t1 = TO_APNT(_conic.t());
    const APNT   u1 = TO_APNT(_conic.u());
    const APNT   v1 = TO_APNT(_conic.v());
    const APNT   w1 = TO_APNT(_conic.w());
    const APNT   r2 = TO_APNT(arc._conic.r());
    const APNT   s2 = TO_APNT(arc._conic.s());
    const APNT   t2 = TO_APNT(arc._conic.t());
    const APNT   u2 = TO_APNT(arc._conic.u());
    const APNT   v2 = TO_APNT(arc._conic.v());
    const APNT   w2 = TO_APNT(arc._conic.w());
    APNT         x, y;
    int          n_ipts = 0;
    APNT         results[16], res;
    int          x_indices[16], y_indices[16], temp;
    bool         is_approx;
    int          k, l;

    k = 0;
    for (i = 0; i < n_xs; i++)
    {
      // For each y, calculate: c[j] = |C1(x[i],y[j])| + |C2(x[i],y[j])|.
      x = TO_APNT(xs[i]);
      for (j = 0; j < n_ys; j++)
      {
	y = TO_APNT(ys[j]);
	results[k] = APNT_ABS(r1*x*x + s1*y*y + t1*x*y + u1*x + v1*y + w1) +
	             APNT_ABS(r2*x*x + s2*y*y + t2*x*y + u2*x + v2*y + w2);
	x_indices[k] = i;
	y_indices[k] = j;
	k++;
      }
    }

    // Perform bubble-sort.
    for (k = 0; k < n_xs*n_ys - 1; k++)
    {
      for (l = k + 1; l < n_xs*n_ys; l++)
      {
	if (results[k] > results[l])
	{
	  res = results[k];
	  results[k] = results[l];
	  results[l] = res;

	  temp = x_indices[k];
	  x_indices[k] = x_indices[l];
	  x_indices[l] = temp;
	  temp = y_indices[k];
	  y_indices[k] = y_indices[l];
	  y_indices[l] = temp;
	}
      }
    }

    // Use the x,y pairs for which the c[j] values are the lowest.
    n_ipts = 0;
    k = 0;
    while (k < n_xs*n_ys && n_ipts < 4)
    {
      if (n_ipts >= min_points &&
	  eps_compare<APNT>(results[k]*results[k], 0) != EQUAL)
      {
	break;
      }

      is_approx = (x_indices[k] >= (n_xs - n_approx_xs)) ||
	          (y_indices[k] >= (n_ys - n_approx_ys));

      ipts[n_ipts] = Point_2 (xs[x_indices[k]], ys[y_indices[k]],
			      (is_approx ?
			       Point_2::Intersection_approx :
			       Point_2::Intersection_exact),
			      _conic_id, arc._conic_id);
      k++;
      n_ipts++;
    }
    
    return (n_ipts);
  }
};

#ifndef NO_OSTREAM_INSERT_CONIC_ARC_2
template <class NT>
std::ostream& operator<< (std::ostream& os, const Conic_arc_2<NT>& arc)
{
  typename Conic_arc_2<NT>::Conic_2 conic = arc.conic();
  typename Conic_arc_2<NT>::Point_2 source = arc.source(),
      target = arc.target();

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
