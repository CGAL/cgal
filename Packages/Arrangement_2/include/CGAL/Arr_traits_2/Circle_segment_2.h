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

#ifndef CGAL_CIRCULAR_ARC_2_H
#define CGAL_CIRCULAR_ARC_2_H

/*! \file
 * Header file for the _Circle_segment_2<Kernel> class.
 */

#include <CGAL/Arr_traits_2/One_root_number.h>
#include <CGAL/Bbox_2.h>
#include <ostream>

CGAL_BEGIN_NAMESPACE

/*! \class
 * Representation of a point whose coordinates are one-root numbers.
 */
template <class NumberType_>
class _One_root_point_2
{
public:

  typedef NumberType_                 NT;
  typedef _One_root_number<NT>        CoordNT;

private:

  CoordNT       _x;
  CoordNT       _y;

public:

  /*! Default constructor. */
  _One_root_point_2 () :
    _x (0),
    _y (0)
  {}

  /*! Constructor of a point with rational coefficients. */
  _One_root_point_2 (const NT& x, const NT& y) :
    _x (x),
    _y (y)
  {}

  /*! Constructor of a point with one-root coefficients. */
  _One_root_point_2 (const CoordNT& x, const CoordNT& y) :
    _x (x),
    _y (y)
  {}

  /*! Get the x-coordinate. */
  const CoordNT& x () const
  {
    return (_x);
  }

  /*! Get the y-coordinate. */
  const CoordNT& y () const
  {
    return (_y);
  }

  /*! Check for equality. */
  bool equals (const _One_root_point_2<NT>& p) const
  {
    return (CGAL::compare (_x, p._x) == EQUAL &&
            CGAL::compare (_y, p._y) == EQUAL);
  }
};

/*!
 * Exporter for conic arcs.
 */
template <class NT>
std::ostream& 
operator<< (std::ostream& os, 
            const _One_root_point_2<NT>& p)
{
  os << CGAL::to_double(p.x()) << ' ' << CGAL::to_double(p.y());
  return (os);
}

/*! \class
 * Representation of a circle, a circular arc or a line segment.
 */
template <class Kernel_>
class _Circle_segment_2
{
public:

  typedef Kernel_                                          Kernel;
  typedef typename Kernel::FT                              NT;
  typedef _One_root_point_2<NT>                            Point_2;
  typedef typename Kernel::Circle_2                        Circle_2;
  typedef typename Kernel::Segment_2                       Segment_2;
  typedef typename Kernel::Line_2                          Line_2;

private:

  Line_2        _line;        // The supporting line (for line segments).
  Circle_2      _circ;        // The supporting circle (for circular arc).
  Point_2       _source;      // The source point.
  Point_2       _target;      // The target point.
  Orientation   _orient;      // The orientation (COLLINEAR for line segments).
  bool          _is_full;     // Whether we have a full circle.

public:

  /*! Default constructor. */
  _Circle_segment_2 () :
    _orient (COLLINEAR),
    _is_full (false)
  {}

  /*! Constructor from a line segment.
   * \param seg The segment.
   */
  _Circle_segment_2 (const Segment_2& seg) :
    _line (seg),
    _source (seg.source().x(), seg.source().y()),
    _target (seg.target().x(), seg.target().y()),
    _orient (COLLINEAR),
    _is_full (false)
  {}

  /*!
   * Constructor of a segment, given a supporting line and two endpoints,
   * which need not necessarily have rational coordinates.
   * \param line The supporting line.
   * \param source The source point.
   * \param target The target point.
   * \pre Both endpoints lie on the supporting line.
   */
  _Circle_segment_2 (const Line_2& line,
                     const Point_2& source, const Point_2& target) :
    _line (line),
    _source (source),
    _target (target),
    _orient (COLLINEAR),
    _is_full (false)
  {
    CGAL_precondition (CGAL::sign (line.a()*source.x() +
                                   line.b()*source.y() +
                                   line.c()) == ZERO);

    CGAL_precondition (CGAL::sign (line.a()*target.x() +
                                   line.b()*target.y() +
                                   line.c()) == ZERO);
  }

  /*! Constructor from a circle. */
  _Circle_segment_2 (const Circle_2& circ) :
    _circ (circ),
    _orient (circ.orientation()),
    _is_full (true)
  {
    CGAL_assertion (_orient != COLLINEAR);
  }

  /*!
   * Constructor of a circular, given a supporting circle and two endpoints,
   * which need not necessarily have rational coordinates. The orientation of
   * the circle determines the orientation of the arc.
   * \param circ The supporting circle.
   * \param source The source point.
   * \param target The target point.
   * \pre Both endpoints lie on the supporting circle.
   */
  _Circle_segment_2 (const Circle_2& circ,
                     const Point_2& source, const Point_2& target) :
    _circ (circ),
    _source (source),
    _target (target),
    _orient (circ.orientation()),
    _is_full (false)
  {
    CGAL_precondition
      (CGAL::compare (CGAL::square (source.x() - circ.center().x()),
                      circ.squared_radius -
                      CGAL::square (source.y() - circ.center().y())) == EQUAL);

    CGAL_precondition
      (CGAL::compare (CGAL::square (target.x() - circ.center().x()),
                      circ.squared_radius -
                      CGAL::square (target.y() - circ.center().y())) == EQUAL);
  }

  /*!
   * Get the orientation of the curve. 
   * \return COLLINEAR in case of a line segment,
   *         CLOCKWISE or COUNTERCLOCKWISE for circular curves.
   */
  Orientation orientation () const
  {
    return (_orient);
  }

  /*!
   * Get the supporting line.
   * \pre The curve orientation is COLLINEAR.
   */
  const Line_2& line () const
  {
    CGAL_precondition (_orient == COLLINEAR);
    return (_line);
  }

  /*!
   * Get the supporting circle.
   * \pre The curve orientation is not COLLINEAR.
   */
  const Circle_2& circle () const
  {
    CGAL_precondition (_orient != COLLINEAR);
    return (_circ);
  }

  /*! Check if the curve is a full circle. */
  bool is_full () const
  {
    return (_is_full);
  }

  /*!
   * Get the source point.
   * \pre The curve is not a full circle.
   */
  const Point_2& source () const
  {
    CGAL_precondition (! _is_full);
    return (_source);
  }

  /*!
   * Get the target point.
   * \pre The curve is not a full circle.
   */
  const Point_2& target () const
  {
    CGAL_precondition (! _is_full);
    return (_target);
  }

  /*!
   * Get the vertical tangency points the arc contains.
   * \param vpts Output: The vertical tagnecy points.
   * \pre The curve is circular.
   * \return The number of points (0, 1, or 2).
   */
  unsigned int vertical_tangency_points (Point_2 *vpts) const
  {
    CGAL_precondition (_orient != COLLINEAR);
    unsigned int  n_vpts = 0;

    if (_is_full)
    {
      // In case of a full circle, create both vertical tangency points:
      const NT&                  x0 = _circ.center().x();
      const NT&                  y0 = _circ.center().y();
      typename Point_2::CoordNT  xv_left (x0, -1, _circ.squared_radius());
      typename Point_2::CoordNT  xv_right (x0, 1, _circ.squared_radius());

      vpts[0] = Point_2 (xv_left, y0);
      vpts[1] = Point_2 (xv_right, y0);
      return (2);
    }

    if (_orient == COUNTERCLOCKWISE)
    {
      // Compute the vertical tangency points for the arc:
      n_vpts = _ccw_vertical_tangency_points (_source, _target, vpts);
    }
    else
    {
      // Compute the vertical tangency points for the opposite arc:
      n_vpts = _ccw_vertical_tangency_points (_target, _source, vpts);

      // Swap their order, if necessary.
      if (n_vpts == 2)
      {
        Point_2   temp = vpts[0];
        vpts[0] = vpts[1];
        vpts[1] = temp;
      }
    }

    return (n_vpts);
  }

private:

  /*!
   * Get the vertical tangency points the arc contains, assuming it is
   * counterclockwise oreinted.
   * \param vpts Output: The vertical tagnecy points.
   * \return The number of points (0, 1, or 2).
   */
  unsigned int _ccw_vertical_tangency_points (const Point_2& src,
                                              const Point_2& trg,
                                              Point_2 *vpts) const
  {
    unsigned int  n_vpts = 0;
    const NT&     x0 = _circ.center().x();
    const NT&     y0 = _circ.center().y();
    int           qs = _quart_index (src);
    int           qt = _quart_index (trg);
  
    if (qs == qt)
    {
      if ((qs == 0 || qs == 1) && CGAL::compare (src.x(), trg.x()) == LARGER)
        // We have an x-monotone arc lying on the upper half of the circle:
        return (0);
      
      if ((qs == 2 || qs == 3) && CGAL::compare (src.x(), trg.x()) == SMALLER)
        // We have an x-monotone arc lying on the lower half of the circle:
        return (0);
    }

    // Make sure the target quarter is larger than the source quarter, by 
    // adding 4 to its index, if necessary.
    if (qt <= qs)
      qt += 4;

    // Start traversing the quarter-planes and collect the vertical tangency
    // points we encounter.
    while (qs < qt)
    {
      if ((qs % 4) == 1)
      {
        // We collect the left tangency point when going from Q[1] to Q[2]:
        if (CGAL::compare (y0, trg.y()) != EQUAL)
        {
          typename Point_2::CoordNT  xv_left (x0, -1, _circ.squared_radius());

          vpts[n_vpts] = Point_2 (xv_left, y0);
          n_vpts++;
        }
      }
      else if ((qs % 4) == 3)
      {
        // We collect the right tangency point when going from Q[3] to Q[0]:
        if (CGAL::compare (y0, trg.y()) != EQUAL)
        {
          typename Point_2::CoordNT  xv_right (x0, 1, _circ.squared_radius());
          
          vpts[n_vpts] = Point_2 (xv_right, y0);
          n_vpts++;
        }
      }

      qs++;
    }

    return (n_vpts);
  }

  /*!
   * Get the index of the quarter-plane containing the given point,
   * where the circle center is considered to be the origin.
   */
  int _quart_index (const Point_2& p) const
  {
    // The plane looks like:
    //
    //      Q[1] :  |   Q[0]:
    //      x <= 0  |   x >  0
    //      y >  0  |   y >= 0
    //    ----------+-----------
    //      Q[2] :  |   Q[3]:
    //      x <  0  |   x >= 0
    //      y <= 0  |   y <  0
    //
    const CGAL::Sign   sign_x = CGAL::sign (p.x() - _circ.center().x());
    const CGAL::Sign   sign_y = CGAL::sign (p.y() - _circ.center().y());

    if (sign_x == POSITIVE)
    {
      return ((sign_y == NEGATIVE) ? 3 : 0);
    }
    else if (sign_x == NEGATIVE)
    {
      return ((sign_y == POSITIVE) ? 1 : 2);
    }

    CGAL_assertion (sign_y != ZERO);
    return ((sign_y == POSITIVE) ? 0 : 2);
  }
};

/*! \class
 * Representation of an x-monotone circular arc.
 */

template <class Kernel_>
class _X_monotone_circle_segment_2
{
public:

  typedef Kernel_                                          Kernel;
  typedef _X_monotone_circle_segment_2<Kernel>             Self;
  typedef typename Kernel::FT                              NT;
  typedef _One_root_point_2<NT>                            Point_2;
  typedef typename Kernel::Circle_2                        Circle_2;
  typedef typename Kernel::Line_2                          Line_2;
  typedef typename Point_2::CoordNT                        CoordNT;

protected:

  NT           _first;       // The x-coordinate of the circle center.
                             // Or: the coefficient of x in the line equation.

  NT           _second;      // The y-coordinate of the circle center.
                             // Or: the coefficient of y in the line equation.

  NT           _third;       // The squared radius of the supporting circle.
                             // Or: the free coefficient in the line equation.
  
  Point_2      _source;      // The source point.
  Point_2      _target;      // The target point.
  Orientation  _orient;      // The orientation of the arc (COLLINEAR for line
                             // segments).
  bool         _dir_right;   // Is the arc directed from left to right.
  bool         _is_vert;     // Is this a vertical segment.

public:

  /*!
   * Default constructor.
   */
  _X_monotone_circle_segment_2 () :
    _first(), 
    _second(),
    _third(),
    _source(), _target(),
    _orient (COLLINEAR),
    _dir_right (false),
    _is_vert (false)
  {}

  /*!
   * Construct an arc from a line segment.
   * \param line The supporting line.
   * \param source The source point.
   * \param target The target point.
   */
  _X_monotone_circle_segment_2 (const Line_2& line,
                                const Point_2& source, const Point_2& target) :
    _first (line.a()),
    _second (line.b()),
    _third (line.c()),
    _source (source), 
    _target(target),
    _orient (COLLINEAR),
    _is_vert (false)
  {
    // Check if the segment is directed left or right:
    Comparison_result   res = CGAL::compare (source.x(), target.x());

    if (res == EQUAL)
    {
      CGAL_precondition (CGAL::sign(_second) == ZERO);

      // We have a vertical segment - compare the points by their
      // y-coordinates:
      _is_vert = true;
      res = CGAL::compare (source.y(), target.y());
    }

    CGAL_precondition (res != EQUAL);
    _dir_right = (res == SMALLER);
  }

  /*! 
   * Construct a circular arc.
   * \param line The supporting line.
   * \param source The source point.
   * \param target The target point.
   * \param orient The orientation of the arc.
   */
  _X_monotone_circle_segment_2 (const Circle_2& circ,
                                const Point_2& source, const Point_2& target,
                                Orientation orient) :
    _first (circ.center().x()),
    _second (circ.center().y()),
    _third (circ.squared_radius()),
    _source (source), 
    _target(target),
    _orient (orient),
    _is_vert (false)
  {
    // Check if the segment is directed left or right:
    Comparison_result   res = CGAL::compare (source.x(), target.x());

    CGAL_precondition (res != EQUAL);
    _dir_right = (res == SMALLER);
  }

  /*! Check if the arc is linear. */
  inline bool is_linear () const
  {
    return (_orient == COLLINEAR);
  }

  /*! Check if the arc is circular. */
  inline bool is_circular () const
  {
    return (_orient != COLLINEAR);
  }

  /*!
   * Get the supporting circle. 
   * \pre The arc is circular.
   */
  Circle_2 supporting_circle () const
  {
    typename Kernel::Point_2  center (x0(), y0());
    return (Circle_2 (center , sqr_r(), _orient));
  }

  /*! Get the source point. */
  inline const Point_2& source () const
  {
    return (_source);
  }

  /*! Get the target point. */
  inline const Point_2& target () const
  {
    return (_target);
  }

  /*! Get the left endpoint of the arc. */
  inline const Point_2& left () const
  {
    return (_dir_right ? _source : _target);
  }

  /*! Get the right endpoint of the arc. */
  inline const Point_2& right () const
  {
    return (_dir_right ? _target : _source);
  }

  /*!
   * Check whether the given point is in the x-range of the arc.
   */
  bool is_in_x_range (const Point_2& p) const
  {
    Comparison_result    res = CGAL::compare (p.x(), left().x());

    if (res == SMALLER)
      return (false);
    else if (res == EQUAL)
      return (true);

    return (CGAL::compare (p.x(), right().x()) != LARGER);
  }

  /*! Check if the arc is a vertical segment. */
  inline bool is_vertical () const
  {
    return (_is_vert);
  }

  /*!
   * Check the position of a given point with respect to the arc.
   */
  Comparison_result point_position (const Point_2& p) const
  {
    if (is_linear())
      return (_line_point_position (p));
    else
      return (_circ_point_position (p));
  }

  /*!
   * Compare the two arcs to the right of their intersection point.
   */
  Comparison_result compare_to_right (const Self& cv, const Point_2& p) const
  {
    if (is_linear())
    {
      if (cv.is_linear())
        return (_lines_compare_to_right (cv, p));
      
      Comparison_result   res = cv._circ_line_compare_to_right (*this, p);
      
      if (res != EQUAL)
        res = (res == SMALLER) ? LARGER : SMALLER;

      return (res);
    }
    else
    {
      if (cv.is_linear())
        return (_circ_line_compare_to_right (cv, p));

      return (_circs_compare_to_right (cv, p));
    }
  }

  /*!
   * Compare the two arcs to the left of their intersection point.
   */
  Comparison_result compare_to_left (const Self& cv, const Point_2& p) const
  {
    if (is_linear())
    {
      if (cv.is_linear())
        return (_lines_compare_to_left (cv, p));
      
      Comparison_result   res = cv._circ_line_compare_to_left (*this, p);
      
      if (res != EQUAL)
        res = (res == SMALLER) ? LARGER : SMALLER;

      return (res);
    }
    else
    {
      if (cv.is_linear())
        return (_circ_line_compare_to_left (cv, p));

      return (_circs_compare_to_left (cv, p));
    }
  }

  /*!
   * Check whether the two arcs have the same supporting curve.
   */
  bool has_same_supporting_curve (const Self& cv) const
  {
    // Make sure that the supporting curves are of the same type.
    if (is_linear() && ! cv.is_linear())
      return (false);

    if (! is_linear() && cv.is_linear())
      return (false);

    // Compare the curve coefficients.
    return (CGAL::compare (_first, cv._first) == EQUAL &&
            CGAL::compare (_second, cv._second) == EQUAL &&
            CGAL::compare (_third, cv._third) == EQUAL);
  }

  /*!
   * Check if the two curves are equal.
   */
  bool equals (const Self& cv) const
  {
    if (! this->has_same_supporting_curve (cv))
      return (false);

    if (is_linear())
    {
      // In case of line segments we can swap the source and target:
      return ((_source.equals (cv._source) && _target.equals (cv._target)) ||
              (_source.equals (cv._target) && _target.equals (cv._source)));
    }

    // Once again, opposite circular arcs are considered to be equal:
    return ((_orient == cv._orient &&
             _source.equals (cv._source) && _target.equals (cv._target)) ||
            (_orient != cv._orient &&
             _source.equals (cv._target) && _target.equals (cv._source)));
  }

  /*!
   * Split the curve at a given point into two sub-arcs.
   */
  void split (const Point_2& p, Self& c1, Self& c2) const
  {
    // Copy the properties of this arc to the sub-arcs.
    c1 = *this;
    c2 = *this;

    // Change the endpoint, such that c1 lies to the right of c2:
    if (_dir_right)
    {
      c1._target = p;
      c2._source = p;
    }
    else
    {
      c1._source = p;
      c2._target = p;
    }

    return;
  }

  /*!
   * Compute the intersections between the two circles.
   */
  template <class OutputIterator>
  OutputIterator intersect (const Self& cv, OutputIterator oi) const
  {
    // First check whether the two arcs have the same supporting curve.
    if (has_same_supporting_curve (cv))
    {
      // Check for overlaps between the two arcs.
      Self    overlap;

      if (_compute_overlap (cv, overlap))
      {
        // There can be just a single overlap between two x-monotone arcs:
        *oi = CGAL::make_object (overlap);
        ++oi;
        return (oi);
      }

      // In case there is not overlap and the supporting curves are the same,
      // there cannot be any intersection points, unless the two arcs share
      // a common end point.
      // Note that in this case we do not define the multiplicity of the
      // intersection points we report.
      unsigned int  mult = 0;
      if (left().equals (cv.left()))
      {
        *oi = CGAL::make_object (std::make_pair (left(), mult));
        ++oi;
      }

      if (right().equals (cv.right()))
      {
        *oi = CGAL::make_object (std::make_pair (right(), mult));
        ++oi;
      }

      return (oi);
    }

    // Compute the intersections between the two arcs.
    if (is_linear())
    {
      if (cv.is_linear())
        return (_lines_intersect (cv, oi));
      
      return (cv._circ_line_intersect (*this, oi, true));
    }
    else
    {
      if (cv.is_linear())
        return (_circ_line_intersect (cv, oi, false));

      return (_circs_intersect (cv, oi));
    }
  }

  /*!
   * Check whether it is possible to merge our arc with the given arc.
   */
  bool can_merge_with (const Self& cv) const
  {
    // In order to merge the two arcs, they should have the same supporting
    // curve.
    if (! _has_same_supporting_curve (cv))
      return (false);

    // Check if the left endpoint of one curve is the right endpoint of the
    // other.
    return (right().equals (cv.left()) ||
            left().equals (cv.right()));
  }

  /*!
   * Merge our arc with the given arc.
   * \pre The two arcs are mergeable.
   */
  void merge (const Self& cv)
  {
    CGAL_precondition (this->can_merge_with (cv));

    // Check if we should extend the arc to the left or to the right.
    if (right().equals (cv.left()))
    {
      // Extend the arc to the right.
      if (_dir_right)
        this->_target = cv.right();
      else
        this->_source = cv.right();
    }
    else
    {
      CGAL_precondition (left().equals (cv.right()));

      // Extend the arc to the left.
      if (_dir_right)
        this->_source = cv.left();
      else
        this->_target = cv.left();
    }

    return;
  }

protected:

  /// \name Accessors for circular arcs.
  //@{

  /*! Get the x-coordinate of the center of the supporting circle. */
  inline const NT& x0 () const
  {
    return (_first);
  }

  /*! Get the y-coordinate of the center of the supporting circle. */
  inline const NT& y0 () const
  {
    return (_second);
  }

  /*! Get the squared radius of the supporting circle. */
  inline const NT& sqr_r () const
  {
    return (_third);
  }

  /*!
   * Check if the circular arc lies on the upper half of the supporting circle.
   */
  inline bool _is_upper () const
  {
    CGAL_precondition (_orient != COLLINEAR);

    return ((_orient == COUNTERCLOCKWISE && !_dir_right) |
            (_orient == CLOCKWISE && _dir_right));
  }
  //@}

  /// \name Accessors for line segments.
  //@{

  /*! Get the coefficient of x in the equation of the supporting line. */
  inline const NT& a () const
  {
    return (_first);
  }

  /*! Get the coefficient of y in the equation of the supporting line. */
  inline const NT& b () const
  {
    return (_second);
  }

  /*! Get the free coefficient in the equation of the supporting line. */
  inline const NT& c () const
  {
    return (_third);
  }
  //@}

  /// \name Auxiliary functions for the point_position predicate.
  //@{

  /*!
   * Check the position of a given point with respect to a line segment.
   */
  Comparison_result _line_point_position (const Point_2& p) const
  {
    // Check if we have a vertical segment.
    Comparison_result    res;

    if (is_vertical())
    {
      // left() is the lower endpoint:
      res = CGAL::compare (p.y(), left().y());

      if (res != LARGER)
        return (res);

      // left() is the upper endpoint:
      res = CGAL::compare (p.y(), right().y());

      if (res != SMALLER)
        return (res);

      // p lies in the interior of the vertical segment:
      return (EQUAL);
    }

    // Compute the y-coordinate of the vertical projection of p onto the
    // supporting line.
    const CoordNT        y_proj = (a()*p.x() + c()) / (-b());

    return (CGAL::compare (p.y(), y_proj));
  }

  /*!
   * Check the position of a given point with respect to a circular arc.
   */
  Comparison_result _circ_point_position (const Point_2& p) const
  {
    Comparison_result   c_res = CGAL::compare (p.y(), y0());

    if (_is_upper())
    {
      // Check if p lies below the "equator" (while the arc lies above it):
      if (c_res == SMALLER)
        return (SMALLER);
    }
    else
    {
      // Check if p lies above the "equator" (while the arc lies below it):
      if (c_res == LARGER)
        return (LARGER);
    }

    // Check if p lies inside the supporting circle, namely we have to check
    // whether (p.x() - x0)^2 + (p.y() - y0)^2 < r^2:
    Comparison_result   res =
                         CGAL::compare (CGAL::square (p.x() - x0()),
                                        sqr_r() - CGAL::square (p.y() - y0()));

    if (res == EQUAL)
      // p lies on the circle:
      return (EQUAL);

    if (_is_upper())
    {
      // If p is inside the circle, it lies below the upper arc:
      return (res);
    }
    else
    {
      // If p is inside the circle, it lies above the lower arc:
      return (res == SMALLER ? LARGER : SMALLER);
    }
  }
  //@}

  /// \name Auxiliary functions for the compare_to_right predicate.
  //@{

  /*!
   * Compare two line segments to the right of their intersection point.
   */
  Comparison_result _lines_compare_to_right (const Self& cv,
                                             const Point_2& p) const
  {
    // Special treatment for vertical segments: a vertical segment is larger
    // than any other non-vertical segment.
    if (is_vertical())
    {
      if (cv.is_vertical())
        return (EQUAL);

      return (LARGER);
    }
    else if (cv.is_vertical())
    {
      return (SMALLER);
    }

    // Compare the slopes: -A1/B1 and -A2/B2. We actually negate the slopes
    // and swap the result.
    return (CGAL::compare (cv.a()/cv.b(), a()/b()));
  }

  /*!
   * Compare a circular arcs (this) and a line segment (cv) to the right of
   * their intersection point.
   */
  Comparison_result _circ_line_compare_to_right (const Self& cv,
                                                 const Point_2& p) const
  {
    // A vertical segment lies above any other circle to the right of p:
    if (cv.is_vertical())
      return (SMALLER);
    
    // We have to compare the slopes of the supporting circles and the
    // supporting line at p:
    //
    //    p.x() - x0(1)           A(2)
    //   ---------------  and  - ------
    //    y0(1) - p.y()           B(2)
    //
    const CGAL::Sign  sign_denom1 = CGAL::sign (y0() - p.y());

    // Check the case of a vertical tangent.
    if (sign_denom1 == ZERO)
    {
      // The arc lies above any line segment if it is an upper arc, or below
      // any segment if it is a lower arc.
      return (_is_upper() ? LARGER : SMALLER);
    }

    // Compare (p.x() - x0(1)) and (A(2)/B(2)*(p.y() - y0(1)).
    // Note that if the denominator is negative, we have to swap the result.
    const bool        swap_res = (sign_denom1 == NEGATIVE);
    Comparison_result slope_res = CGAL::compare (p.x() - x0(),
                                                 (p.y() - y0())*cv.a()/cv.b());

    if (slope_res != EQUAL)
    {
      if (swap_res)
        // Swap the comparison result, if necessary:
        slope_res = (slope_res == SMALLER) ? LARGER : SMALLER;
      
      return (slope_res);
    }

    // In this case we have a tangency point at p. If the circular arc is an
    // upper arc, it must lie below the tangent line, and if it is a lower arc
    // it must lie above the tangent line.
    return (_is_upper() ? SMALLER : LARGER);
  }

  /*!
   * Compare two circular arcs to the right of their intersection point.
   */
  Comparison_result _circs_compare_to_right (const Self& cv,
                                             const Point_2& p) const
  {
    // We have to compare the slopes of the two supporting circles at p:
    //
    //    p.x() - x0(1)         p.x() - x0(2)
    //   ---------------  and  ---------------
    //    y0(1) - p.y()         y0(2) - p.y()
    //
    const CGAL::Sign  sign_numer1 = CGAL::sign (p.x() - x0());
    const CGAL::Sign  sign_denom1 = CGAL::sign (y0() - p.y());
    const CGAL::Sign  sign_numer2 = CGAL::sign (p.x() - cv.x0());
    const CGAL::Sign  sign_denom2 = CGAL::sign (cv.y0() - p.y());

    // Check the case of vertical tangents.
    if (sign_denom1 == ZERO)
    {
      if (sign_denom2 == ZERO)
      {
        if (_is_upper())
        {
          if (cv._is_upper())
          {
            // The two circles have a vertical tangent:
            // The one with a larger radius is above the other.
            return (CGAL::compare (sqr_r(), cv.sqr_r()));
          }
          else
          {
            // The other curve is directed downwards:
            return (LARGER);
          }
        }
        else
        {
          if (cv._is_upper())
          {
            // The other curve is directed upwards:
            return (SMALLER);
          }
          else
          {
            // The two circles have a vertical tangent:
            // The one with a smaller radius is above the other.
            return (CGAL::compare (cv.sqr_r(), sqr_r()));
          }
        }
      }

      // The other arc does not have a vertical tangent.
      return (_is_upper() ? LARGER : SMALLER);
    }
    else if (sign_denom2 == ZERO)
    {
      return (cv._is_upper() ? SMALLER : LARGER);
    }

    // Try to act according to the slope signs.
    CGAL::Sign   sign_slope1;
    CGAL::Sign   sign_slope2;

    if (sign_numer1 == sign_denom1)
      sign_slope1 = POSITIVE;
    else if (sign_numer1 == ZERO)
      sign_slope1 = ZERO;
    else
      sign_slope1 = NEGATIVE;

    if (sign_numer2 == sign_denom2)
      sign_slope2 = POSITIVE;
    else if (sign_numer2 == ZERO)
      sign_slope2 = ZERO;
    else
      sign_slope2 = NEGATIVE;

    if ((sign_slope1 == POSITIVE && sign_slope2 != POSITIVE) ||
        (sign_slope1 == ZERO && sign_slope2 == NEGATIVE))
      return (LARGER);

    if ((sign_slope2 == POSITIVE && sign_slope1 != POSITIVE) ||
        (sign_slope2 == ZERO && sign_slope1 == NEGATIVE))
      return (SMALLER);

    // Compare the slopes of the two tangents to the circles.
    Comparison_result  slope_res;
    
    if (sign_slope1 == ZERO && sign_slope2 == ZERO)
    {
      // Special case were both circles have a horizontal tangent:
      slope_res = EQUAL;
    }
    else
    {
      // Actually compare the slopes.    
      const bool    swap_res = (sign_denom1 != sign_denom2);
      const CoordNT A = (cv.y0() - y0())*p.x() + (y0()*cv.x0() - cv.y0()*x0());
      const CoordNT B = (cv.x0() - x0())*p.y();
     
      slope_res = CGAL::compare (A, B);

      if (slope_res != EQUAL && swap_res)
      {
        // Swap the comparison result, if necessary:
        slope_res = (slope_res == SMALLER) ? LARGER : SMALLER;
      }
    }

    // In case the two circles have different tangent slopes at p:
    if (slope_res != EQUAL)
      return (slope_res);

    // In this case we have a tangency point at p.
    if (_is_upper())
    {
      if (cv._is_upper())
      {
        // The circle with a larger radius is above the other.
        return (CGAL::compare (sqr_r(), cv.sqr_r()));
      }
      else
      {
        // The other curve is above our curve:
        return (SMALLER);
      }
    }
    else
    {
      if (cv._is_upper())
      {
        // Out curve is above the other curve:
        return (LARGER);
      }
      else
      {
        // The circle with a smaller radius is above the other.
        return (CGAL::compare (cv.sqr_r(), sqr_r()));
      }
    }
  }
  //@}

  /// \name Auxiliary functions for the compare_to_left predicate.
  //@{

  /*!
   * Compare two line segments to the left of their intersection point.
   */
  Comparison_result _lines_compare_to_left (const Self& cv,
                                             const Point_2& p) const
  {
    // Special treatment for vertical segments: a vertical segment is smaller
    // than any other non-vertical segment.
    if (is_vertical())
    {
      if (cv.is_vertical())
        return (EQUAL);

      return (SMALLER);
    }
    else if (cv.is_vertical())
    {
      return (LARGER);
    }

    // Compare the slopes: -A1/B1 and -A2/B2 and swap the result.
    //  We actually negate the slopes and compare them.
    return (CGAL::compare (a()/b(), cv.a()/cv.b()));
  }

  /*!
   * Compare a circular arcs (this) and a line segment (cv) to the left of
   * their intersection point.
   */
  Comparison_result _circ_line_compare_to_left (const Self& cv,
                                                const Point_2& p) const
  {
    // A vertical segment lies below any other circle to the left of p:
    if (cv.is_vertical())
      return (LARGER);
    
    // We have to compare the slopes of the supporting circles and the
    // supporting line at p, and return the swapped result:
    //
    //    p.x() - x0(1)           A(2)
    //   ---------------  and  - ------
    //    y0(1) - p.y()           B(2)
    //
    const CGAL::Sign  sign_denom1 = CGAL::sign (y0() - p.y());

    // Check the case of a vertical tangent.
    if (sign_denom1 == ZERO)
    {
      // The arc lies above any line segment if it is an upper arc, or below
      // any segment if it is a lower arc.
      return (_is_upper() ? LARGER : SMALLER);
    }

    // Compare (p.x() - x0(1)) and (A(2)/B(2)*(p.y() - y0(1)).
    // Note that if the denominator is negative, we have to swap the result.
    const bool        swap_res = (sign_denom1 == NEGATIVE);
    Comparison_result slope_res = CGAL::compare (p.x() - x0(),
                                                 (p.y() - y0())*cv.a()/cv.b());

    if (slope_res != EQUAL)
    {
      if (swap_res)
        // Swap the comparison result, if necessary:
        slope_res = (slope_res == SMALLER) ? LARGER : SMALLER;
      
      // Swap at any case to get the position to the left:
      return ((slope_res == SMALLER) ? LARGER : SMALLER);
    }

    // In this case we have a tangency point at p. If the circular arc is an
    // upper arc, it must lie below the tangent line, and if it is a lower arc
    // it must lie above the tangent line.
    return (_is_upper() ? SMALLER : LARGER);
  }

  /*!
   * Compare the two arcs to the left of their intersection point.
   */
  Comparison_result _circs_compare_to_left (const Self& cv, 
                                            const Point_2& p) const
  {
    // We have to compare the slopes of the two supporting circles at p:
    //
    //    p.x() - x0(1)         p.x() - x0(2)
    //   ---------------  and  ---------------
    //    y0(1) - p.y()         y0(2) - p.y()
    //
    // Eventually, we should take the opposite result.
    const CGAL::Sign  sign_numer1 = CGAL::sign (p.x() - x0());
    const CGAL::Sign  sign_denom1 = CGAL::sign (y0() - p.y());
    const CGAL::Sign  sign_numer2 = CGAL::sign (p.x() - cv.x0());
    const CGAL::Sign  sign_denom2 = CGAL::sign (cv.y0() - p.y());

    // Check the case of vertical tangents.
    if (sign_denom1 == ZERO)
    {
      if (sign_denom2 == ZERO)
      {
        if (_is_upper())
        {
          if (cv._is_upper())
          {
            // The two circles have a vertical tangent:
            // The one with a larger radius is above the other.
            return (CGAL::compare (sqr_r(), cv.sqr_r()));
          }
          else
          {
            // The other curve is directed downwards:
            return (LARGER);
          }
        }
        else
        {
          if (cv._is_upper())
          {
            // The other curve is directed upwards:
            return (SMALLER);
          }
          else
          {
            // The two circles have a vertical tangent:
            // The one with a smaller radius is above the other.
            return (CGAL::compare (cv.sqr_r(), sqr_r()));
          }
        }
      }

      // The other arc does not have a vertical tangent.
      return (_is_upper() ? LARGER : SMALLER);
    }
    else if (sign_denom2 == ZERO)
    {
      return (cv._is_upper() ? SMALLER : LARGER);
    }

    // Try to act according to the slope signs.
    CGAL::Sign   sign_slope1;
    CGAL::Sign   sign_slope2;

    if (sign_numer1 == sign_denom1)
      sign_slope1 = POSITIVE;
    else if (sign_numer1 == ZERO)
      sign_slope1 = ZERO;
    else
      sign_slope1 = NEGATIVE;

    if (sign_numer2 == sign_denom2)
      sign_slope2 = POSITIVE;
    else if (sign_numer2 == ZERO)
      sign_slope2 = ZERO;
    else
      sign_slope2 = NEGATIVE;

    if ((sign_slope1 == POSITIVE && sign_slope2 != POSITIVE) ||
        (sign_slope1 == ZERO && sign_slope2 == NEGATIVE))
      return (SMALLER);

    if ((sign_slope2 == POSITIVE && sign_slope1 != POSITIVE) ||
        (sign_slope2 == ZERO && sign_slope1 == NEGATIVE))
      return (LARGER);

    // Compare the slopes of the two tangents to the circles.
    Comparison_result  slope_res;
    
    if (sign_slope1 == ZERO && sign_slope2 == ZERO)
    {
      // Special case were both circles have a horizontal tangent:
      slope_res = EQUAL;
    }
    else
    {
      // Actually compare the slopes.    
      const bool    swap_res = (sign_denom1 != sign_denom2);
      const CoordNT A = (cv.y0() - y0())*p.x() + (y0()*cv.x0() - cv.y0()*x0());
      const CoordNT B = (cv.x0() - x0())*p.y();
     
      slope_res = CGAL::compare (A, B);

      if (slope_res != EQUAL && swap_res)
      {
        // Swap the comparison result, if necessary:
        slope_res = (slope_res == SMALLER) ? LARGER : SMALLER;
      }
    }

    // In case the two circles have different tangent slopes at p, return
    // the opposite of the slope result (since the slope result is the
    // comparison result to the right of the intersection point):
    if (slope_res != EQUAL)
      return ((slope_res == SMALLER) ? LARGER : SMALLER);

    // In this case we have a tangency point at p.
    if (_is_upper())
    {
      if (cv._is_upper())
      {
        // The circle with a larger radius is above the other.
        return (CGAL::compare (sqr_r(), cv.sqr_r()));
      }
      else
      {
        // The other curve is above our curve:
        return (SMALLER);
      }
    }
    else
    {
      if (cv._is_upper())
      {
        // Out curve is above the other curve:
        return (LARGER);
      }
      else
      {
        // The circle with a smaller radius is above the other.
        return (CGAL::compare (cv.sqr_r(), sqr_r()));
      }
    }
  }
  //@}

  /// \name Auxiliary functions for computing intersections.
  //@{

  /*!
   * Compute the intersections between two line segments.
   */
  template <class OutputIterator>
  OutputIterator _lines_intersect (const Self& cv, OutputIterator oi) const
  {
    // The intersection of the lines:
    //   a1*x + b1*y + c1 = 0   and   a2*x + b2*y + c2 = 0 
    // Is given by:
    //
    //      b1*c2 - c1*b2     c1*a2 - a1*c2
    //   ( --------------- , --------------- )
    //      a1*b2 - b1*a2     a1*b2 - b1*a2
    //
    unsigned int  mult = 1;
    const NT      denom = a()*cv.b() - b()*cv.a();

    if (CGAL::sign(denom) == ZERO)
    {
      // The supporting lines are parallel - no intersections.
      return (oi);
    }

    const NT      x_numer = b()*cv.c() - c()*cv.b();
    const NT      y_numer = c()*cv.a() - a()*cv.c();
    Point_2       p (x_numer / denom, y_numer / denom);

    // Check if the point lies on both segments. If so, we found an
    // intersection point of multiplicity 1.
    if (this->_is_between_endpoints (p) &&
        cv._is_between_endpoints (p))
    {
      *oi = CGAL::make_object (std::make_pair (p, mult));
      ++oi;
    }

    return (oi);
  }

  /*!
   * Compute the intersections between a circular arc (this) and a line
   * segment (cv).
   */
  template <class OutputIterator>
  OutputIterator _circ_line_intersect (const Self& cv, OutputIterator oi,
                                       bool roles_swapped) const
  {
    // Intersect the supporting circle and the supporting line.
    unsigned int    mult;
    Point_2         ps[2];
    unsigned int    n_ps = _intersect_supporting_circ_line (cv, ps);
    unsigned int    k;

    if (n_ps == 0)
      return (oi);
    
    if (n_ps == 1)
    {
      // We found a single tangency point, with multiplicity 2:
      mult = 2;
    }
    else
    {
      CGAL_assertion (n_ps == 2);

      // We found two intersection points with multiplicity 1 each:
      mult = 1;
    }

    // Report just the points that lies on both arcs.
    for (k = 0; k < n_ps; k++)
    {
      // Note that we check the first curve sent to intersect() first.
      // This is why we keep track whether the roles of the curves have
      // been swapped.
      if (! roles_swapped)
      {
        if (this->_is_between_endpoints (ps[k]) &&
            cv._is_between_endpoints (ps[k]))
        {
          *oi = CGAL::make_object (std::make_pair (ps[k], mult));
          ++oi;
        }
      }
      else
      {
        if (cv._is_between_endpoints (ps[k]) &&
            this->_is_between_endpoints (ps[k]))
        {
          *oi = CGAL::make_object (std::make_pair (ps[k], mult));
          ++oi;
        }
      }
    }

    return (oi);
  }

  /*!
   * Compute the intersections between two circular arcs.
   */
  template <class OutputIterator>
  OutputIterator _circs_intersect (const Self& cv, OutputIterator oi) const
  {
    // Intersect the two supporting circles.
    unsigned int    mult;
    Point_2         ps[2];
    unsigned int    n_ps = _intersect_supporting_circles (cv, ps);
    unsigned int    k;

    if (n_ps == 0)
      return (oi);
    
    if (n_ps == 1)
    {
      // We found a single tangency point, with multiplicity 2:
      mult = 2;
    }
    else
    {
      CGAL_assertion (n_ps == 2);

      // We found two intersection points with multiplicity 1 each:
      mult = 1;
    }

    // Report just the points that lies on both arcs.
    for (k = 0; k < n_ps; k++)
    {
      if (this->_is_between_endpoints (ps[k]) &&
          cv._is_between_endpoints (ps[k]))
      {
        *oi = CGAL::make_object (std::make_pair (ps[k], mult));
        ++oi;
      }
    }

    return (oi);
  }

  /*!
   * Compute the intersection points between the supporting circle of (*this)
   * and the supporting line of cv.
   */
  unsigned int _intersect_supporting_circ_line (const Self& cv,
                                                Point_2* ps) const
  {
    // First check the special cases of vertical and horizontal lines.
    if (cv.is_vertical())
    {
      // The equation of the vertical line is x = -c / a.
      // The y-coordinates of the intersection points are:
      //   y =  y0 +/- sqrt(r^2 - (x - x0)^2)
      //
      const NT   vx = -cv.c() / cv.a();
      const NT   vdisc = sqr_r() - CGAL::square (vx - x0());
      CGAL::Sign sign_vdisc = CGAL::sign (vdisc);
      
      if (sign_vdisc == NEGATIVE)
      {
        // The circle and the vertical line do not intersect.
        return (0);
      }
      else if (sign_vdisc == ZERO)
      {
        // A single tangency point, given by:
        ps[0] = Point_2 (vx, y0());
        return (1);
      }
      
      // Compute the two intersection points:
      ps[0] = Point_2 (CoordNT (vx),
                       CoordNT (y0(), -1, vdisc));
      ps[1] = Point_2 (CoordNT (vx),
                       CoordNT (y0(), 1, vdisc));
      return (2);
    }
    else if (CGAL::sign (cv.a()) == ZERO)
    {
      // The equation of the vertical line is y = -c / b.
      // The y-coordinates of the intersection points are:
      //   x =  x0 +/- sqrt(r^2 - (y - y0)^2)
      //
      const NT   hy = -cv.c() / cv.b();
      const NT   hdisc = sqr_r() - CGAL::square (hy - y0());
      CGAL::Sign sign_hdisc = CGAL::sign (hdisc);
      
      if (sign_hdisc == NEGATIVE)
      {
        // The circle and the vertical line do not intersect.
        return (0);
      }
      else if (sign_hdisc == ZERO)
      {
        // A single tangency point, given by:
        ps[0] = Point_2 (x0(), hy);
        return (1);
      }
      
      // Compute the two intersection points:
      ps[0] = Point_2 (CoordNT (x0(), -1, hdisc),
                       CoordNT (hy));
      ps[1] = Point_2 (CoordNT (x0(), 1, hdisc),
                       CoordNT (hy));
      return (2);
    }

    // Compute the squared distance between the line and the circle center,
    // inducing the discriminant of the quadratic equations we have to solve.
    const NT   line_factor = CGAL::square(cv.a()) + CGAL::square(cv.b());
    const NT   disc = line_factor*sqr_r() -
                      CGAL::square(cv.a()*x0() + cv.b()*y0() + cv.c());
    CGAL::Sign sign_disc = CGAL::sign (disc);

    if (sign_disc == NEGATIVE)
      // The circle and the line do not intersect:
      return (0);
 
    // Compare the square-free part of the solution:
    const NT   aux = cv.b()*x0() - cv.a()*y0();
    const NT   x_base = (aux - cv.a()*cv.c()) / line_factor;
    const NT   y_base = (-aux - cv.b()*cv.c()) / line_factor;
    Point_2    p1, p2;

    if (sign_disc == ZERO)
    {
      // A single tangency point, given by:
      ps[0] = Point_2 (x_base, y_base);
      return (1);
    }

    // We have two intersection points, whose coordinates are one-root numbers.
    bool       minus_root_first = (CGAL::sign(cv.b()) == POSITIVE);
    const NT   x_root_coeff = cv.b() / line_factor;
    const NT   y_root_coeff = cv.a() / line_factor;

    if (minus_root_first)
    {
      ps[0] = Point_2 (CoordNT (x_base, -x_root_coeff, disc),
                       CoordNT (y_base, y_root_coeff, disc));
      ps[1] = Point_2 (CoordNT (x_base, x_root_coeff, disc),
                       CoordNT (y_base, -y_root_coeff, disc));
    }
    else
    {
      ps[0] = Point_2 (CoordNT (x_base, x_root_coeff, disc),
                       CoordNT (y_base, -y_root_coeff, disc));
      ps[1] = Point_2 (CoordNT (x_base, -x_root_coeff, disc),
                       CoordNT (y_base, y_root_coeff, disc));
    }

    return (2);
  }

  /*!
   * Compute the intersection points between the supporting circles of the
   * two arcs.
   */
  unsigned int _intersect_supporting_circles (const Self& cv,
                                              Point_2* ps) const
  {
    // Compute the squared distance between the circle centers, inducing the
    // discriminant of the quadratic equations we have to solve.
    const NT   diff_x = cv.x0() - x0();
    const NT   diff_y = cv.y0() - y0();    
    const NT   sqr_dist = CGAL::square(diff_x) + CGAL::square(diff_y);
    const NT   diff_sqr_rad = sqr_r() - cv.sqr_r();
    const NT   disc = 2*sqr_dist*(sqr_r() + cv.sqr_r()) -
                      (CGAL::square(diff_sqr_rad) + CGAL::square(sqr_dist));
    CGAL::Sign sign_disc = CGAL::sign (disc);

    if (sign_disc == NEGATIVE)
      // The two circles do not intersect.
      return (0);
 
    // Compare the square-free part of the solution:
    const NT   x_base = ((x0() + cv.x0()) + diff_x*diff_sqr_rad/sqr_dist) / 2;
    const NT   y_base = ((y0() + cv.y0()) + diff_y*diff_sqr_rad/sqr_dist) / 2;
    Point_2    p1, p2;

    if (sign_disc == ZERO)
    {
      // A single tangency point, given by:
      ps[0] = Point_2 (x_base, y_base);
      return (1);
    }

    // We have two intersection points, whose coordinates are one-root numbers.
    CGAL::Sign sign_diff_y = CGAL::sign (diff_y);
    bool       minus_root_first;

    if (sign_diff_y == ZERO)
      minus_root_first = (CGAL::sign (diff_x) == NEGATIVE);
    else
      minus_root_first = (sign_diff_y == POSITIVE);

    const NT   x_root_coeff = diff_y / (2 * sqr_dist);
    const NT   y_root_coeff = diff_x / (2 * sqr_dist);

    if (minus_root_first)
    {
      ps[0] = Point_2 (CoordNT (x_base, -x_root_coeff, disc),
                       CoordNT (y_base, y_root_coeff, disc));
      ps[1] = Point_2 (CoordNT (x_base, x_root_coeff, disc),
                       CoordNT (y_base, -y_root_coeff, disc));
    }
    else
    {
      ps[0] = Point_2 (CoordNT (x_base, x_root_coeff, disc),
                       CoordNT (y_base, -y_root_coeff, disc));
      ps[1] = Point_2 (CoordNT (x_base, -x_root_coeff, disc),
                       CoordNT (y_base, y_root_coeff, disc));
    }

    return (2);
  }

  /*!
   * Check if the given point lies on the arc.
   * \pre p lies on the supporting curve.
   */
  bool _is_between_endpoints (const Point_2& p) const
  {
    if (is_linear())
    {
      if (is_vertical())
      {
        // Check if the point is in the y-range of the arc.
        // Note that left() is the lower endpoint and right() is the upper
        // endpoint of the segment in this case. 
        Comparison_result    res = CGAL::compare (p.y(), left().y());

        if (res == SMALLER)
          return (false);
        else if (res == EQUAL)
          return (true);

        return (CGAL::compare (p.y(), right().y()) != LARGER);
      }

      // For non-vertical segments, it is sufficient to check if the point
      // is in the x-range of the arc.
      return (this->is_in_x_range (p));
    }

    // The supporting curve is a circle:
    // Check whether p lies on the upper or on the lower part of the circle.
    Comparison_result   c_res = CGAL::compare (p.y(), y0());

    if ((_is_upper() && c_res == SMALLER) ||
        (! _is_upper() && c_res == LARGER))
    {
      // The point lies on the other half of the circle:
      return (false);
    }

    // Check if the point is in the x-range of the arc.
    return (this->is_in_x_range (p));
  }

  /*!
   * Check if the given point lies in the interior of the arc.
   * \pre p lies on the supporting curve.
   */
  bool _is_strictly_between_endpoints (const Point_2& p) const
  {
    if (p.equals (_source) || p.equals (_target))
      return (false);
    
    return (_is_between_endpoints (p));
  }

  /*!
   * Compute the overlap with a given arc having the same supporting curve.
   * \param cv The given arc.
   * \param overlap Output: The overlapping arc (if any).
   * \return Whether we found an overlap.
   */
  bool _compute_overlap (const Self& cv, Self& overlap) const
  {
    // Check if the two arcs are identical.
    if (is_linear())
    {
      // In case of line segments we can swap the source and target:
      if (((_source.equals (cv._source) && _target.equals (cv._target)) ||
           (_source.equals (cv._target) && _target.equals (cv._source))))
      {
        overlap = cv;
        return (true);
      }
    }
    else
    {
      if ((_orient == cv._orient &&
           _source.equals (cv._source) && _target.equals (cv._target)) ||
          (_orient != cv._orient &&
           _source.equals (cv._target) && _target.equals (cv._source)))
      {
        overlap = cv;
        return (true);
      }
    }

    // Check for other overlaps:
    if (_is_strictly_between_endpoints (cv.left()))
    {
      if (_is_strictly_between_endpoints (cv.right()))
      {
        // Case 1 - *this:     +----------->
        //             cv:       +=====>
        overlap = cv;
        return (true);
      }
      else
      {
        // Case 2 - *this:     +----------->
        //             cv:               +=====>
        overlap = *this;

        if (overlap._dir_right)
          overlap._source = cv.left();
        else
          overlap._target = cv.left();

        return (true);
      }
    }
    else if (_is_strictly_between_endpoints (cv.right()))
    {
      // Case 3 - *this:     +----------->
      //             cv:   +=====>
      overlap = *this;

      if (overlap._dir_right)
        overlap._target = cv.right();
      else
        overlap._source = cv.right();

      return (true);
    }
    else if (cv._is_between_endpoints (_source) &&
             cv._is_between_endpoints (_target) &&
             (cv._is_strictly_between_endpoints (_source) ||
              cv._is_strictly_between_endpoints (_target)))
    {
      // Case 4 - *this:     +----------->
      //             cv:   +================>
      overlap = *this;
      return (true);
    }

    // If we reached here, there are no overlaps:
    return (false);
  }

  //@}
};

/*!
 * Exporter for circular arcs (or line segments).
 */
template <class Kernel>
std::ostream& 
operator<< (std::ostream& os, 
            const _X_monotone_circle_segment_2<Kernel> & arc)
{
  if (! arc.is_linear())
    os << "(" << arc.supporting_circle() << ") ";

  os << "[" << arc.source() << " --> " << arc.target() << "]" << std::endl;
  return (os);
}

CGAL_END_NAMESPACE

#endif
