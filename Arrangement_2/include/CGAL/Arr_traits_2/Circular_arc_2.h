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
// $URL$
// $Id$ $Date$
// 
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_CIRCULAR_ARC_2_H
#define CGAL_CIRCULAR_ARC_2_H

/*! \file
 * Header file for the _Circular_arc_2<Kernel> class.
 */

#include <CGAL/Arr_traits_2/One_root_number.h>
#include <CGAL/Bbox_2.h>
#include <ostream>

CGAL_BEGIN_NAMESPACE

/*! \class
 * Representation of a point on a circular arc.
 */
template <class NumberType_>
class _Circle_point_2
{
public:

  typedef NumberType_                 NT;
  typedef _One_root_number<NT>        CoordNT;

private:

  CoordNT     _x;
  CoordNT     _y;

public:

  /*! Default constructor. */
  _Circle_point_2 () :
    _x (0),
    _y (0)
  {}

  /*! Constructor of a point with rational coefficients. */
  _Circle_point_2 (const NT& x, const NT& y) :
    _x (x),
    _y (y)
  {}

  /*! Constructor of a point with one-root coefficients. */
  _Circle_point_2 (const CoordNT& x, const CoordNT& y) :
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
  bool equals (const _Circle_point_2<NT>& p) const
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
            const _Circle_point_2<NT>& p)
{
  os << CGAL::to_double(p.x()) << ' ' << CGAL::to_double(p.y());
  return (os);
}

/*! \class
 * Representation of a circle with a rational radius.
 */
template <class Kernel_>
class _Rat_circle_2
{
public:

  typedef Kernel_                                          Kernel;
  typedef typename Kernel::FT                              NT;
  typedef typename Kernel::Point_2                         Point_2;

private:

  Point_2       _c;
  NT            _r;

public:

  /*! Default constructor. */
  _Rat_circle_2 () :
    _c (),
    _r (0)
  {}

  /*! Constructor of a circle given its center and radius. */
  _Rat_circle_2 (const Point_2& center, const NT& radius) :
    _c (center),
    _r (radius)
  {
    CGAL_precondition (CGAL::sign(radius) != NEGATIVE);
  }

  /*! Get the circle center. */
  const Point_2& center () const
  {
    return (_c);
  }

  /*! Get the radius. */
  const NT& radius () const
  {
    return (_r);
  }
};

/*! \class
 * Representation of an x-monotone circular arc.
 */

template <class Kernel_>
class _Circular_arc_2
{
public:

  typedef Kernel_                                          Kernel;
  typedef _Circular_arc_2<Kernel>                          Self;
  typedef typename Kernel::FT                              NT;
  typedef _Rat_circle_2<Kernel>                            Circle_2;
  typedef _Circle_point_2<NT>                              Point_2;

protected:

  typedef typename Kernel::Point_2                         Rat_point_2;
  typedef typename Point_2::CoordNT                        CoordNT;

protected:

  NT           _x0;        // The x-coordinate of the circle center.
  NT           _y0;        // The y-coordinate of the circle center.
  NT           _sqr_rad;   // The squared radius of the supporting circle.
  Point_2      _source;    // The source point.
  Point_2      _target;    // The target point.
  bool         _is_upper;  // Is the arc on the upper half of the circle.
  bool         _dir_right; // Is the arc directed from left to right.

public:

  /*!
   * Default constructor.
   */
  _Circular_arc_2 () :
    _x0(), _y0(),
    _sqr_rad(),
    _source(), _target(),
    _is_upper(false),
    _dir_right(false)
  {}

  /*! 
   * Construct a conic which corresponds to half a circle.
   * \param circ The circle.
   * \param is_upper Does the arc correspond to the upper or to the lower
   *                 half of the circle.
   */
  _Circular_arc_2 (const Circle_2& circ, bool is_upper) :
    _x0 (circ.center().x()),
    _y0 (circ.center().y()),
    _sqr_rad (circ.radius() * circ.radius()),
    _is_upper (is_upper)
  {
    CGAL_precondition (CGAL::sign (circ.radius()) == POSITIVE);

    // The two endpoints correspond to (x0 + rad, y0) and (x0 - rad, y0):
    if (is_upper)
    {
      // Upper arc, directed counterclockwise from "3 o'clock" to "9 o'clock".
      _source = Point_2 (_x0 + circ.radius(), _y0);
      _target = Point_2 (_x0 - circ.radius(), _y0);
      _dir_right = false;
    }
    else
    {
      // Lower arc, directed counterclockwise from "9 o'clock" to "3 o'clock".
      _source = Point_2 (_x0 - circ.radius(), _y0);
      _target = Point_2 (_x0 + circ.radius(), _y0);
      _dir_right = true;
    }
  }

  /*! Get the circle center. */
  Rat_point_2 center () const
  {
    return (Rat_point_2 (_x0, _y0));
  }

  /*! Get the squared radius. */
  const NT& squared_radius () const
  {
    return (_sqr_rad);
  }

  /*! Get the source point. */
  const Point_2& source () const
  {
    return (_source);
  }

  /*! Get the target point. */
  const Point_2& target () const
  {
    return (_target);
  }

  /*! Get the left endpoint of the arc. */
  const Point_2& left () const
  {
    return (_dir_right ? _source : _target);
  }

  /*! Get the right endpoint of the arc. */
  const Point_2& right () const
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

  /*!
   * Check the position of a given point with respect to the arc.
   */
  Comparison_result point_position (const Point_2& p) const
  {
    Comparison_result   c_res = CGAL::compare (p.y(), _y0);

    if (_is_upper)
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
                         CGAL::compare (CGAL::square (p.x() - _x0),
                                        _sqr_rad - CGAL::square (p.y() - _y0));

    if (res == EQUAL)
      // p lies on the circle:
      return (EQUAL);

    if (_is_upper)
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

  /*!
   * Compare the two arcs to the right of their intersection point.
   */
  Comparison_result compare_to_right (const Self& cv, const Point_2& p) const
  {
    // We have to compare the slopes of the two supporting circles at p:
    //
    //    p.x() - x0(1)         p.x() - x0(2)
    //   ---------------  and  ---------------
    //    y0(1) - p.y()         y0(2) - p.y()
    //
    const CGAL::Sign  sign_numer1 = CGAL::sign (p.x() - _x0);
    const CGAL::Sign  sign_denom1 = CGAL::sign (_y0 - p.y());
    const CGAL::Sign  sign_numer2 = CGAL::sign (p.x() - cv._x0);
    const CGAL::Sign  sign_denom2 = CGAL::sign (cv._y0 - p.y());

    // Check the case of vertical tangents.
    if (sign_denom1 == ZERO)
    {
      if (sign_denom2 == ZERO)
      {
        if (_is_upper)
        {
          if (cv._is_upper)
          {
            // The two circles have a vertical tangent:
            // The one with a larger radius is above the other.
            return (CGAL::compare (_sqr_rad, cv._sqr_rad));
          }
          else
          {
            // The other curve is directed downwards:
            return (LARGER);
          }
        }
        else
        {
          if (cv._is_upper)
          {
            // The other curve is directed upwards:
            return (SMALLER);
          }
          else
          {
            // The two circles have a vertical tangent:
            // The one with a smaller radius is above the other.
            return (CGAL::compare (cv._sqr_rad, _sqr_rad));
          }
        }
      }

      // The other arc does not have a vertical tangent.
      return (_is_upper ? LARGER : SMALLER);
    }
    else if (sign_denom2 == ZERO)
    {
      return (cv._is_upper ? SMALLER : LARGER);
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
      const CoordNT A = (cv._y0 - _y0)*p.x() + (_y0*cv._x0 - cv._y0*_x0);
      const CoordNT B = (cv._x0 - _x0)*p.y();
     
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
    if (_is_upper)
    {
      if (cv._is_upper)
      {
        // The circle with a larger radius is above the other.
        return (CGAL::compare (_sqr_rad, cv._sqr_rad));
      }
      else
      {
        // The other curve is above our curve:
        return (SMALLER);
      }
    }
    else
    {
      if (cv._is_upper)
      {
        // Out curve is above the other curve:
        return (LARGER);
      }
      else
      {
        // The circle with a smaller radius is above the other.
        return (CGAL::compare (cv._sqr_rad, _sqr_rad));
      }
    }
  }

  /*!
   * Compare the two arcs to the left of their intersection point.
   */
  Comparison_result compare_to_left (const Self& cv, const Point_2& p) const
  {
    // We have to compare the slopes of the two supporting circles at p:
    //
    //    p.x() - x0(1)         p.x() - x0(2)
    //   ---------------  and  ---------------
    //    y0(1) - p.y()         y0(2) - p.y()
    //
    // Eventually, we should take the opposite result.
    const CGAL::Sign  sign_numer1 = CGAL::sign (p.x() - _x0);
    const CGAL::Sign  sign_denom1 = CGAL::sign (_y0 - p.y());
    const CGAL::Sign  sign_numer2 = CGAL::sign (p.x() - cv._x0);
    const CGAL::Sign  sign_denom2 = CGAL::sign (cv._y0 - p.y());

    // Check the case of vertical tangents.
    if (sign_denom1 == ZERO)
    {
      if (sign_denom2 == ZERO)
      {
        if (_is_upper)
        {
          if (cv._is_upper)
          {
            // The two circles have a vertical tangent:
            // The one with a larger radius is above the other.
            return (CGAL::compare (_sqr_rad, cv._sqr_rad));
          }
          else
          {
            // The other curve is directed downwards:
            return (LARGER);
          }
        }
        else
        {
          if (cv._is_upper)
          {
            // The other curve is directed upwards:
            return (SMALLER);
          }
          else
          {
            // The two circles have a vertical tangent:
            // The one with a smaller radius is above the other.
            return (CGAL::compare (cv._sqr_rad, _sqr_rad));
          }
        }
      }

      // The other arc does not have a vertical tangent.
      return (_is_upper ? LARGER : SMALLER);
    }
    else if (sign_denom2 == ZERO)
    {
      return (cv._is_upper ? SMALLER : LARGER);
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
      const CoordNT A = (cv._y0 - _y0)*p.x() + (_y0*cv._x0 - cv._y0*_x0);
      const CoordNT B = (cv._x0 - _x0)*p.y();
     
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
    if (_is_upper)
    {
      if (cv._is_upper)
      {
        // The circle with a larger radius is above the other.
        return (CGAL::compare (_sqr_rad, cv._sqr_rad));
      }
      else
      {
        // The other curve is above our curve:
        return (SMALLER);
      }
    }
    else
    {
      if (cv._is_upper)
      {
        // Out curve is above the other curve:
        return (LARGER);
      }
      else
      {
        // The circle with a smaller radius is above the other.
        return (CGAL::compare (cv._sqr_rad, _sqr_rad));
      }
    }
  }

  /*!
   * Check whether the two arcs have the sam supporting circle.
   */
  bool has_same_supporting_circle (const Self& cv) const
  {
    return (CGAL::compare (_x0, cv._x0) == EQUAL &&
            CGAL::compare (_y0, cv._y0) == EQUAL &&
            CGAL::compare (_sqr_rad, cv._sqr_rad) == EQUAL);
  }

  /*!
   * Check if the two curves are equal.
   */
  bool equals (const Self& cv) const
  {
    return (this->has_same_supporting_circle (cv) &&
            _is_upper == cv._is_upper &&
            (_dir_right == cv._dir_right &&
             _source.equals (cv._source) && _target.equals (cv._target)) ||
            (_dir_right != cv._dir_right &&
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
    unsigned int    mult;

    if (this->has_same_supporting_circle (cv))
    {
      // RWRW: To be done!
      mult = 0;
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

    // Intersect the two supporting circles.
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

protected:

  /*!
   * Compute the intersection points between the supporting circles of the
   * two arcs.
   */
  unsigned int _intersect_supporting_circles (const Self& cv,
                                              Point_2* ps) const
  {
    // Compute the squared distance between the circle centers.
    const NT   diff_x = cv._x0 - _x0;
    const NT   diff_y = cv._y0 - _y0;    
    const NT   sqr_dist = CGAL::square(diff_x) + CGAL::square(diff_y);
    const NT   diff_sqr_rad = _sqr_rad - cv._sqr_rad;
    const NT   disc = 2*sqr_dist*(_sqr_rad + cv._sqr_rad) -
                      (CGAL::square(diff_sqr_rad) + CGAL::square(sqr_dist));
    CGAL::Sign sign_disc = CGAL::sign (disc);

    if (sign_disc == NEGATIVE)
      // The two circles do not intersect.
      return (0);
 
    // Compare the square-free part of the solution:
    const NT   x_base = ((_x0 + cv._x0) + diff_x*diff_sqr_rad / sqr_dist) / 2;
    const NT   y_base = ((_y0 + cv._y0) + diff_y*diff_sqr_rad / sqr_dist) / 2;
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
   * \pre p lies on the supporting circle.
   */
  bool _is_between_endpoints (const Point_2& p) const
  {
    // Check whether p lies on the upper or on the lower part of the circle.
    Comparison_result   c_res = CGAL::compare (p.y(), _y0);

    if ((_is_upper && c_res == SMALLER) ||
        (! _is_upper && c_res == LARGER))
    {
      // The point lies on the other half of the circle:
      return (false);
    }

    // Check if the point is in the x-range of the arc.
    return (this->is_in_x_range (p));
  }
};

/*!
 * Exporter for conic arcs.
 */
template <class Kernel>
std::ostream& 
operator<< (std::ostream& os, 
            const _Circular_arc_2<Kernel> & arc)
{
  os << "(" << arc.center() << " , " << arc.squared_radius() << ") ["
     << arc.source() << " -> " << arc.target() << "]" << std::endl;
  return (os);
}

CGAL_END_NAMESPACE

#endif
