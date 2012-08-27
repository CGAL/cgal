// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Ron Wein        <wein@post.tau.ac.il>

#ifndef CGAL_HYPERBOLIC_ARC_2_H
#define CGAL_HYPERBOLIC_ARC_2_H

/*! \file
 * Header file for the _Hyperbolic_arc_2<Kernel, Filter> class.
 */
#include <CGAL/Arr_traits_2/Circle_segment_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Handle_for.h>
#include <list>
#include <map>
#include <ostream>

namespace CGAL {

/*! \class
 * Representation of an x-monotone hyprbolic arc, defined by:
 *   y = a*x^2 + b*x + c, for x <= x_min < = x_max.
 * Note that a, b and c are rational numbers, while x_min and x_max
 * may be one-root numbers.
 */
template <class Kernel_, bool Filter_>
class _Hyperbolic_arc_2
{
public:

  typedef Kernel_                                          Kernel;
  typedef _Hyperbolic_arc_2<Kernel, Filter_>               Self;
  typedef typename Kernel::FT                              NT;
  typedef _One_root_point_2<NT, Filter_>                   Point_2;
  typedef typename Point_2::CoordNT                        CoordNT;

protected:

  NT           _a;           // The coefficient of x^2.
  NT           _b;           // The coefficient of x.
  NT           _c;           // The free coefficient.
  
  Point_2      _source;      // The source point.
  Point_2      _target;      // The target point.

public:

  /*!
   * Default constructor.
   */
  _Hyperbolic_arc_2 () :
    _a(), _b(), _c(),
    _source(), _target()
  {}

  /*!
   * Construct an arc.
   * \param a, b, c The coefficients of the supporting hyperbola.
   * \param x_min, x_max Define the x-range of the hyperbolic arc.
   */
  _Hyperbolic_arc_2 (const NT& a, const NT& b, const NT& c,
                     const CoordNT& x_min, const CoordNT& x_max) :
    _a (a),
    _b (b),
    _c (c)
  {
    CGAL_precondition (CGAL::compare (x_min, x_max) == CGAL::SMALLER);

    // Compute the endpoints.
    _source = Point_2 (x_min, _get_y (x_min));
    _target = Point_2 (x_max, _get_y (x_max));
  }

  /*!
   * Construct a segment arc from two kernel points
   * \param source the source point.
   * \ param target the target point.
   * \pre The source is lexicographically smaller than the target.
   */
  _Hyperbolic_arc_2 (const typename Kernel::Point_2& source,
                     const typename Kernel::Point_2& target) :
    _source(source.x(), source.y()),
    _target(target.x(), target.y())
  {
    CGAL_precondition (CGAL::compare(source.x(), target.x()) == CGAL::SMALLER);
 
    // Set the coefficients of the supporting curve
    typename Kernel::Line_2 line(source, target);
    _a = 0;
    _b = - line.a() / line.b();
    _c = - line.c() / line.b();    
  }

  /*! Check if the arc is linear. */
  inline bool is_linear () const
  {
    return (CGAL::sign(_a) == CGAL::ZERO);
  }

  /*! Check if the arc is hyperbolic. */
  inline bool is_hyperbolic () const
  {
    return (CGAL::sign(_a) != CGAL::ZERO);
  }

  /*!
   * Get the coefficients of the supporting curve.
   */
  const NT& a () const
  {
    return (_a);
  }

  const NT& b () const
  {
    return (_b);
  }

  const NT& c () const
  {
    return (_c);
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
    return (_source);
  }

  /*! Get the right endpoint of the arc. */
  inline const Point_2& right () const
  {
    return (_target);
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
    return (CGAL::compare (p.y(), _get_y (p.x())));
  }

  /*!
   * Compare the two arcs to the right of their intersection point.
   */
  Comparison_result compare_to_right (const Self& cv, const Point_2& p) const
  {
    // Compute the first-order derivatives of both curves at p. 
    const CoordNT      der1 = p.x() * 2 * _a + _b;
    const CoordNT      der2 = p.x() * 2 * cv._a + cv._b;
    Comparison_result  res = CGAL::compare (der1, der2);

    // In case of inequality, return the comparison reult.
    if (res != CGAL::EQUAL)
      return (res);

    // In case of equality, compare the second-order derivatives.
    return (CGAL::compare (_a, cv._a));
  }

  /*!
   * Compare the two arcs to the left of their intersecton point.
   */
  Comparison_result compare_to_left (const Self& cv, const Point_2& p) const
  {
    // Compute the first-order derivatives of both curves at p. 
    const CoordNT      der1 = p.x() * 2 * _a + _b;
    const CoordNT      der2 = p.x() * 2 * cv._a + cv._b;
    Comparison_result  res = CGAL::compare (der1, der2);

    // In case of inequality, negate the comparison reult.
    if (res != CGAL::EQUAL)
      return (CGAL::opposite (res));

    // In case of equality, compare the second-order derivatives.
    return (CGAL::compare (_a, cv._a));
  }

  /*!
   * Check whether the two arcs have the same supporting curve.
   */
  bool has_same_supporting_curve (const Self& cv) const
  {
    return (CGAL::compare (_a, cv._a) == CGAL::EQUAL &&
            CGAL::compare (_b, cv._b) == CGAL::EQUAL &&
            CGAL::compare (_c, cv._c) == CGAL::EQUAL);
  }

  /*!
   * Check if the two curves are equal.
   */
  bool equals (const Self& cv) const
  {
    if (! this->has_same_supporting_curve (cv))
      return (false);

    return (left().equals (cv.left()) && right().equals (cv.right()));
  }

  /*!
   * Split the curve at a given point into two sub-arcs.
   */
  void split (const Point_2& p, Self& c1, Self& c2) const
  {
    // Copy the properties of this arc to the sub-arcs.
    c1 = *this;
    c2 = *this;

    // Change the endpoints, such that c1 lies to the right of c2:
    c1._target = p;
    c2._source = p;

    return;
  }

  /*!
   * Compute the intersections between the two arcs or segments.
   */
  template <class OutputIterator>
  OutputIterator intersect (const Self& cv, OutputIterator oi) const
  {
    // Solve the quadratic equation A*x^2 + B*x + C = 0 in order to find
    // the x-coordinates of the intersection points, where:
    const NT       A = _a - cv._a;
    const NT       B = _b - cv._b;
    const NT       C = _c - cv._c;
    Point_2        p;
    unsigned int   mult;

    // Check if we have a linear equation.
    if (CGAL::sign (A) == CGAL::ZERO)
    {
      if (CGAL::sign (B) == ZERO)
      {
        if (CGAL::sign (C) == ZERO)
        {
          // Here we have to handle overlaps!
          //CGAL_error();
        }

        return (oi);
      }

      // We have a single (rational) intersection point.
      const NT    x = -C / B;
      const NT    y = (x * _a + _b) * x + _c;

      p = Point_2 (CoordNT(x), CoordNT(y));
      mult = 1;

      *oi = CGAL::make_object (std::make_pair (p, mult));
      ++oi;

      return (oi);
    }

    // In this case we have to solve a quadratic equation.
    const NT       disc = B*B - 4*A*C;
    CGAL::Sign     sign_disc = CGAL::sign (disc);

    if (sign_disc == CGAL::NEGATIVE)
      // No intersection:
      return (oi);

    const NT       _1_over_2A = 1 / (2*A);

    if (sign_disc == CGAL::ZERO)
    {
      // We have a single tangency point with rational coordinates.
      const NT    x = -B * _1_over_2A;
      const NT    y = (x * _a + _b) * x + _c;

      p = Point_2 (CoordNT(x), CoordNT(y));
      mult = 2;

      *oi = CGAL::make_object (std::make_pair (p, mult));
      ++oi;

      return (oi);
    }

    // In this case we have two solutions, given by:
    CoordNT         xs[2];
    int             k;

    if (CGAL::sign (A) == CGAL::POSITIVE)
    {
      xs[0] = CoordNT (-B * _1_over_2A, - _1_over_2A, disc);
      xs[1] = CoordNT (-B * _1_over_2A, _1_over_2A, disc);
    }
    else
    {
      xs[0] = CoordNT (-B * _1_over_2A, _1_over_2A, disc);
      xs[1] = CoordNT (-B * _1_over_2A, - _1_over_2A, disc);
    }

    for (k = 0; k < 2; k++)
    {
      // Check if the x-coordinate is in the x-range of both arcs.
      if ((CGAL::compare (xs[k], left().x()) != SMALLER &&
           CGAL::compare (xs[k], right().x()) != LARGER) &&
          (CGAL::compare (xs[k], cv.left().x()) != SMALLER &&
           CGAL::compare (xs[k], cv.right().x()) != LARGER))
      {
        p = Point_2 (xs[k], _get_y (xs[k]));
        mult = 1;

        *oi = CGAL::make_object (std::make_pair (p, mult));
        ++oi;
      }
    }

    return (oi);
  }

  Bbox_2 bbox() const
  {
    double x_min = CGAL::to_double (left().x());
    double x_max = CGAL::to_double (right().x());   
    double y_min = CGAL::to_double (left().y()); 
    double y_max = CGAL::to_double (right().y());

    if(y_min > y_max)
      std::swap(y_min, y_max);
    
    return Bbox_2(x_min, y_min, x_max, y_max);
  }

private:

  /*!
   * Compute the y-coordiate at a given x-coordinate.
   */
  CoordNT _get_y (const CoordNT& x) const
  {
    if (x.is_rational())
      return ((_a * x.alpha() + _b) * x.alpha() + _c);

    const CoordNT   z1 = _a * CGAL::square(x);
    const CoordNT   z2 = _b * x;

    return (CoordNT (z1.alpha() + z2.alpha() + _c,
                     z1.beta() + z2.beta(),
                     z1.gamma()));
  }
};

/*!
 * Exporter for circular arcs (or line segments).
 */
template <class Kernel, bool Filter>
std::ostream& 
operator<< (std::ostream& os, 
            const _Hyperbolic_arc_2<Kernel, Filter> & arc)
{
  if (! arc.is_linear())
    os << "(" << arc.a() << "*x^2 + " << arc.b() << "*x + " << arc.c() << ") ";

  os << "[" << arc.source() << " --> " << arc.target() << "]" << std::endl;
  return (os);
}

} //namespace CGAL

#endif
