// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_PREDICATES_ON_LINES_2_H
#define CGAL_CARTESIAN_PREDICATES_ON_LINES_2_H

#include <CGAL/Cartesian/Point_2.h>
#include <CGAL/Cartesian/Line_2.h>
#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

template < class K >
inline
bool
equal_line(const LineC2<K> &l1, const LineC2<K> &l2)
{
  return equal_lineC2(l1.a(), l1.b(), l1.c(), l2.a(), l2.b(), l2.c());
}

template < class K >
inline
Comparison_result
compare_x(const PointC2<K> &p,
          const LineC2<K> &l,
          const LineC2<K> &h)
{
  return K().compare_x_2_object()(p, l, h);
}

template < class K >
inline
Comparison_result
compare_x(const LineC2<K> &l,
          const LineC2<K> &h1,
          const LineC2<K> &h2)
{
  return K().compare_x_2_object()(l, h1, h2);
}

template < class K >
inline
Comparison_result
compare_x(const LineC2<K> &l1,
          const LineC2<K> &h1,
          const LineC2<K> &l2,
          const LineC2<K> &h2)
{
  return K().compare_x_2_object()(l1, l2, h1, h2);
}

template < class K >
inline
Comparison_result
compare_y(const PointC2<K> &p,
          const LineC2<K> &l1,
          const LineC2<K> &l2)
{
  return K().compare_y_2_object()(p, l1, l2);
}

template < class K >
inline
Comparison_result
compare_y(const LineC2<K> &l1,
          const LineC2<K> &l2,
          const LineC2<K> &h1,
          const LineC2<K> &h2)
{
  return K().compare_y_2_object()(l1, l2, h1, h2);
}

template < class K >
inline
Comparison_result
compare_y(const LineC2<K> &l,
          const LineC2<K> &h1,
          const LineC2<K> &h2)
{
  return K().compare_y_2_object()(l, h1, h2);
}  

template < class K >
inline
Comparison_result
compare_y_at_x(const PointC2<K> &p, const LineC2<K> &h)
{
  return K().compare_y_at_x_2_object()(p, h);
}  

template < class K >
inline
Comparison_result
compare_y_at_x(const PointC2<K> &p,
               const LineC2<K> &h1,
               const LineC2<K> &h2)
{
  return K().compare_y_at_x_2_object()(p, h1, h2);
}

template < class K >
inline
Comparison_result
compare_y_at_x(const LineC2<K> &l1,
               const LineC2<K> &l2,
               const LineC2<K> &h)
{
  return K().compare_y_at_x_2_object()(l1, l2, h);
}

template < class K >
inline
Comparison_result
compare_y_at_x(const LineC2<K> &l1,
               const LineC2<K> &l2,
               const LineC2<K> &h1,
               const LineC2<K> &h2)
{
  return K().compare_y_at_x_2_object()(l1, l2, h1, h2);
}

template < class K >
inline
Comparison_result
compare_x_at_y(const PointC2<K> &p, const LineC2<K> &h)
{
  return K().compare_x_at_y_2_object()(p, h);
}

template < class K >
inline
Comparison_result
compare_x_at_y(const PointC2<K> &p,
               const LineC2<K> &h1,
               const LineC2<K> &h2)
{
  return K().compare_x_at_y_2_object()(p, h1, h2);
}

template < class K >
inline
Comparison_result
compare_x_at_y(const LineC2<K> &l1,
               const LineC2<K> &l2,
               const LineC2<K> &h)
{
  return K().compare_x_at_y_2_object()(l1, l2, h);
}

template < class K >
inline
Comparison_result
compare_x_at_y(const LineC2<K> &l1,
               const LineC2<K> &l2,
               const LineC2<K> &h1,
               const LineC2<K> &h2)
{
  return K().compare_x_at_y_2_object()(l1, l2, h1, h2);
}

template < class K >
inline
Comparison_result
compare_slopes(const LineC2<K> &l1,
               const LineC2<K> &l2)
{
  return K().compare_slope_2_object()(l1, l2);
}

template < class K >
inline
Oriented_side
side_of_oriented_line(const LineC2<K> &l,
                      const PointC2<K> &p)
{
  return side_of_oriented_lineC2(l.a(), l.b(), l.c(), p.x(), p.y());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PREDICATES_ON_LINES_2_H
