// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/predicates_on_lines_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

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
