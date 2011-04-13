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

#include <CGAL/Cartesian/redefine_names_2.h>
#include <CGAL/cartesian_classes.h>
#include <CGAL/Cartesian/Point_2.h>
#include <CGAL/Cartesian/Line_2.h>
#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
bool
equal_line(const LineC2<R CGAL_CTAG> &l1, const LineC2<R CGAL_CTAG> &l2)
{
  return equal_lineC2(l1.a(), l1.b(), l1.c(), l2.a(), l2.b(), l2.c());
}

template < class R >
inline
Comparison_result
compare_x(const PointC2<R CGAL_CTAG> &p,
          const LineC2<R CGAL_CTAG> &l,
          const LineC2<R CGAL_CTAG> &h)
{
  return compare_xC2(p.x(), l.a(), l.b(), l.c(), h.a(), h.b(), h.c());
}

template < class R >
inline
Comparison_result
compare_x(const LineC2<R CGAL_CTAG> &l,
          const LineC2<R CGAL_CTAG> &h1,
          const LineC2<R CGAL_CTAG> &h2)
{
  return compare_xC2(l.a(), l.b(), l.c(), h1.a(), h1.b(), h1.c(),
                                          h2.a(), h2.b(), h2.c());
}

template < class R >
inline
Comparison_result
compare_x(const LineC2<R CGAL_CTAG> &l1,
          const LineC2<R CGAL_CTAG> &h1,
          const LineC2<R CGAL_CTAG> &l2,
          const LineC2<R CGAL_CTAG> &h2)
{
  return compare_xC2(l1.a(), l1.b(), l1.c(), h1.a(), h1.b(), h1.c(),
                     l2.a(), l2.b(), l2.c(), h2.a(), h2.b(), h2.c());
}

template < class R >
inline
Comparison_result
compare_y(const PointC2<R CGAL_CTAG> &p,
          const LineC2<R CGAL_CTAG> &l1,
          const LineC2<R CGAL_CTAG> &l2)
{
  return compare_xC2(p.y(), l1.b(), l1.a(), l1.c(), l2.b(), l2.a(), l2.c());
}

template < class R >
inline
Comparison_result
compare_y(const LineC2<R CGAL_CTAG> &l1,
          const LineC2<R CGAL_CTAG> &l2,
          const LineC2<R CGAL_CTAG> &h1,
          const LineC2<R CGAL_CTAG> &h2)
{
  return compare_xC2(l1.b(), l1.a(), l1.c(), l2.b(), l2.a(), l2.c(),
                     h1.b(), h1.a(), h1.c(), h2.b(), h2.a(), h2.c());
}

template < class R >
inline
Comparison_result
compare_y(const LineC2<R CGAL_CTAG> &l,
          const LineC2<R CGAL_CTAG> &h1,
          const LineC2<R CGAL_CTAG> &h2)
{
  return compare_xC2(l.b(), l.a(), l.c(), h1.b(), h1.a(), h1.c(),
                     l.b(), l.a(), l.c(), h2.b(), h2.a(), h2.c());
}

template < class R >
inline
Comparison_result
compare_y_at_x(const PointC2<R CGAL_CTAG> &p, const LineC2<R CGAL_CTAG> &h)
{
  return compare_y_at_xC2(p.x(), p.y(), h.a(), h.b(), h.c());
}

template < class R >
inline
Comparison_result
compare_y_at_x(const PointC2<R CGAL_CTAG> &p,
               const LineC2<R CGAL_CTAG> &h1,
               const LineC2<R CGAL_CTAG> &h2)
{
  return compare_y_at_xC2(p.x(), h1.a(), h1.b(), h1.c(),
	                         h2.a(), h2.b(), h2.c());
}

template < class R >
inline
Comparison_result
compare_y_at_x(const LineC2<R CGAL_CTAG> &l1,
               const LineC2<R CGAL_CTAG> &l2,
               const LineC2<R CGAL_CTAG> &h)
{
  return compare_y_at_xC2(l1.a(), l1.b(), l1.c(), l2.a(), l2.b(), l2.c(),
                          h.a(), h.b(), h.c());
}

template < class R >
inline
Comparison_result
compare_y_at_x(const LineC2<R CGAL_CTAG> &l1,
               const LineC2<R CGAL_CTAG> &l2,
               const LineC2<R CGAL_CTAG> &h1,
               const LineC2<R CGAL_CTAG> &h2)
{
  return compare_y_at_xC2(l1.a(), l1.b(), l1.c(), l2.a(), l2.b(), l2.c(),
                          h1.a(), h1.b(), h1.c(), h2.a(), h2.b(), h2.c());
}

template < class R >
inline
Comparison_result
compare_x_at_y(const PointC2<R CGAL_CTAG> &p, const LineC2<R CGAL_CTAG> &h)
{
  return compare_y_at_xC2(p.y(), p.x(), h.b(), h.a(), h.c());
}

template < class R >
inline
Comparison_result
compare_x_at_y(const PointC2<R CGAL_CTAG> &p,
               const LineC2<R CGAL_CTAG> &h1,
               const LineC2<R CGAL_CTAG> &h2)
{
  return compare_y_at_xC2(p.y(), h1.b(), h1.a(), h1.c(),
	                         h2.b(), h2.a(), h2.c());
}

template < class R >
inline
Comparison_result
compare_x_at_y(const LineC2<R CGAL_CTAG> &l1,
               const LineC2<R CGAL_CTAG> &l2,
               const LineC2<R CGAL_CTAG> &h)
{
  return compare_y_at_xC2(l1.b(), l1.a(), l1.c(), l2.b(), l2.a(), l2.c(),
                          h.b(), h.a(), h.c());
}

template < class R >
inline
Comparison_result
compare_x_at_y(const LineC2<R CGAL_CTAG> &l1,
               const LineC2<R CGAL_CTAG> &l2,
               const LineC2<R CGAL_CTAG> &h1,
               const LineC2<R CGAL_CTAG> &h2)
{
  return compare_y_at_xC2(l1.b(), l1.a(), l1.c(), l2.b(), l2.a(), l2.c(),
                          h1.b(), h1.a(), h1.c(), h2.b(), h2.a(), h2.c());
}

template < class R >
inline
Oriented_side
side_of_oriented_line(const LineC2<R CGAL_CTAG> &l,
                      const PointC2<R CGAL_CTAG> &p)
{
  return side_of_oriented_lineC2(l.a(), l.b(), l.c(), p.x(), p.y());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PREDICATES_ON_LINES_2_H
