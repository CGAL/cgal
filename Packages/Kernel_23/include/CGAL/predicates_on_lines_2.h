// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : predicates_on_lines_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_PREDICATES_ON_LINES_2_H
#define CGAL_PREDICATES_ON_LINES_2_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#ifndef CGAL_PREDICATES_ON_LINESH2_H
#include <CGAL/predicates_on_linesH2.h>
#endif // CGAL_PREDICATES_ON_LINESH2_H
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#ifndef CGAL_PREDICATES_ON_LINESC2_H
#include <CGAL/Cartesian/predicates_on_lines_2.h>
#endif // CGAL_PREDICATES_ON_LINESC2_H
#endif // CGAL_CARTESIAN_H

#ifdef CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/SimpleCartesian/predicates_on_linesS2.h>
#endif // CGAL_SIMPLE_CARTESIAN_H


#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H
#ifndef CGAL_PREDICATES_ON_POINTS_2_H
#include <CGAL/predicates_on_points_2.h>
#endif // CGAL_PREDICATES_ON_POINTS_2_H
#ifndef CGAL_LINE_2_H
#include <CGAL/Line_2.h>
#endif // CGAL_LINE_2_H

CGAL_BEGIN_NAMESPACE


template < class R >
inline
Comparison_result
compare_x(const Point_2<R> &p,
               const Line_2<R> &l1,
               const Line_2<R> &l2)
{
  typedef typename  R::Point_2_base  RPoint_2;
  typedef typename  R::Line_2_base  RLine_2 ;
  return compare_x((const RPoint_2&)p,
                   (const RLine_2&)l1,
                   (const RLine_2&)l2);
}

template < class R >
inline
Comparison_result
compare_x(const Line_2<R> &l1,
          const Line_2<R> &l2,
          const Line_2<R> &h1,
          const Line_2<R> &h2)
{
  typedef typename  R::Line_2_base  RLine_2 ;
  return compare_x((const RLine_2&)l1, (const RLine_2&)l2,
                   (const RLine_2&)h1, (const RLine_2&)h2);
}


template < class R >
inline
Comparison_result
compare_x(const Line_2<R> &l,
          const Line_2<R> &h1,
          const Line_2<R> &h2)
{
  typedef typename  R::Line_2_base  RLine_2 ;
  return compare_x((const RLine_2&)l, (const RLine_2&)h1,
                   (const RLine_2&)l, (const RLine_2&)h2);
}

template < class R >
inline
Comparison_result
compare_y(const Point_2<R> &p,
          const Line_2<R> &l1,
          const Line_2<R> &l2)
{
  typedef typename  R::Point_2_base  RPoint_2;
  typedef typename  R::Line_2_base  RLine_2 ;
  return compare_y((const RPoint_2&)p,
                   (const RLine_2&)l1,
                   (const RLine_2&)l2);
}

template < class R >
inline
Comparison_result
compare_y(const Line_2<R> &l1,
          const Line_2<R> &l2,
          const Line_2<R> &h1,
          const Line_2<R> &h2)
{
  typedef typename  R::Line_2_base  RLine_2 ;
  return compare_y((const RLine_2&)l1, (const RLine_2&)l2,
                   (const RLine_2&)h1, (const RLine_2&)h2);
}

template < class R >
inline
Comparison_result
compare_y(const Line_2<R> &l,
          const Line_2<R> &h1,
          const Line_2<R> &h2)
{
  typedef typename  R::Line_2_base  RLine_2 ;
  return compare_y((const RLine_2&)l, (const RLine_2&)h1,
                   (const RLine_2&)l, (const RLine_2&)h2);
}


template < class R >
inline
Comparison_result
compare_y_at_x(const Point_2<R> &p, const Line_2<R> &h)
{
  typedef typename  R::Point_2_base  RPoint_2;
  typedef typename  R::Line_2_base  RLine_2 ;
  return compare_y_at_x((const RPoint_2&)p,
                        (const RLine_2&)h);
}

template < class R >
inline
Comparison_result
compare_y_at_x(const Line_2<R> &l1,
                    const Line_2<R> &l2,
                    const Line_2<R> &h)
{
  typedef typename  R::Line_2_base  RLine_2 ;
  return compare_y_at_x((const RLine_2&)l1,
                        (const RLine_2&)l2,
                        (const RLine_2&)h) ;
}

template < class R >
inline
Comparison_result
compare_y_at_x(const Point_2<R> &p,
               const Line_2<R> &h1,
               const Line_2<R> &h2)
{
  CGAL_kernel_precondition( (! h1.is_vertical()) && (! h2.is_vertical()) );
  typedef typename  R::Point_2_base  RPoint_2;
  typedef typename  R::Line_2_base  RLine_2 ;
  return compare_y_at_x((const RPoint_2&)p,
                             (const RLine_2&)h1,
                             (const RLine_2&)h2);
}

template < class R >
inline
Comparison_result
compare_y_at_x(const Line_2<R> &l1,
               const Line_2<R> &l2,
               const Line_2<R> &h1,
               const Line_2<R> &h2)
{
  CGAL_kernel_precondition( (! h1.is_vertical()) && (! h2.is_vertical()) );
  typedef typename  R::Line_2_base  RLine_2 ;
  return compare_y_at_x((const RLine_2&)l1,
                        (const RLine_2&)l2,
                        (const RLine_2&)h1,
                        (const RLine_2&)h2);
}

CGAL_END_NAMESPACE


#endif  // CGAL_PREDICATES_ON_LINES_2_H
