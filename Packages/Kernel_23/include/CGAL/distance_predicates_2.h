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
// file          : distance_predicates_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_DISTANCE_PREDICATES_2_H
#define CGAL_DISTANCE_PREDICATES_2_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/distance_predicatesH2.h>
#endif

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/Cartesian/distance_predicates_2.h>
#endif

#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>

CGAL_BEGIN_NAMESPACE


template < class R >
inline
Comparison_result
cmp_dist_to_point(const Point_2<R>& p,
                  const Point_2<R>& q,
                  const Point_2<R>& r)
{
  typedef typename  R::Point_2_base  RPoint_2;
  return cmp_dist_to_point((const RPoint_2& )p,
                           (const RPoint_2& )q,
                           (const RPoint_2& )r );
}

template < class R >
inline
bool
has_larger_dist_to_point(const Point_2<R>& p,
                         const Point_2<R>& q,
                         const Point_2<R>& r)
{
  typedef typename  R::Point_2_base  RPoint_2;
  return has_larger_dist_to_point((const RPoint_2& )p,
                                  (const RPoint_2& )q,
                                  (const RPoint_2& )r );
}

template < class R >
inline
bool
has_smaller_dist_to_point(const Point_2<R>& p,
                               const Point_2<R>& q,
                               const Point_2<R>& r)
{
  typedef typename  R::Point_2_base  RPoint_2;
  return has_smaller_dist_to_point((const RPoint_2& )p,
                                        (const RPoint_2& )q,
                                        (const RPoint_2& )r );
}

template < class R >
inline
Comparison_result
cmp_signed_dist_to_line(const Line_2<R>&  l,
                        const Point_2<R>& p,
                        const Point_2<R>& q)
{
  typedef typename  R::Point_2_base  RPoint_2;
  typedef typename  R::Line_2_base  RLine_2 ;
  return cmp_signed_dist_to_line((const RLine_2& )l,
                                      (const RPoint_2& )p,
                                      (const RPoint_2& )q );
}

template < class R >
inline
bool
has_larger_signed_dist_to_line(const Line_2<R>&  l,
                                    const Point_2<R>& p,
                                    const Point_2<R>& q)
{
  typedef typename  R::Point_2_base  RPoint_2;
  typedef typename  R::Line_2_base  RLine_2 ;
  return has_larger_signed_dist_to_line((const RLine_2& )l,
                                             (const RPoint_2& )p,
                                             (const RPoint_2& )q );
}

template < class R >
inline
bool
has_smaller_signed_dist_to_line(const Line_2<R>&  l,
                                     const Point_2<R>& p,
                                     const Point_2<R>& q)
{
  typedef typename  R::Point_2_base  RPoint_2;
  typedef typename  R::Line_2_base  RLine_2 ;
  return has_smaller_signed_dist_to_line((const RLine_2& )l,
                                              (const RPoint_2& )p,
                                              (const RPoint_2& )q );
}

template < class R >
inline
Comparison_result
cmp_signed_dist_to_line(const Point_2<R>& p,
                             const Point_2<R>& q,
                             const Point_2<R>& r,
                             const Point_2<R>& s)
{
  typedef typename  R::Point_2_base  RPoint_2;
  return cmp_signed_dist_to_line((const RPoint_2& )p,
                                      (const RPoint_2& )q,
                                      (const RPoint_2& )r,
                                      (const RPoint_2& )s );
}

template < class R >
inline
bool
has_smaller_signed_dist_to_line(const Point_2<R>& p,
                                     const Point_2<R>& q,
                                     const Point_2<R>& r,
                                     const Point_2<R>& s)
{
  typedef typename  R::Point_2_base  RPoint_2;
  return has_smaller_signed_dist_to_line((const RPoint_2& )p,
                                              (const RPoint_2& )q,
                                              (const RPoint_2& )r,
                                              (const RPoint_2& )s );
}

template < class R >
inline
bool
has_larger_signed_dist_to_line(const Point_2<R>& p,
                                    const Point_2<R>& q,
                                    const Point_2<R>& r,
                                    const Point_2<R>& s)
{
  typedef typename  R::Point_2_base  RPoint_2;
  return has_larger_signed_dist_to_line((const RPoint_2& )p,
                                             (const RPoint_2& )q,
                                             (const RPoint_2& )r,
                                             (const RPoint_2& )s );
}

CGAL_END_NAMESPACE

#endif //CGAL_DISTANCE_PREDICATES_2_H
