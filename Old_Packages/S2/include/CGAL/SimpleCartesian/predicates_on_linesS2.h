// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 16
//
// source        : webS2/S2.lw
// file          : include/CGAL/SimpleCartesian/predicates_on_linesS2.h
// package       : S2 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.6
// revision_date : 27 Jun 2000
// author(s)     : Stefan Schirra
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================


#ifndef CGAL_PREDICATES_ON_LINESS2_H
#define CGAL_PREDICATES_ON_LINESS2_H

#include <CGAL/SimpleCartesian/simple_cartesian_classes.h>
#include <CGAL/SimpleCartesian/LineS2.h>
#include <CGAL/SimpleCartesian/predicates_on_pointsS2.h>

CGAL_BEGIN_NAMESPACE


template < class FT >
Comparison_result compare_x(const PointS2<FT> &p,
                            const LineS2<FT> &l1,
                            const LineS2<FT> &l2)
{
  return compare_xC2(p.x(), l1.a(),l1.b(),l1.c(),l2.a(),l2.b(),l2.c());
}


template < class FT >
Comparison_result compare_x(const LineS2<FT> &l1,
                            const LineS2<FT> &l2,
                            const LineS2<FT> &h1,
                            const LineS2<FT> &h2)
{
  return compare_xC2(l1.a(),l1.b(),l1.c(),l2.a(),l2.b(),l2.c(),
                     h1.a(),h1.b(),h1.c(),h2.a(),h2.b(),h2.c());
}



template < class FT >
Comparison_result compare_y(const PointS2<FT> &p,
                            const LineS2<FT> &l1,
                            const LineS2<FT> &l2)
{
  return compare_xC2(p.y(), l1.b(),l1.a(),l1.c(),l2.b(),l2.a(),l2.c());
}


template < class FT >
Comparison_result compare_y(const LineS2<FT> &l1,
                            const LineS2<FT> &l2,
                            const LineS2<FT> &h1,
                            const LineS2<FT> &h2)
{
  return compare_xC2(l1.b(),l1.a(),l1.c(),l2.b(),l2.a(),l2.c(),
                     h1.b(),h1.a(),h1.c(),h2.b(),h2.a(),h2.c());
}
template < class FT >
Comparison_result compare_y_at_x(const PointS2<FT> &p, const LineS2<FT> &h)
{
  return compare_y_at_xC2(p.x(),p.y(),h.a(),h.b(),h.c());
}


template < class FT >
Comparison_result compare_y_at_x(const PointS2<FT> &p,
                                 const LineS2<FT> &h1,
                                 const LineS2<FT> &h2)
{
  return compare_y_at_xC2(p.x(), h1.a(),h1.b(),h1.c(),h2.a(),h2.b(),h2.c());
}


template < class FT >
Comparison_result compare_y_at_x(const LineS2<FT> &l1,
                                 const LineS2<FT> &l2,
                                 const LineS2<FT> &h)
{
  return compare_y_at_xC2(l1.a(),l1.b(),l1.c(),l2.a(),l2.b(),l2.c(),
                          h.a(),h.b(),h.c());
}


template < class FT >
Comparison_result compare_y_at_x(const LineS2<FT> &l1,
                                 const LineS2<FT> &l2,
                                 const LineS2<FT> &h1,
                                 const LineS2<FT> &h2)
{
  return compare_y_at_xC2(l1.a(),l1.b(),l1.c(),l2.a(),l2.b(),l2.c(),
                          h1.a(),h1.b(),h1.c(),h2.a(),h2.b(),h2.c());
}


CGAL_END_NAMESPACE

#endif  // CGAL_PREDICATES_ON_LINESS2_H
