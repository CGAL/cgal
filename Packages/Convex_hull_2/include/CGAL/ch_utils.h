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
// release_date  : 2000, August 03
//
// file          : ch_utils.h
// package       : Convex_hull (3.3)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// source        : convex_hull_2.lw
// revision      : 3.3
// revision_date : 03 Aug 2000
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CH_UTILS_H
#define CGAL_CH_UTILS_H

#ifndef CGAL_CONFIG_H
#include <CGAL/config.h>
#endif // CGAL_CONFIG_H
#include <CGAL/ch_assertions.h>

#define CGAL_CH_USE_ARGUMENT(arg)  (void)(arg)


CGAL_BEGIN_NAMESPACE
template <class Point, class BinaryPredicate>
class ch_Binary_predicate_reversor
{
public:
  ch_Binary_predicate_reversor() {}
  ch_Binary_predicate_reversor( const BinaryPredicate& p) : bp(p) {}

  bool operator() (const Point& p1, const Point& p2) const
       { return bp(p2,p1); }

private:
  BinaryPredicate  bp;
};
CGAL_END_NAMESPACE

#endif // CGAL_CH_UTILS_H

