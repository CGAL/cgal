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
// release_date  : 
//
// file          : include/CGAL/ch_utils.h
// package       : Convex_hull_2 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CH_UTILS_H
#define CGAL_CH_UTILS_H

#include <CGAL/config.h>
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

