// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        : predicate_classes_3.fw
// file          : predicate_classes_3.h
// revision      : 2.4
// revision_date : 24 Aug 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#ifndef CGAL_PREDICATES_CLASSES_3_H
#define CGAL_PREDICATES_CLASSES_3_H
#include <CGAL/predicates_on_points_3.h>

CGAL_BEGIN_NAMESPACE

template < class Point>
class Less_xyz
{
public:
  bool
  operator()( const Point& p, const Point& q)
  { return lexicographically_xyz_smaller(p,q); }
};

CGAL_END_NAMESPACE


#endif // CGAL_PREDICATES_CLASSES_3_H
