// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : include/CGAL/ch_predicate_classes_3.h
// package       : $CGAL_Package: Convex_hull_3 $
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
// chapter       : Convex Hulls and Extreme Points
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: predicates used for 3D quickhull and convexity checking
// ============================================================================

#ifndef CGAL_PREDICATES_CLASSES_3_H
#define CGAL_PREDICATES_CLASSES_3_H


#include <CGAL/distance_predicates_3.h>

CGAL_BEGIN_NAMESPACE

template <class Plane, class Point>
class Less_signed_dist_to_plane_3
{
 public:
  Less_signed_dist_to_plane_3( const Plane& p) : _p(p) {}
  bool operator()( const Point& q, const Point& r)
       { return has_smaller_signed_dist_to_plane( _p,q,r); }
 private:
  Plane _p;
};


CGAL_END_NAMESPACE


#endif // CGAL_PREDICATES_CLASSES_3_H
