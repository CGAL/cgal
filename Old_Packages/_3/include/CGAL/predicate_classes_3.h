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
// file          : predicate_classes_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_PREDICATES_CLASSES_3_H
#define CGAL_PREDICATES_CLASSES_3_H

CGAL_BEGIN_NAMESPACE

/*
template < class Point>
class Less_xyz
{
public:
  bool
  operator()( const Point& p, const Point& q)
  { return lexicographically_xyz_smaller(p,q); }
};

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
*/

CGAL_END_NAMESPACE

#endif // CGAL_PREDICATES_CLASSES_3_H
