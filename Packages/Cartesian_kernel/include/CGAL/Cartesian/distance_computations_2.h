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
// file          : include/CGAL/Cartesian/distance_computations_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_DISTANCE_COMPUTATIONS_2_H
#define CGAL_CARTESIAN_DISTANCE_COMPUTATIONS_2_H

#include <CGAL/Cartesian/redefine_names_2.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
typename R::FT
squared_distance( PointC2<R CGAL_CTAG> const& p,
                  PointC2<R CGAL_CTAG> const& q)
{
  return squared_distanceC2(p.x(),p.y(),q.x(),q.y());
}

template < class R >
inline
typename R::FT
scaled_distance_to_line( LineC2<R CGAL_CTAG> const& l,
                         PointC2<R CGAL_CTAG> const& p)
{
  return scaled_distance_to_lineC2(l.a(),l.b(),l.c(),p.x(),p.y());
}

template < class R >
inline
typename R::FT
scaled_distance_to_line( PointC2<R CGAL_CTAG> const& p,
                         PointC2<R CGAL_CTAG> const& q,
                         PointC2<R CGAL_CTAG> const& r)
{
  return scaled_distance_to_lineC2(p.x(),p.y(),q.x(),q.y(),r.x(),r.y());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_DISTANCE_COMPUTATIONS_2_H
