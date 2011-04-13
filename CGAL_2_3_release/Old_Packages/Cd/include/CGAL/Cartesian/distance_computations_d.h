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
// file          : include/CGAL/Cartesian/distance_computations_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_DISTANCE_COMPUTATIONS_D_H
#define CGAL_CARTESIAN_DISTANCE_COMPUTATIONS_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/Cartesian/Point_d.h>
#include <CGAL/Cartesian/Plane_d.h>
#include <CGAL/constructions/kernel_ftCd.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
typename R::FT
squared_distance(const PointCd<R CGAL_CTAG> &p,
                 const PointCd<R CGAL_CTAG> &q)
{
  CGAL_kernel_precondition( p.dimension() == q.dimension() );
  return squared_distanceCd(p.begin(),p.end(),q.begin());
}

template < class R >
inline
typename R::FT
scaled_distance_to_plane(const PlaneCd<R CGAL_CTAG> &h,
                         const PointCd<R CGAL_CTAG> &p)
{
  return scaled_distance_to_planeCd(h.begin(),h.end(),p.begin(),p.end());
}

template < class R, class PointIterator >
inline
typename R::FT
scaled_distance_to_plane(const PointIterator &first,
                         const PointIterator &last,
                         const PointCd<R CGAL_CTAG> &p)
{
  return scaled_distance_to_planeCd(first,last,p.begin(),p.end(),R());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_DISTANCE_COMPUTATIONS_D_H
