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
// file          : include/CGAL/Cartesian/constructions_on_planes_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_CONSTRUCTIONS_ON_PLANES_D_H
#define CGAL_CARTESIAN_CONSTRUCTIONS_ON_PLANES_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/Cartesian/Point_d.h>
#include <CGAL/Cartesian/Plane_d.h>
#include <CGAL/constructions/kernel_ftCd.h>
#include <CGAL/predicates/kernel_ftCd.h>

CGAL_BEGIN_NAMESPACE

template < class R, class PointIterator >
CGAL_KERNEL_LARGE_INLINE
PlaneCd<R CGAL_CTAG>
plane_from_points(int dim,
                  const PointIterator &first, const PointIterator &last,
		  const R &r = R())
{
  CGAL_kernel_precondition(last-first == dim);
  PointIterator i;
  for (i=first; i!=last; ++i)
    CGAL_kernel_precondition( (*i).dimension() == dim );
  PlaneCd<R CGAL_CTAG> h(dim);
  plane_from_pointsCd(dim, first, last, h.begin(), r);
  return h;
}

template < class R, class PointIterator >
CGAL_KERNEL_LARGE_INLINE
PlaneCd<R CGAL_CTAG>
plane_from_points(int dim,
                  const PointIterator &first, const PointIterator &last,
		  const PlaneCd<R CGAL_CTAG>::Point_d &p, Oriented_side o,
		  const R &r = R())
{
  CGAL_kernel_precondition( o != ON_ORIENTED_BOUNDARY );
  PlaneCd<R CGAL_CTAG> h = plane_from_points(dim,first,last,r);
  if (side_of_oriented_plane(h,p) == o) return h;
  return h.opposite();
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PlaneCd<R CGAL_CTAG>
plane_from_point_direction(const PointCd<R CGAL_CTAG>& p,
                           const DirectionCd<R CGAL_CTAG>& d)
{
  CGAL_kernel_precondition( p.dimension() == d.dimension() );
  PlaneCd<R CGAL_CTAG> h(p.dimension());
  plane_from_point_directionCd(p.dimension(), p.begin(), p.end(),
                               d.begin(), d.end(), h.begin(), R());
  return h; 
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointCd<R CGAL_CTAG>
point_on_plane(const PlaneCd<R CGAL_CTAG>& h, int i = 0)
{
  PointCd<R CGAL_CTAG> p(h.dimension());
  point_on_planeCd(h.dimension(), h.begin(), h.end(), i, p.begin(), R());
  return p;
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointCd<R CGAL_CTAG>
projection_plane(const PointCd<R CGAL_CTAG>& p,
                 const PlaneCd<R CGAL_CTAG>& h)
{
  CGAL_kernel_precondition( p.dimension() == h.dimension() );
  PointCd<R CGAL_CTAG> q(h.dimension());
  projection_planeCd(h.dimension(), h.begin(), h.end(),
                     p.begin(), p.end(), q.begin(), R());
  return q;
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CONSTRUCTIONS_ON_PLANES_D_H
