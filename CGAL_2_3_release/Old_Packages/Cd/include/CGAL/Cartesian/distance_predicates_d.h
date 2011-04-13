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
// file          : include/CGAL/Cartesian/distance_predicates_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_DISTANCE_PREDICATES_D_H
#define CGAL_CARTESIAN_DISTANCE_PREDICATES_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/Cartesian/Point_d.h>
#include <CGAL/Cartesian/Plane_d.h>
#include <CGAL/predicates/kernel_ftCd.h>

CGAL_BEGIN_NAMESPACE

template < class R >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
cmp_dist_to_point(const PointCd<R CGAL_CTAG> &p,
                  const PointCd<R CGAL_CTAG> &q,
                  const PointCd<R CGAL_CTAG> &r)
{
  return cmp_dist_to_pointCd(p.begin(),p.end(),
                             q.begin(),q.end(),
                             r.begin(),r.end());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_dist_to_point(const PointCd<R CGAL_CTAG> &p,
                         const PointCd<R CGAL_CTAG> &q,
                         const PointCd<R CGAL_CTAG> &r)
{
  return has_larger_dist_to_pointCd(p.begin(),p.end(),
                                    q.begin(),q.end(),
                                    r.begin(),r.end());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
has_smaller_dist_to_point(const PointCd<R CGAL_CTAG> &p,
                          const PointCd<R CGAL_CTAG> &q,
                          const PointCd<R CGAL_CTAG> &r)
{
  return has_smaller_dist_to_pointCd(p.begin(),p.end(),
                                     q.begin(),q.end(),
                                     r.begin(),r.end());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
cmp_signed_dist_to_plane(const PlaneCd<R CGAL_CTAG> &h,
                         const PointCd<R CGAL_CTAG> &p,
                         const PointCd<R CGAL_CTAG> &q)
{
  return cmp_signed_dist_to_directionCd(h.begin(),h.end(),
                                        p.begin(),p.end(),
                                        q.begin(),q.end());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_signed_dist_to_plane(const PlaneCd<R CGAL_CTAG> &h,
                                const PointCd<R CGAL_CTAG> &p,
                                const PointCd<R CGAL_CTAG> &q)
{
  return has_larger_signed_dist_to_directionCd(h.begin(),h.end(),
                                               p.begin(),p.end(),
                                               q.begin(),q.end());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
has_smaller_signed_dist_to_plane(const PlaneCd<R CGAL_CTAG> &h,
                                 const PointCd<R CGAL_CTAG> &p,
                                 const PointCd<R CGAL_CTAG> &q)
{
  return has_smaller_signed_dist_to_directionCd(h.begin(),h.end(),
                                                p.begin(),p.end(),
                                                q.begin(),q.end());
}

template < class R, class ForwardIterator >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
cmp_signed_dist_to_plane(const ForwardIterator &first,
                         const ForwardIterator &last, 
                         const PointCd<R CGAL_CTAG> &p,
                         const PointCd<R CGAL_CTAG> &q)
{
  return cmp_signed_dist_to_planeCd(first, last,
                                    p.begin(),p.end(),
                                    q.begin(),q.end());
}

template < class R, class ForwardIterator >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_signed_dist_to_plane(const ForwardIterator &first,
                                const ForwardIterator &last, 
                                const PointCd<R CGAL_CTAG> &p,
                                const PointCd<R CGAL_CTAG> &q)
{
  return has_larger_signed_dist_to_planeCd(first, last,
                                           p.begin(),p.end(),
                                           q.begin(),q.end());
}

template < class R, class ForwardIterator >
CGAL_KERNEL_LARGE_INLINE
bool
has_smaller_signed_dist_to_plane(const ForwardIterator &first,
                                 const ForwardIterator &last, 
                                 const PointCd<R CGAL_CTAG> &p,
                                 const PointCd<R CGAL_CTAG> &q)
{
  return has_smaller_signed_dist_to_planeCd(first, last,
                                            p.begin(),p.end(),
                                            q.begin(),q.end());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_DISTANCE_PREDICATES_3_H
