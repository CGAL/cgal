// ==========================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// --------------------------------------------------------------------------
//

// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/distance_predicates_3.h
// source        : include/CGAL/Cartesian/distance_predicates_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================


#ifndef DISTANCE_PREDICATES_3_H
#define DISTANCE_PREDICATES_3_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#include <CGAL/Cartesian/redefine_names_3.h>
#endif

#ifndef CGAL_CARTESIAN_POINT_3_H
#include <CGAL/Cartesian/Point_3.h>
#endif
#ifndef CGAL_CARTESIAN_PLANE_3_H
#include <CGAL/Cartesian/Plane_3.h>
#endif
#ifndef CGAL_PREDICATES_KERNEL_FTC3_H
#include <CGAL/predicates/kernel_ftC3.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
cmp_dist_to_point(const PointC3<R CGAL_CTAG> &p,
                  const PointC3<R CGAL_CTAG> &q,
                  const PointC3<R CGAL_CTAG> &r)
{
  return cmp_dist_to_pointC3(p.x(),p.y(),p.z(),
                             q.x(),q.y(),q.z(),
                             r.x(),r.y(),r.z());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_dist_to_point(const PointC3<R CGAL_CTAG> &p,
                         const PointC3<R CGAL_CTAG> &q,
                         const PointC3<R CGAL_CTAG> &r)
{
  return has_larger_dist_to_pointC3(p.x(),p.y(),p.z(),
                                    q.x(),q.y(),q.z(),
                                    r.x(),r.y(),r.z());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
has_smaller_dist_to_point(const PointC3<R CGAL_CTAG> &p,
                          const PointC3<R CGAL_CTAG> &q,
                          const PointC3<R CGAL_CTAG> &r)
{
  return has_smaller_dist_to_pointC3(p.x(),p.y(),p.z(),
                                     q.x(),q.y(),q.z(),
                                     r.x(),r.y(),r.z());
}
template < class R >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
cmp_signed_dist_to_plane(const PlaneC3<R CGAL_CTAG> &h,
                         const PointC3<R CGAL_CTAG> &p,
                         const PointC3<R CGAL_CTAG> &q)
{
  return cmp_signed_dist_to_directionC3(h.a(),h.b(),h.c(),
                                        p.x(),p.y(),p.z(),
                                        q,x(),q.y(),q.z());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_signed_dist_to_plane(const PlaneC3<R CGAL_CTAG> &h,
                                const PointC3<R CGAL_CTAG> &p,
                                const PointC3<R CGAL_CTAG> &q)
{
  return has_larger_signed_dist_to_directionC3(h.a(),h.b(),h.c(),
                                               p.x(),p.y(),p.z(),
                                               q,x(),q.y(),q.z());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
has_smaller_signed_dist_to_plane(const PlaneC3<R CGAL_CTAG> &h,
                                 const PointC3<R CGAL_CTAG> &p,
                                 const PointC3<R CGAL_CTAG> &q)
{
  return has_smaller_signed_dist_to_directionC3(h.a(),h.b(),h.c(),
                                                p.x(),p.y(),p.z(),
                                                q,x(),q.y(),q.z());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
cmp_signed_dist_to_plane(const PointC3<R CGAL_CTAG> &hp,
                         const PointC3<R CGAL_CTAG> &hq,
                         const PointC3<R CGAL_CTAG> &hr,
                         const PointC3<R CGAL_CTAG> &p,
                         const PointC3<R CGAL_CTAG> &q)
{
  return cmp_signed_dist_to_planeC3(hp.x(),hp.y(),hp.z(),
                                    hq.x(),hq.y(),hq.z(),
                                    hr.x(),hr.y(),hr.z(),
                                    p.x(),p.y(),p.z(),
                                    q,x(),q.y(),q.z());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_signed_dist_to_plane(const PointC3<R CGAL_CTAG> &hp,
                                const PointC3<R CGAL_CTAG> &hq,
                                const PointC3<R CGAL_CTAG> &hr,
                                const PointC3<R CGAL_CTAG> &p,
                                const PointC3<R CGAL_CTAG> &q)
{
  return has_larger_signed_dist_to_planeC3(hp.x(),hp.y(),hp.z(),
                                           hq.x(),hq.y(),hq.z(),
                                           hr.x(),hr.y(),hr.z(),
                                           p.x(),p.y(),p.z(),
                                           q,x(),q.y(),q.z());
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
has_smaller_signed_dist_to_plane(const PointC3<R CGAL_CTAG> &hp,
                                 const PointC3<R CGAL_CTAG> &hq,
                                 const PointC3<R CGAL_CTAG> &hr,
                                 const PointC3<R CGAL_CTAG> &p,
                                 const PointC3<R CGAL_CTAG> &q)
{
  return has_smaller_signed_dist_to_planeC3(hp.x(),hp.y(),hp.z(),
                                            hq.x(),hq.y(),hq.z(),
                                            hr.x(),hr.y(),hr.z(),
                                            p.x(),p.y(),p.z(),
                                            q,x(),q.y(),q.z());
}

CGAL_END_NAMESPACE

#endif // DISTANCE_PREDICATES_3_H
