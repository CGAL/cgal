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
// release_date  : 2000, October 15
//
// source        : webS3/S3.lw
// file          : include/CGAL/SimpleCartesian/distance_predicatesS3.h
// package       : S3 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.7
// revision_date : 15 Oct 2000
// author(s)     : Stefan Schirra <Stefan.Schirra@@mpi-sb.mpg.de>
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================

#ifndef DISTANCE_PREDICATESS3_H
#define DISTANCE_PREDICATESS3_H

#include <CGAL/SimpleCartesian/PointS3.h>
#include <CGAL/SimpleCartesian/PlaneS3.h>
#include <CGAL/predicates/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
cmp_dist_to_point(const PointS3<FT>& p,
                  const PointS3<FT>& q,
                  const PointS3<FT>& r)
{
  return cmp_dist_to_pointC3(p.x(),p.y(),p.z(),
                             q.x(),q.y(),q.z(),
                             r.x(),r.y(),r.z());
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_dist_to_point(const PointS3<FT>& p,
                         const PointS3<FT>& q,
                         const PointS3<FT>& r)
{
  return has_larger_dist_to_pointC3(p.x(),p.y(),p.z(),
                                    q.x(),q.y(),q.z(),
                                    r.x(),r.y(),r.z());
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
bool
has_smaller_dist_to_point(const PointS3<FT>& p,
                          const PointS3<FT>& q,
                          const PointS3<FT>& r)
{
  return has_smaller_dist_to_pointC3(p.x(),p.y(),p.z(),
                                     q.x(),q.y(),q.z(),
                                     r.x(),r.y(),r.z());
}
template < class FT >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
cmp_signed_dist_to_plane(const PlaneS3<FT>& h,
                         const PointS3<FT>& p,
                         const PointS3<FT>& q)
{
  return cmp_signed_dist_to_directionC3(h.a(),h.b(),h.c(),
                                        p.x(),p.y(),p.z(),
                                        q,x(),q.y(),q.z());
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_signed_dist_to_plane(const PlaneS3<FT>& h,
                                const PointS3<FT>& p,
                                const PointS3<FT>& q)
{
  return has_larger_signed_dist_to_directionC3(h.a(),h.b(),h.c(),
                                               p.x(),p.y(),p.z(),
                                               q,x(),q.y(),q.z());
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
bool
has_smaller_signed_dist_to_plane(const PlaneS3<FT>& h,
                                 const PointS3<FT>& p,
                                 const PointS3<FT>& q)
{
  return has_smaller_signed_dist_to_directionC3(h.a(),h.b(),h.c(),
                                                p.x(),p.y(),p.z(),
                                                q,x(),q.y(),q.z());
}
template < class FT >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
cmp_signed_dist_to_plane(const PointS3<FT>& hp,
                         const PointS3<FT>& hq,
                         const PointS3<FT>& hr,
                         const PointS3<FT>& p,
                         const PointS3<FT>& q)
{
  return cmp_signed_dist_to_planeC3(hp.x(),hp.y(),hp.z(),
                                    hq.x(),hq.y(),hq.z(),
                                    hr.x(),hr.y(),hr.z(),
                                    p.x(),p.y(),p.z(),
                                    q,x(),q.y(),q.z());
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_signed_dist_to_plane(const PointS3<FT>& hp,
                                const PointS3<FT>& hq,
                                const PointS3<FT>& hr,
                                const PointS3<FT>& p,
                                const PointS3<FT>& q)
{
  return has_larger_signed_dist_to_planeC3(hp.x(),hp.y(),hp.z(),
                                           hq.x(),hq.y(),hq.z(),
                                           hr.x(),hr.y(),hr.z(),
                                           p.x(),p.y(),p.z(),
                                           q,x(),q.y(),q.z());
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
bool
has_smaller_signed_dist_to_plane(const PointS3<FT>& hp,
                                 const PointS3<FT>& hq,
                                 const PointS3<FT>& hr,
                                 const PointS3<FT>& p,
                                 const PointS3<FT>& q)
{
  return has_smaller_signed_dist_to_planeC3(hp.x(),hp.y(),hp.z(),
                                            hq.x(),hq.y(),hq.z(),
                                            hr.x(),hr.y(),hr.z(),
                                            p.x(),p.y(),p.z(),
                                            q,x(),q.y(),q.z());
}


CGAL_END_NAMESPACE

#endif // DISTANCE_PREDICATESS3_H
