// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_DISTANCE_PREDICATES_3_H
#define CGAL_CARTESIAN_DISTANCE_PREDICATES_3_H

#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/Cartesian/Plane_3.h>
#include <CGAL/predicates/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template < class K >
inline
Comparison_result
compare_distance_to_point(const PointC3<K> &p,
                          const PointC3<K> &q,
                          const PointC3<K> &r)
{
  return K().compare_distance_3_object()(p, q, r);
}

template < class K >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_distance_to_point(const PointC3<K> &p,
                             const PointC3<K> &q,
                             const PointC3<K> &r)
{
  return has_larger_dist_to_pointC3(p.x(), p.y(), p.z(),
                                    q.x(), q.y(), q.z(),
                                    r.x(), r.y(), r.z());
}

template < class K >
inline
bool
has_smaller_distance_to_point(const PointC3<K> &p,
                              const PointC3<K> &q,
                              const PointC3<K> &r)
{
  return K().less_distance_to_point_3_object()(p, q, r);
}

template < class K >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
compare_signed_distance_to_plane(const PlaneC3<K> &h,
                                 const PointC3<K> &p,
                                 const PointC3<K> &q)
{ // FIXME : probably not compiled by the test-suite.
  return cmp_signed_dist_to_directionC3(h.a(), h.b(), h.c(),
                                        p.x(), p.y(), p.z(),
                                        q.x(), q.y(), q.z());
}

template < class K >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_signed_distance_to_plane(const PlaneC3<K> &h,
                                    const PointC3<K> &p,
                                    const PointC3<K> &q)
{ // FIXME : probably not compiled by the test-suite.
  return has_larger_signed_dist_to_directionC3(h.a(), h.b(), h.c(),
                                               p.x(), p.y(), p.z(),
                                               q.x(), q.y(), q.z());
}

template < class K >
inline
bool
has_smaller_signed_distance_to_plane(const PlaneC3<K> &h,
                                     const PointC3<K> &p,
                                     const PointC3<K> &q)
{ 
  return K().less_signed_distance_to_plane_3_object()(h, p, q);
}

template < class K >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
compare_signed_distance_to_plane(const PointC3<K> &hp,
                                 const PointC3<K> &hq,
                                 const PointC3<K> &hr,
                                 const PointC3<K> &p,
                                 const PointC3<K> &q)
{ // FIXME : probably not compiled by the test-suite.
  return cmp_signed_dist_to_planeC3(hp.x(), hp.y(), hp.z(),
                                    hq.x(), hq.y(), hq.z(),
                                    hr.x(), hr.y(), hr.z(),
                                    p.x(),  p.y(),  p.z(),
                                    q.x(),  q.y(),  q.z());
}

template < class K >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_signed_distance_to_plane(const PointC3<K> &hp,
                                    const PointC3<K> &hq,
                                    const PointC3<K> &hr,
                                    const PointC3<K> &p,
                                    const PointC3<K> &q)
{ // FIXME : probably not compiled by the test-suite.
  return has_larger_signed_dist_to_planeC3(hp.x(), hp.y(), hp.z(),
                                           hq.x(), hq.y(), hq.z(),
                                           hr.x(), hr.y(), hr.z(),
                                           p.x(),  p.y(),  p.z(),
                                           q.x(),  q.y(),  q.z());
}

template < class K >
CGAL_KERNEL_LARGE_INLINE
bool
has_smaller_signed_distance_to_plane(const PointC3<K> &hp,
                                     const PointC3<K> &hq,
                                     const PointC3<K> &hr,
                                     const PointC3<K> &p,
                                     const PointC3<K> &q)
{ // FIXME : probably not compiled by the test-suite.
  return has_smaller_signed_dist_to_planeC3(hp.x(), hp.y(), hp.z(),
                                            hq.x(), hq.y(), hq.z(),
                                            hr.x(), hr.y(), hr.z(),
                                            p.x(),  p.y(),  p.z(),
                                            q.x(),  q.y(),  q.z());
}

#ifndef CGAL_NO_DEPRECATED_CODE
template < class K >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
cmp_dist_to_point(const PointC3<K> &p,
                  const PointC3<K> &q,
                  const PointC3<K> &r)
{
  return cmp_dist_to_pointC3(p.x(), p.y(), p.z(),
                             q.x(), q.y(), q.z(),
                             r.x(), r.y(), r.z());
}

template < class K >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_dist_to_point(const PointC3<K> &p,
                         const PointC3<K> &q,
                         const PointC3<K> &r)
{
  return has_larger_dist_to_pointC3(p.x(), p.y(), p.z(),
                                    q.x(), q.y(), q.z(),
                                    r.x(), r.y(), r.z());
}

template < class K >
CGAL_KERNEL_LARGE_INLINE
bool
has_smaller_dist_to_point(const PointC3<K> &p,
                          const PointC3<K> &q,
                          const PointC3<K> &r)
{
  return has_smaller_dist_to_pointC3(p.x(), p.y(), p.z(),
                                     q.x(), q.y(), q.z(),
                                     r.x(), r.y(), r.z());
}

template < class K >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
cmp_signed_dist_to_plane(const PlaneC3<K> &h,
                         const PointC3<K> &p,
                         const PointC3<K> &q)
{ // FIXME : probably not compiled by the test-suite.
  return cmp_signed_dist_to_directionC3(h.a(), h.b(), h.c(),
                                        p.x(), p.y(), p.z(),
                                        q.x(), q.y(), q.z());
}

template < class K >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_signed_dist_to_plane(const PlaneC3<K> &h,
                                const PointC3<K> &p,
                                const PointC3<K> &q)
{ // FIXME : probably not compiled by the test-suite.
  return has_larger_signed_dist_to_directionC3(h.a(), h.b(), h.c(),
                                               p.x(), p.y(), p.z(),
                                               q.x(), q.y(), q.z());
}

template < class K >
CGAL_KERNEL_LARGE_INLINE
bool
has_smaller_signed_dist_to_plane(const PlaneC3<K> &h,
                                 const PointC3<K> &p,
                                 const PointC3<K> &q)
{ // FIXME : probably not compiled by the test-suite.
  return has_smaller_signed_dist_to_directionC3(h.a(), h.b(), h.c(),
                                                p.x(), p.y(), p.z(),
                                                q.x(), q.y(), q.z());
}

template < class K >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
cmp_signed_dist_to_plane(const PointC3<K> &hp,
                         const PointC3<K> &hq,
                         const PointC3<K> &hr,
                         const PointC3<K> &p,
                         const PointC3<K> &q)
{ // FIXME : probably not compiled by the test-suite.
  return cmp_signed_dist_to_planeC3(hp.x(), hp.y(), hp.z(),
                                    hq.x(), hq.y(), hq.z(),
                                    hr.x(), hr.y(), hr.z(),
                                    p.x(),  p.y(),  p.z(),
                                    q.x(),  q.y(),  q.z());
}

template < class K >
CGAL_KERNEL_LARGE_INLINE
bool
has_larger_signed_dist_to_plane(const PointC3<K> &hp,
                                const PointC3<K> &hq,
                                const PointC3<K> &hr,
                                const PointC3<K> &p,
                                const PointC3<K> &q)
{ // FIXME : probably not compiled by the test-suite.
  return has_larger_signed_dist_to_planeC3(hp.x(), hp.y(), hp.z(),
                                           hq.x(), hq.y(), hq.z(),
                                           hr.x(), hr.y(), hr.z(),
                                           p.x(),  p.y(),  p.z(),
                                           q.x(),  q.y(),  q.z());
}

template < class K >
CGAL_KERNEL_LARGE_INLINE
bool
has_smaller_signed_dist_to_plane(const PointC3<K> &hp,
                                 const PointC3<K> &hq,
                                 const PointC3<K> &hr,
                                 const PointC3<K> &p,
                                 const PointC3<K> &q)
{ // FIXME : probably not compiled by the test-suite.
  return has_smaller_signed_dist_to_planeC3(hp.x(), hp.y(), hp.z(),
                                            hq.x(), hq.y(), hq.z(),
                                            hr.x(), hr.y(), hr.z(),
                                            p.x(),  p.y(),  p.z(),
                                            q.x(),  q.y(),  q.z());
}

#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_DISTANCE_PREDICATES_3_H
