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

#ifndef CGAL_CARTESIAN_DISTANCE_PREDICATES_2_H
#define CGAL_CARTESIAN_DISTANCE_PREDICATES_2_H

#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

template <class K>
inline
bool
has_larger_distance_to_point(const PointC2<K>& p,
                             const PointC2<K>& q,
                             const PointC2<K>& r)
{
  return has_larger_dist_to_pointC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

template <class K>
inline
Comparison_result
compare_signed_distance_to_line(const LineC2<K>&  l,
                                const PointC2<K>& p,
                                const PointC2<K>& q)
{
  return cmp_signed_dist_to_directionC2(l.a(), l.b(), p.x(), p.y(),
                                        q.x(), q.y());
}

template <class K>
inline
bool
has_larger_signed_distance_to_line(const LineC2<K>&  l,
                                   const PointC2<K>& p,
                                   const PointC2<K>& q)
{
  return has_larger_signed_dist_to_directionC2(l.a(), l.b(), p.x(), p.y(),
                                               q.x(), q.y());
}

template <class K>
inline
bool
has_smaller_signed_distance_to_line(const LineC2<K>&  l,
                                    const PointC2<K>& p,
                                    const PointC2<K>& q)
{
  return has_smaller_signed_dist_to_directionC2(l.a(), l.b(), p.x(), p.y(),
                                                q.x(), q.y());
}

template <class K>
inline
Comparison_result
compare_signed_distance_to_line(const PointC2<K>& p,
                                const PointC2<K>& q,
                                const PointC2<K>& r,
                                const PointC2<K>& s)
{
  return cmp_signed_dist_to_lineC2(p.x(), p.y(), q.x(), q.y(),
                                   r.x(), r.y(), s.x(), s.y());
}

template <class K>
inline
bool
has_smaller_signed_distance_to_line(const PointC2<K>& p,
                                    const PointC2<K>& q,
                                    const PointC2<K>& r,
                                    const PointC2<K>& s)
{
  return has_smaller_signed_dist_to_lineC2(p.x(), p.y(), q.x(), q.y(),
                                           r.x(), r.y(), s.x(), s.y());
}

template <class K>
inline
bool
has_larger_signed_distance_to_line(const PointC2<K>& p,
                                   const PointC2<K>& q,
                                   const PointC2<K>& r,
                                   const PointC2<K>& s)
{
  return has_larger_signed_dist_to_lineC2(p.x(), p.y(), q.x(), q.y(),
                                          r.x(), r.y(), s.x(), s.y());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_DISTANCE_PREDICATES_2_H
