// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
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

#ifndef CGAL_CARTESIAN_POINT_CONSTRUCTIONS_3_H
#define CGAL_CARTESIAN_POINT_CONSTRUCTIONS_3_H

#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/constructions/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template < class K >
inline
typename K::Point_3
midpoint(const PointC3<K> &p,
         const PointC3<K> &q)
{
  return K().construct_midpoint_3_object()(p, q);
}

template < class K >
inline
PointC3<K>
centroid(const PointC3<K> &p,
         const PointC3<K> &q,
         const PointC3<K> &r,
         const PointC3<K> &s)
{
  return K().construct_centroid_3_object()(p, q, r, s);
}

template < class K >
inline
PointC3<K>
centroid(const PointC3<K> &p,
         const PointC3<K> &q,
         const PointC3<K> &r)
{
  return K().construct_centroid_3_object()(p, q, r);
}

template < class K >
inline
PointC3<K>
circumcenter(const PointC3<K> &p,
             const PointC3<K> &q,
             const PointC3<K> &r,
             const PointC3<K> &s)
{
  return K().construct_circumcenter_3_object()(p, q, r, s);
}

template < class K >
inline
PointC3<K>
circumcenter(const PointC3<K> &p,
             const PointC3<K> &q,
             const PointC3<K> &r)
{
  return K().construct_circumcenter_3_object()(p, q, r);
}

template <class K>
CGAL_KERNEL_LARGE_INLINE
PointC3<K>
point_on_line(int i, const LineC3<K> &l)
{
  typename K::FT x, y, z;
  point_on_lineC3(l.point().x(), l.point().y(), l.point().z(),
                  l.direction().dx(), l.direction().dy(), l.direction().dz(),
                  i, x, y, z);
  return PointC3<K>(x, y, z);
}

template <class K>
CGAL_KERNEL_LARGE_INLINE
PointC3<K>
projection_line(const PointC3<K> &p, const LineC3<K> &l)
{
  typename K::FT x, y, z;
  projection_lineC3(p.x(), p.y(), p.z(),
		    l.point().x(), l.point().y(), l.point().z(),
                    l.direction().dx(), l.direction().dy(), l.direction().dz(),
                    x, y, z);
  return PointC3<K>(x, y, z);
}

template <class K>
CGAL_KERNEL_LARGE_INLINE
PointC3<K>
point_on_plane(const PlaneC3<K> &p)
{
  typename K::FT x, y, z;
  point_on_planeC3(p.a(), p.b(), p.c(), p.d(), x, y, z);
  return PointC3<K>(x, y, z);
}

template <class K>
CGAL_KERNEL_LARGE_INLINE
PointC3<K>
projection_plane(const PointC3<K> &p,
                 const PlaneC3<K> &h)
{
  typename K::FT x, y, z;
  projection_planeC3(h.a(), h.b(), h.c(), h.d(),
                     p.x(), p.y(), p.z(),
                     x, y, z);
  return PointC3<K>(x, y, z);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_POINT_CONSTRUCTIONS_3_H
