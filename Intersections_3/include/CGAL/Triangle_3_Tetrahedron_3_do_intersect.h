// Copyright (c) 2005  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Nico Kruithof

#ifndef CGAL_TRIANGLE_3_TETRAHEDRON_3_DO_INTERSECT_H
#define CGAL_TRIANGLE_3_TETRAHEDRON_3_DO_INTERSECT_H

#include <CGAL/Triangle_3_Triangle_3_do_intersect.h>
#include <CGAL/internal/Intersections_3/Iso_cuboid_3_Triangle_3_do_intersect.h>

namespace CGAL {

  template <class K>
  class Tetrahedron_3;

namespace internal {

// This code is not optimized:
  template <class K, class Bounded>
typename K::Boolean
do_intersect_tetrahedron_bounded(const typename Bounded &tr,
                                 const typename K::Tetrahedron_3 &tet,
                                 const typename K::Point_3 p,
                                 const K & k)
{
    typedef typename K::Triangle_3 Triangle;

    CGAL_kernel_precondition( ! k.is_degenerate_3_object() (tr) );
    CGAL_kernel_precondition( ! k.is_degenerate_3_object() (tet) );

    if (do_intersect(tr, Triangle(tet[0], tet[1], tet[2]), k)) return true;
    if (do_intersect(tr, Triangle(tet[0], tet[1], tet[3]), k)) return true;
    if (do_intersect(tr, Triangle(tet[0], tet[2], tet[3]), k)) return true;
    if (do_intersect(tr, Triangle(tet[1], tet[2], tet[3]), k)) return true;

    return k.has_on_bounded_side_3_object()(tet, p);
}


template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3 &tet,
	     const typename K::Triangle_3 &tr,
	     const K & k)
{
  return do_intersect_tetrahedron_bounded(tr, tet, tr[0], k);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Triangle_3 &tr,
             const typename K::Tetrahedron_3 &tet,
	     const K & k)
{
  return do_intersect_tetrahedron_bounded(tr, tet, tr[0], k);
}

  
template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3 &tet,
	     const typename K::Segment_3 &seg,
	     const K & k)
{
  return do_intersect_tetrahedron_bounded(seg, tet, seg.source(), k);
}
    
template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Segment_3 &seg,
             const typename K::Tetrahedron_3 &tet,
	     const K & k)
{
  return do_intersect_tetrahedron_bounded(seg, tet, seg.source(), k);
}
  

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3 &tet,
             const typename K::Iso_cuboid_3 &ic,
	     const K & k)
{
  return do_intersect_tetrahedron_bounded(ic, tet, ic[0], k);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Iso_cuboid_3 &ic,
             const typename K::Tetrahedron_3 &tet,
	     const K & k)
{
  return do_intersect_tetrahedron_bounded(ic, tet, ic[0], k);
}

  
template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3 &tet,
             const typename K::Sphere_3 &sp,
	     const K & k)
{
  return do_intersect_tetrahedron_bounded(sp, tet, sp.center(), k);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Sphere_3 &sp,
             const typename K::Tetrahedron_3 &tet,
	     const K & k)
{
  return do_intersect_tetrahedron_bounded(sp, tet, sp.center(), k);
}

  
template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3 &tet,
	     const typename K::Tetrahedron_3 &sp,
	     const K & k)
{
  return do_intersect_tetrahedron_bounded(sp, tet, tet[0], k);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Tetrahedron_3 &tet,
	     const typename Bbox_3 &bb,
	     const K & k)
{
  return do_intersect_tetrahedron_bounded(bb, tet, typename K::Point_3(bb.xmin(), bb.ymin(), bb.zmin()), k);
}

  template <class K>
inline
typename K::Boolean
do_intersect(const typename Bbox_3 &bb,
             const typename K::Tetrahedron_3 &tet,
	     const K & k)
{
  return do_intersect_tetrahedron_bounded(bb, tet, typename K::Point_3(bb.xmin(), bb.ymin(), bb.zmin()), k);
}
  
} // namespace internal

CGAL_DO_INTERSECT_FUNCTION(Triangle_3, Tetrahedron_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Segment_3, Tetrahedron_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Iso_cuboid_3, Tetrahedron_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Sphere_3, Tetrahedron_3, 3)
CGAL_DO_INTERSECT_FUNCTION_SELF(Tetrahedron_3, 3)


template<typename K>
bool do_intersect(const CGAL::Tetrahedron_3<K>& a,
                  const CGAL::Bbox_3& b) {
  return K().do_intersect_3_object()(a, b);
}

  
template<typename K>
bool do_intersect(const CGAL::Bbox_3& b,
                  const CGAL::Tetrahedron_3<K>& a) {
  return K().do_intersect_3_object()(a, b);
}

} //namespace CGAL

#endif // CGAL_TRIANGLE_3_TETRAHEDRON_3_DO_INTERSECT_H
