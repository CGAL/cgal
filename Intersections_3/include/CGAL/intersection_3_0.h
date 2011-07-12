// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Philipp Möller

#ifndef CGAL_INTERSECTION_3_0_H
#define CGAL_INTERSECTION_3_0_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Kernel_traits.h>

namespace CGAL {

  template <class A, class B>
  inline
  Object
  intersection(const A& a,
               const B& b, typename Intersection_traits_3< typename CGAL::Kernel_traits<A>::Kernel, A, B>::result_type* d = 0)
  {
    (void)d;
    typedef typename CGAL::Kernel_traits<A>::Kernel::Intersect_3 Intersect;
    return Intersect()(a, b);
  }

  // the special plane_3 function
  template <class K>
  inline
  Object 
  intersection(const Plane_3<K> &plane1, const Plane_3<K> &plane2,
               const Plane_3<K> &plane3)
  {
    return typename K::Intersect_3()(plane1, plane2, plane3);
  }

  // intersections with Bbox_3 for which no kernel traits exist
  template<class A>
  inline Object
  intersection(const A& a, const CGAL::Bbox_3& b) {
    typedef typename CGAL::Kernel_traits<A>::Kernel::Intersect_3 Intersect;
    return Intersect()(a, b);
  }

  template<class A>
  inline Object
  intersection(const CGAL::Bbox_3& b, const A& a) {
    typedef typename CGAL::Kernel_traits<A>::Kernel::Intersect_3 Intersect;
    return Intersect()(a, b);
  }

  template <class A, class B>
  inline bool 
  do_intersect(const A& p1, const B& p2, 
               typename Intersection_traits_3< typename CGAL::Kernel_traits<A>::Kernel, A, B>::result_type* d = 0)
  {
    (void)d;
    typedef typename CGAL::Kernel_traits<A>::Kernel::Do_intersect_3 Do_intersect;
    return Do_intersect()(p1, p2);
  }

  template<class A>
  inline bool
  do_intersect(const A& a, const CGAL::Bbox_3& b) {
    typedef typename CGAL::Kernel_traits<A>::Kernel::Do_intersect_3 Do_intersect;
    return Do_intersect()(a, b);
  }

  template<class A>
  inline bool
  do_intersect(const CGAL::Bbox_3& b, const A& a) {
    typedef typename CGAL::Kernel_traits<A>::Kernel::Do_intersect_3 Do_intersect;
    return Do_intersect()(a, b);
  }
} // CGAL


#endif /* _INTERSECTION_3_0_H_ */
