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
  // Intersection_traits.h provides the defaults, we add overloads for all special cases

  // the special plane_3 function
  template <class K>
  inline 
  #if CGAL_INTERSECTION_VERSION < 2
  CGAL::Object
  #else
  typename boost::optional< boost::variant< Point_3<K>, Line_3<K>, Plane_3<K> > >
  #endif
  intersection(const Plane_3<K> &plane1, const Plane_3<K> &plane2,
               const Plane_3<K> &plane3)
  {
    return K().intersect_3_object()(plane1, plane2, plane3);
  }

  // intersections with Bbox_3 for which no kernel traits exist
  template<class A>
  inline 
  #if CGAL_INTERSECTION_VERSION < 2
  CGAL::Object
  #else
  typename IT< typename CGAL::Kernel_traits<A>::Kernel, A, CGAL::Bbox_3 >::result_type
  #endif
  intersection(const A& a, const CGAL::Bbox_3& b) {
    typedef typename CGAL::Kernel_traits<A>::Kernel Kernel;
    return Kernel().intersect_3_object()(a, b);
  }

  template<class A>
  inline 
  #if CGAL_INTERSECTION_VERSION < 2
  CGAL::Object
  #else
  typename IT< typename CGAL::Kernel_traits<A>::Kernel, A, CGAL::Bbox_3 >::result_type
  #endif
  intersection(const CGAL::Bbox_3& b, const A& a) {
    typedef typename CGAL::Kernel_traits<A>::Kernel Kernel;
    return Kernel().intersect_3_object()(a, b);
  }

  template<class A>
  inline bool
  do_intersect(const A& a, const CGAL::Bbox_3& b) {
    typedef typename CGAL::Kernel_traits<A>::Kernel Kernel;
    return Kernel().do_intersect_3_object()(a, b);
  }

  template<class A>
  inline bool
  do_intersect(const CGAL::Bbox_3& b, const A& a) {
    typedef typename CGAL::Kernel_traits<A>::Kernel Kernel;
    return Kernel().do_intersect_3_object()(a, b);
  }
} // CGAL


#endif /* _INTERSECTION_3_0_H_ */
