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

#ifndef CGAL_INTERSECTION_2_0_H
#define CGAL_INTERSECTION_2_0_H

#include <CGAL/Object.h>
#include <CGAL/Intersection_traits_2.h>
#include <CGAL/Kernel_traits.h>

namespace CGAL {
  
  // The global functions do_intersect and intersection
  
  // SFINAE is used to distinguish between a 2-dimensional and a
  // 3-dimensional object using the Intersection_traits_2/3
  template <class A, class B>
  inline
  Object
  intersection(const A& a,
               const B& b, typename Intersection_traits_2< typename CGAL::Kernel_traits<A>::Kernel, A, B>::result_type* d = 0)
  {
    (void)d;
    typedef typename CGAL::Kernel_traits<A>::Kernel::Intersect_2 Intersect;
    return Intersect()(a, b);
  }

  template <class A, class B>
  inline bool 
  do_intersect(const A& p1, const B& p2, 
               typename Intersection_traits_2< typename CGAL::Kernel_traits<A>::Kernel, A, B>::result_type* d = 0)
  {
    (void)d;
    typedef typename CGAL::Kernel_traits<A>::Kernel::Do_intersect_2 Do_intersect;
    return Do_intersect()(p1, p2);
  }

}

#endif /* _INTERSECTION_2_0_H_ */
