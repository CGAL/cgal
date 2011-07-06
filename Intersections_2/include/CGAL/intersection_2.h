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
// $URL$
// $Id$
// 
//
// Author(s)     : Geert-Jan Giezeman



#ifndef CGAL_INTERSECTION_2_H
#define CGAL_INTERSECTION_2_H

#include <CGAL/Intersection_traits_2.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/intersection_2_1.h>
#include <CGAL/intersection_2_2.h>
#include <CGAL/intersection_2_3.h>

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
  typedef typename CGAL::Kernel_traits<A>::Kernel K;
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(a, b);
}

template <class A, class B>
inline bool 
do_intersect(const A& p1, const B& p2, 
             typename Intersection_traits_2< typename CGAL::Kernel_traits<A>::Kernel, A, B>::result_type* d = 0)
{
  (void)d;
  typedef typename CGAL::Kernel_traits<A>::Kernel K;
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(p1, p2);
}

}



#endif // CGAL_INTERSECTION_2_H
