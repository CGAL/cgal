// Copyright (c) 2003  Utrecht University (The Netherlands),
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
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_KERNEL_GLOBAL_FUNCTIONS_2_H
#define CGAL_KERNEL_GLOBAL_FUNCTIONS_2_H

// Generic functions calling the kernel functor.

#include <CGAL/user_classes.h>

CGAL_BEGIN_NAMESPACE

template <typename K>
inline
typename K::Line_2
bisector(const typename K::Point_2 &p,
         const typename K::Point_2 &q, const K &k)
{
  return k.construct_bisector_2_object()(p, q);
}

template <typename K>
inline
typename K::Line_2
bisector(const typename K::Line_2 &l1,
         const typename K::Line_2 &l2, const K &k)
{
  return k.construct_bisector_2_object()(l1, l2);
}

template <typename K>
inline
typename K::Line_2
bisector(const Point_2<K> &p, const Point_2<K> &q)
{
  return bisector(p, q, K());
}

template <typename K>
inline
typename K::Line_2
bisector(const Line_2<K> &l1, const Line_2<K> &l2)
{
  return bisector(l1, l2, K());
}

template <typename K>
inline
bool
parallel(const typename CGAL_WRAP(K)::Line_2 &l1,
         const typename CGAL_WRAP(K)::Line_2 &l2, const K &k)
{
  return k.are_parallel_2_object()(l1, l2);
}

template <typename K>
inline
bool
parallel(const typename CGAL_WRAP(K)::Ray_2 &r1,
         const typename CGAL_WRAP(K)::Ray_2 &r2, const K &k)
{
  return k.are_parallel_2_object()(r1, r2);
}

template <typename K>
inline
bool
parallel(const typename CGAL_WRAP(K)::Segment_2 &s1,
         const typename CGAL_WRAP(K)::Segment_2 &s2, const K &k)
{
  return k.are_parallel_2_object()(s1, s2);
}

CGAL_END_NAMESPACE

#endif  // CGAL_KERNEL_GLOBAL_FUNCTIONS_2_H
