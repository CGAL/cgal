// Copyright (c) 2003
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_KERNEL_GLOBAL_FUNCTIONS_H
#define CGAL_KERNEL_GLOBAL_FUNCTIONS_H

// Generic functions calling the kernel functor.
// Independent of the dimension.

#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Kernel/hash_functions.h>

namespace CGAL {

template <class T1, class T2, class T3>
inline
Comparison_result
compare_distance(const T1 &o1,
                 const T2 &o2,
                 const T3 &o3)
{
  typedef typename Kernel_traits<T1>::Kernel K;
  return internal::compare_distance(o1, o2, o3, K());
}

template <class T1, class T2, class T3, class T4>
inline
Comparison_result
compare_distance(const T1 &o1,
                 const T2 &o2,
                 const T3 &o3,
                 const T4 &o4)
{
  typedef typename Kernel_traits<T1>::Kernel K;
  return internal::compare_distance(o1, o2, o3, o4, K());
}

template <typename O>
inline
bool
parallel(const O &o1, const O &o2)
{
  typedef typename Kernel_traits<O>::Kernel K;
  return internal::parallel(o1, o2, K());
}

} //namespace CGAL

#endif  // CGAL_KERNEL_GLOBAL_FUNCTIONS_H
