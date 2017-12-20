// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_KERNEL_INTERSECTIONS_H
#define CGAL_CIRCULAR_KERNEL_INTERSECTIONS_H 

#include <CGAL/license/Circular_kernel_2.h>


#include <vector>

#include <CGAL/Circle_2.h>
#include <CGAL/Circular_arc_2.h>
#include <CGAL/Line_arc_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Circular_arc_point_2.h>
#include <CGAL/Circular_kernel_2/Intersection_traits.h>

namespace CGAL {

#define CGAL_CIRCULAR_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_(A,B) \
template < class OutputIterator, class K > \
OutputIterator \
intersection(const A <K> &c1, const B <K> &c2, OutputIterator res) \
{ \
  return typename K::Intersect_2()(c1, c2, res); \
} \
namespace internal { \
  template <class K> \
  inline \
  bool \
  do_intersect(const typename K::A &c1, const typename K::B &c2, const K&) \
  { \
    std::vector< typename CK2_Intersection_traits<K, typename K::A, typename K::B>::type > res; \
    typename K::Intersect_2()(c1,c2,std::back_inserter(res)); \
		return !res.empty(); \
  } \
} \
template <class K> \
inline \
bool \
do_intersect(const A <K> &c1, const B <K> &c2) \
{ \
  return typename K::Do_intersect_2()(c1, c2); \
}

// Circle_2 Circle_2 already has its do_intersect
// so it needs only the global intersection
template < class OutputIterator, class K >
OutputIterator
intersection(const Circle_2 <K> &c1, const Circle_2 <K> &c2, OutputIterator res)
{
  return typename K::Intersect_2()(c1, c2, res);
}

template < class OutputIterator, class K >
OutputIterator
intersection(const Circle_2 <K> &c1, const Line_2 <K> &c2, OutputIterator res)
{
  return typename K::Intersect_2()(c1, c2, res);
}

template < class OutputIterator, class K >
OutputIterator
intersection(const Line_2 <K> &c1, const Circle_2 <K> &c2, OutputIterator res)
{
  return typename K::Intersect_2()(c1, c2, res);
}

template < class OutputIterator, class K >
OutputIterator
intersection(const Line_2 <K> &c1, const Line_2 <K> &c2, OutputIterator res)
{
  return typename K::Intersect_2()(c1, c2, res);
}

CGAL_CIRCULAR_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_(Circular_arc_2, Circular_arc_2)
CGAL_CIRCULAR_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_(Line_arc_2, Line_arc_2)
CGAL_CIRCULAR_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_(Line_arc_2, Circle_2)
CGAL_CIRCULAR_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_(Circle_2, Line_arc_2)
CGAL_CIRCULAR_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_(Circular_arc_2, Circle_2)
CGAL_CIRCULAR_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_(Circle_2, Circular_arc_2)
CGAL_CIRCULAR_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_(Line_arc_2, Circular_arc_2)
CGAL_CIRCULAR_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_(Circular_arc_2, Line_arc_2)
CGAL_CIRCULAR_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_(Line_2, Circular_arc_2)
CGAL_CIRCULAR_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_(Line_2, Line_arc_2)
CGAL_CIRCULAR_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_(Circular_arc_2, Line_2)
CGAL_CIRCULAR_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_(Line_arc_2, Line_2)


} //namespace CGAL

#undef CGAL_CIRCULAR_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_
#endif // CGAL_CIRCULAR_KERNEL_INTERSECTIONS_H
