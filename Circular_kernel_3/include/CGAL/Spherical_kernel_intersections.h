// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a
// STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_INTERNAL_FUNCTIONS_ON_INTERSECTION_3_H
#define CGAL_SPHERICAL_KERNEL_INTERNAL_FUNCTIONS_ON_INTERSECTION_3_H

#include <CGAL/license/Circular_kernel_3.h>


#include <CGAL/Sphere_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Circle_3.h>
#include <CGAL/Circular_arc_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Line_arc_3.h>
#include <CGAL/Circular_arc_point_3.h>

namespace CGAL {

#define CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_JUST_INTERSECTION_2_(A,B) \
template < class OutputIterator, class K > \
OutputIterator \
intersection(const A <K> &c1, const B <K> &c2, OutputIterator res) \
{ \
  return typename K::Intersect_3()(c1, c2, res); \
}

#define CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(A,B) \
template < class OutputIterator, class K > \
OutputIterator \
intersection(const A <K> &c1, const B <K> &c2, OutputIterator res) \
{ \
  return typename K::Intersect_3()(c1, c2, res); \
} \
template <class K> \
inline \
bool \
do_intersect(const A <K> &c1, const B <K> &c2) \
{ \
  return typename K::Do_intersect_3()(c1, c2); \
}

#define CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_3_(A,B,C) \
template < class OutputIterator, class K > \
OutputIterator \
intersection(const A <K> &c1, const B <K> &c2, const C <K> &c3, OutputIterator res) \
{ \
  return typename K::Intersect_3()(c1, c2, c3, res); \
} \
template <class K> \
inline \
bool \
do_intersect(const A <K> &c1, const B <K> &c2, const C <K> &c3) \
{ \
  return typename K::Do_intersect_3()(c1, c2, c3); \
}

CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_JUST_INTERSECTION_2_(Sphere_3, Line_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_JUST_INTERSECTION_2_(Line_3, Sphere_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_3_(Sphere_3, Sphere_3, Sphere_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_3_(Sphere_3, Sphere_3, Plane_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_3_(Plane_3, Sphere_3, Sphere_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_3_(Plane_3, Plane_3, Sphere_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_3_(Sphere_3, Plane_3, Plane_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Circle_3, Plane_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Plane_3, Circle_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Circle_3, Sphere_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Sphere_3, Circle_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Circle_3, Circle_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Circle_3, Line_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Line_3, Circle_3)

CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Line_arc_3, Line_arc_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Line_3, Line_arc_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Line_arc_3, Line_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Circle_3, Line_arc_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Line_arc_3, Circle_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Sphere_3, Line_arc_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Line_arc_3, Sphere_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Plane_3, Line_arc_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Line_arc_3, Plane_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Circular_arc_3, Circular_arc_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Line_3, Circular_arc_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Circular_arc_3, Line_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Circle_3, Circular_arc_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Circular_arc_3, Circle_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Sphere_3, Circular_arc_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Circular_arc_3, Sphere_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Plane_3, Circular_arc_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Circular_arc_3, Plane_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Circular_arc_3, Line_arc_3)
CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_(Line_arc_3, Circular_arc_3)

} //namespace CGAL

#undef CGAL_SPHERICAL_KERNEL_MACRO_GLOBAL_FUNCTION_INTERSECTION_2_
#endif // CGAL_SPHERICAL_KERNEL_INTERNAL_FUNCTIONS_ON_INTERSECTION_3_H
