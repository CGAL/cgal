// Copyright (c) 2025
// GeometryFactory (France),
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Léo Valque

#ifndef CGAL_KERNEL_GLOBAL_FUNCTIONS_DISTANCE_3_H
#define CGAL_KERNEL_GLOBAL_FUNCTIONS_DISTANCE_3_H

// Distance functions calling the kernel functor.

#define CGAL_SQUARED_DISTANCE_FUNCTION(A, B)                            \
template <class K>                                                      \
inline                                                                  \
typename K::FT                                                          \
squared_distance(const A<K>& a, const B<K>& b)                          \
{                                                                       \
  return K().compute_squared_distance_3_object()(a, b);                 \
}                                                                       \
template <class K>                                                      \
inline                                                                  \
typename K::FT                                                          \
squared_distance(const B<K>& a, const A<K>& b)                          \
{                                                                       \
  return K().compute_squared_distance_3_object()(b, a);                 \
}

#define CGAL_SQUARED_DISTANCE_FUNCTION_SELF(A)                          \
template <class K>                                                      \
inline                                                                  \
typename K::FT                                                          \
squared_distance(const A<K>& a, const A<K>& b)                          \
{                                                                       \
  return K().compute_squared_distance_3_object()(a, b);                 \
}

#define CGAL_COMPARE_SQUARED_DISTANCE_FUNCTION(A, B)                    \
template <class K>                                                      \
inline                                                                  \
typename K::Comparison_result                                           \
compare_squared_distance(const A<K>& a,                                 \
                         const B<K>& b,                                 \
                         const typename K::FT& d2)                      \
{                                                                       \
  return K().compare_squared_distance_3_object()(a, b, d2);             \
}                                                                       \
template <class K>                                                      \
inline                                                                  \
typename K::Comparison_result                                           \
compare_squared_distance(const B<K>& b,                                 \
                         const A<K>& a,                                 \
                         const typename K::FT& d2)                      \
{                                                                       \
  return K().compare_squared_distance_3_object()(b, a, d2);             \
}

#define CGAL_COMPARE_SQUARED_DISTANCE_FUNCTION_SELF(A)                  \
template <class K>                                                      \
inline                                                                  \
typename K::Comparison_result                                           \
compare_squared_distance(const A<K>& a,                                 \
                         const A<K>& b,                                 \
                         const typename K::FT& d2)                      \
{                                                                       \
  return K().compare_squared_distance_3_object()(a, b, d2);             \
}

#define CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION_SELF(A)      \
CGAL_SQUARED_DISTANCE_FUNCTION_SELF(A)                                  \
CGAL_COMPARE_SQUARED_DISTANCE_FUNCTION_SELF(A)

#define CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION(A, B)        \
CGAL_SQUARED_DISTANCE_FUNCTION(A, B)                                    \
CGAL_COMPARE_SQUARED_DISTANCE_FUNCTION(A, B)

namespace CGAL {

CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION_SELF(Point_3)
CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION_SELF(Plane_3)
CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION_SELF(Segment_3)
CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION_SELF(Line_3)
CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION_SELF(Triangle_3)

CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION(Point_3, Line_3)
CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION(Point_3, Ray_3)
CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION(Point_3, Segment_3)
CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION(Point_3, Plane_3)
CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION(Point_3, Triangle_3)
CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION(Point_3, Tetrahedron_3)

CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION(Segment_3, Ray_3)

CGAL_COMPUTE_AND_COMPARE_SQUARED_DISTANCE_FUNCTION(Line_3, Plane_3)

} //namespace CGAL

#endif  // CGAL_KERNEL_GLOBAL_FUNCTIONS_DISTANCE_3_H