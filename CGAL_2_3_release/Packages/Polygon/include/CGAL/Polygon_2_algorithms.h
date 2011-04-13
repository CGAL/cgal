// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-0.9-I-06 $
// release_date  : $CGAL_Date: 1998/03/11 $
//
// file          : include/CGAL/Polygon_2_algorithms.h
// source        : 
// revision      : 1.8a
// revision_date : 13 Mar 1998
// author(s)     : Wieger Wesselink <wieger@cs.ruu.nl>
//
// coordinator   : Utrecht University
//
// ============================================================================

#ifndef CGAL_POLYGON_2_ALGORITHMS_H
#define CGAL_POLYGON_2_ALGORITHMS_H

#include <CGAL/basic.h>

#include <CGAL/enum.h>
#include <CGAL/Bbox_2.h>

#include <CGAL/polygon_assertions.h>

#ifdef CGAL_REP_CLASS_DEFINED
 #include <CGAL/Polygon_traits_2.h>
#endif // CGAL_REP_CLASS_DEFINED

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------//
//                  algorithms for sequences of 2D points
//-----------------------------------------------------------------------//

template <class ForwardIterator, class Traits>
ForwardIterator left_vertex_2(ForwardIterator first,
                                   ForwardIterator last,
                                   const Traits& traits);

template <class ForwardIterator, class Traits>
ForwardIterator right_vertex_2(ForwardIterator first,
                                    ForwardIterator last,
                                    const Traits& traits);

template <class ForwardIterator, class Traits>
ForwardIterator top_vertex_2(ForwardIterator first,
                                  ForwardIterator last,
                                  const Traits& traits);

template <class ForwardIterator, class Traits>
ForwardIterator bottom_vertex_2(ForwardIterator first,
                                     ForwardIterator last,
                                     const Traits& traits);

template <class InputIterator>
Bbox_2 bbox_2(InputIterator first, InputIterator last);


template <class ForwardIterator, class Traits>
void 
area_2( ForwardIterator first, ForwardIterator last,
   	typename Traits::FT &result,
        const Traits& traits)
{
   typedef typename Traits::FT FT;
   result = FT(0);
   // check if the polygon is empty
   if (first == last) return;
   ForwardIterator second = first; ++second;
   // check if the polygon has only one point
   if (second == last) return;
   typename Traits::Compute_area_2 compute_area_2 =
            traits.compute_area_2_object();
   typename Traits::Construct_triangle_2 construct_triangle_2 =
            traits.construct_triangle_2_object();
   ForwardIterator third = second;
   while (++third != last) {
	result = result + compute_area_2(
                    construct_triangle_2(*first, *second, *third));
	second = third;
   }
   result = result / FT(2);
}

template <class ForwardIterator, class Traits>
typename Traits::FT 
polygon_area_2( ForwardIterator first, ForwardIterator last,
        const Traits& traits)
{
   typedef typename Traits::FT FT;
   FT result = FT(0);
   // check if the polygon is empty
   if (first == last) return result;
   ForwardIterator second = first; ++second;
   // check if the polygon has only one point
   if (second == last) return result;
   typename Traits::Compute_area_2 compute_area_2 =
            traits.compute_area_2_object();
   typename Traits::Construct_triangle_2 construct_triangle_2 =
            traits.construct_triangle_2_object();
   ForwardIterator third = second;
   while (++third != last) {
	result = result + compute_area_2(
                    construct_triangle_2(*first, *second, *third));
	second = third;
   }
   result = result / FT(2);
   return result;
}

template <class ForwardIterator, class Traits>
bool is_convex_2(ForwardIterator first,
                      ForwardIterator last,
                      const Traits& traits);

template <class ForwardIterator, class Traits>
bool is_simple_2(ForwardIterator first,
                      ForwardIterator last,
                      const Traits& traits);

// In the following two functions we would like to use Traits::Point_2 instead
// of Point, but this is not allowed by g++ 2.7.2.

template <class ForwardIterator, class Point, class Traits>
Oriented_side oriented_side_2(ForwardIterator first,
                                        ForwardIterator last,
                                        const Point& point,
                                        const Traits& traits);

template <class ForwardIterator, class Point, class Traits>
Bounded_side bounded_side_2(ForwardIterator first,
                                      ForwardIterator last,
                                      const Point& point,
                                      const Traits& traits);

template <class ForwardIterator, class Traits>
Orientation orientation_2(ForwardIterator first,
                                    ForwardIterator last,
                                    const Traits& traits);

//-----------------------------------------------------------------------//
//                         implementation
//-----------------------------------------------------------------------//

#ifdef CGAL_REP_CLASS_DEFINED
template <class ForwardIterator, class R>
inline
ForwardIterator left_vertex_2_aux(ForwardIterator first,
                                       ForwardIterator last,
                                       const Point_2<R>)
{
  return left_vertex_2(first, last, R());
}

template <class ForwardIterator>
inline
ForwardIterator left_vertex_2(ForwardIterator first,
                                   ForwardIterator last)
{
  return left_vertex_2_aux(first, last, *first);
}

template <class ForwardIterator, class R>
inline
ForwardIterator right_vertex_2_aux(ForwardIterator first,
                                        ForwardIterator last,
                                        const Point_2<R>&)
{
  return right_vertex_2(first, last, Polygon_traits_2<R>());
}

template <class ForwardIterator>
inline
ForwardIterator right_vertex_2(ForwardIterator first,
                                    ForwardIterator last)
{
  return right_vertex_2_aux(first, last, *first);
}

template <class ForwardIterator, class R>
inline
ForwardIterator top_vertex_2_aux(ForwardIterator first,
                                      ForwardIterator last,
                                      const Point_2<R>&)
{
  return top_vertex_2(first, last, Polygon_traits_2<R>());
}

template <class ForwardIterator>
inline
ForwardIterator top_vertex_2(ForwardIterator first,
                                  ForwardIterator last)
{
  return top_vertex_2_aux(first, last, *first);
}

template <class ForwardIterator, class R>
inline
ForwardIterator bottom_vertex_2_aux(ForwardIterator first,
                                         ForwardIterator last,
                                         const Point_2<R>&)
{
  return bottom_vertex_2(first, last, Polygon_traits_2<R>());
}

template <class ForwardIterator>
inline
ForwardIterator bottom_vertex_2(ForwardIterator first,
                                     ForwardIterator last)
{
  return bottom_vertex_2_aux(first, last, *first);
}

template <class ForwardIterator, class Numbertype, class R>
inline
void area_2_aux(ForwardIterator first,
                     ForwardIterator last,
                     Numbertype& result,
                     const Point_2<R>&)
{
  area_2(first, last, result, Polygon_traits_2<R>());
}

template <class ForwardIterator, class Numbertype>
inline
void area_2(ForwardIterator first,
                 ForwardIterator last,
                 Numbertype& result)
{
  area_2_aux(first, last, result, *first);
}

template <class ForwardIterator, class R>
inline
bool is_convex_2_aux(ForwardIterator first,
                          ForwardIterator last,
                          const Point_2<R>&)
{
  return is_convex_2(first, last, Polygon_traits_2<R>());
}

template <class ForwardIterator>
inline
bool is_convex_2(ForwardIterator first,
                      ForwardIterator last)
{
  return is_convex_2_aux(first, last, *first);
}

template <class ForwardIterator, class R>
inline
bool is_simple_2_aux(ForwardIterator first,
                          ForwardIterator last,
                          const Point_2<R>&)
{
  return is_simple_2(first, last, R());
}

template <class ForwardIterator>
inline
bool is_simple_2(ForwardIterator first,
                      ForwardIterator last)
{
  return is_simple_2_aux(first, last, *first);
}

template <class ForwardIterator, class R>
inline
Oriented_side oriented_side_2(ForwardIterator first,
                                        ForwardIterator last,
                                        const Point_2<R>& point)
{
  return oriented_side_2(first, last, point, Polygon_traits_2<R>());
}

template <class ForwardIterator, class R>
inline
Bounded_side bounded_side_2(ForwardIterator first,
                                      ForwardIterator last,
                                      const Point_2<R>& point)
{
  return bounded_side_2(first, last, point, Polygon_traits_2<R>());
}

template <class ForwardIterator, class R>
inline
Orientation orientation_2_aux(ForwardIterator first,
                                        ForwardIterator last,
                                        const Point_2<R>&)
{
  return orientation_2(first, last, Polygon_traits_2<R>());
}

template <class ForwardIterator>
inline
Orientation orientation_2(ForwardIterator first,
                                    ForwardIterator last)
{
  return orientation_2_aux(first, last, *first);
}
#endif // CGAL_REP_CLASS_DEFINED

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION 
#include <CGAL/Polygon_2_algorithms.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif // CGAL_POLYGON_2_ALGORITHMS_H


