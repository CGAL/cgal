// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-0.9-I-06 $
// release_date  : $CGAL_Date: 1998/03/11 $
//
// file          : include/CGAL/Polygon_2_algorithms.C
// source        :
// revision      : 1.8a
// revision_date : 13 Mar 1998
// author(s)     : Wieger Wesselink <wieger@cs.ruu.nl>
//
// coordinator   : Utrecht University
//
// ======================================================================

#include "CGAL/Polygon_2_algorithms.h"
#include "CGAL/Polygon_2_simplicity.h"
#include <cstdlib>
#include <algorithm>
#include <set>
#include <vector>

CGAL_BEGIN_NAMESPACE


//-----------------------------------------------------------------------//
//                          is_simple_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy_2
//      Traits::less_xy_2
//      Traits::Orientation_2
//      Traits::orientation_2
//      Traits::Point_2

template <class ForwardIterator, class Traits>
bool is_simple_2(ForwardIterator first,
                      ForwardIterator last,
                      const Traits& traits)
{
  return is_simple_polygon(first, last, traits);
}

//-----------------------------------------------------------------------//
//                          left_vertex_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy

template <class ForwardIterator, class Traits>
ForwardIterator left_vertex_2(ForwardIterator first,
                                   ForwardIterator last,
                                   const Traits&)
{
  CGAL_polygon_precondition(first != last);

  typedef typename Traits::Less_xy Less_xy;
  return std::min_element(first, last, Less_xy());
}

//-----------------------------------------------------------------------//
//                          right_vertex_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy

template <class ForwardIterator, class Traits>
ForwardIterator right_vertex_2(ForwardIterator first,
                                    ForwardIterator last,
                                    const Traits&)
{
  CGAL_polygon_precondition(first != last);

  typedef typename Traits::Less_xy Less_xy;
  return std::max_element(first, last, Less_xy());
}

//-----------------------------------------------------------------------//
//                          top_vertex_2
//-----------------------------------------------------------------------//
// uses Traits::Less_yx

template <class ForwardIterator, class Traits>
ForwardIterator top_vertex_2(ForwardIterator first,
                                  ForwardIterator last,
                                  const Traits&)
{
  CGAL_polygon_precondition(first != last);

  typedef typename Traits::Less_yx Less_yx;
  return std::max_element(first, last, Less_yx());
}

//-----------------------------------------------------------------------//
//                          bottom_vertex_2
//-----------------------------------------------------------------------//
// uses Traits::Less_yx

template <class ForwardIterator, class Traits>
ForwardIterator bottom_vertex_2(ForwardIterator first,
                                     ForwardIterator last,
                                     const Traits&)
{
  CGAL_polygon_precondition(first != last);

  typedef typename Traits::Less_yx Less_yx;
  return std::min_element(first, last, Less_yx());
}

//-----------------------------------------------------------------------//
//                          bbox_2
//-----------------------------------------------------------------------//

template <class InputIterator>
Bbox_2 bbox_2(InputIterator first, InputIterator last)
{
  CGAL_polygon_precondition(first != last);
  Bbox_2 result = (*first).bbox();

  while (++first != last)
    result = result + (*first).bbox();

  return result;
}

//-----------------------------------------------------------------------//
//                          area_2
//-----------------------------------------------------------------------//
// uses Traits::determinant_2

// template <class ForwardIterator, class FT, class Traits>
// void area_2(ForwardIterator first,
//                 ForwardIterator last,
//                 FT& result,
//                 const Traits& traits)
//{

//-----------------------------------------------------------------------//
//                          is_convex_2
//-----------------------------------------------------------------------//
// uses Traits::lexicographically_xy_smaller
//      Traits::orientation

template <class ForwardIterator, class Traits>
bool is_convex_2(ForwardIterator first,
                      ForwardIterator last,
                      const Traits& traits)
{
  ForwardIterator previous = first;
  if (previous == last) return true;

  ForwardIterator current = previous; ++current;
  if (current == last) return true;

  ForwardIterator next = current; ++next;
  if (next == last) return true;

  // initialization
  bool HasClockwiseTriples = false;
  bool HasCounterClockwiseTriples = false;
  bool Order = traits.lexicographically_xy_smaller(*previous, *current);
  int NumOrderChanges = 0;

  do {
    switch (traits.orientation(*previous, *current, *next)) {
      case CLOCKWISE:
        HasClockwiseTriples = true;
        break;
      case COUNTERCLOCKWISE:
        HasCounterClockwiseTriples = true;
        break;
      default:
	;
    }

    bool NewOrder = traits.lexicographically_xy_smaller(*current, *next);
    if (Order != NewOrder) NumOrderChanges++;

    if (NumOrderChanges > 2) {
#ifdef CGAL_POLYGON_DEBUG
cout << "too many order changes: not convex!" << endl;
#endif
      return false;
    }

    if (HasClockwiseTriples && HasCounterClockwiseTriples) {
#ifdef CGAL_POLYGON_DEBUG
cout << "polygon not locally convex!" << endl;
#endif
      return false;
    }

    previous = current;
    current = next;
    ++next;
    if (next == last) next = first;
    Order = NewOrder;
  }
  while (previous != first);

  return true;
}

//-----------------------------------------------------------------------//
//                          oriented_side_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy
//      Traits::compare_x
//      Traits::compare_y
//      Traits::determinant_2
//      Traits::orientation
//      Traits::sign

template <class ForwardIterator, class Point, class Traits>
Oriented_side oriented_side_2(ForwardIterator first,
                                        ForwardIterator last,
                                        const Point& point,
                                        const Traits& traits)
{
  Oriented_side result;

  Orientation o = orientation_2(first, last, traits);
  CGAL_polygon_assertion(o != COLLINEAR);

  Bounded_side b = bounded_side_2(first, last, point, traits);
  switch (b) {
    case ON_BOUNDARY:
      result = ON_ORIENTED_BOUNDARY;
      break;

    case ON_BOUNDED_SIDE:
      result = (o == CLOCKWISE) ?
               ON_NEGATIVE_SIDE :
               ON_POSITIVE_SIDE;
      break;

    case ON_UNBOUNDED_SIDE:
      result = (o == CLOCKWISE) ?
               ON_POSITIVE_SIDE :
               ON_NEGATIVE_SIDE;
      break;
  }

  return result;
}

//-----------------------------------------------------------------------//
//                          bounded_side_2
//-----------------------------------------------------------------------//
// uses Traits::compare_x
//      Traits::compare_y
//      Traits::determinant_2
//      Traits::sign
//
// returns ON_BOUNDED_SIDE, ON_BOUNDARY or ON_UNBOUNDED_SIDE

template <class ForwardIterator, class Point, class Traits>
Bounded_side bounded_side_2(ForwardIterator first,
                                      ForwardIterator last,
                                      const Point& point,
                                      const Traits& traits)
{
  ForwardIterator current = first;
  if (current == last) return ON_UNBOUNDED_SIDE;

  ForwardIterator next = current; ++next;
  if (next == last) return ON_UNBOUNDED_SIDE;

  bool IsInside = false;
  Comparison_result CompareCurrent = traits.compare_y(*current, point);

  do // check if the segment (current,next) intersects
     // the ray { (t,y) | t >= point.x() }
  {
    Comparison_result CompareNext = traits.compare_y(*next, point);

    switch (CompareCurrent) {
      case SMALLER:
        switch (CompareNext) {
          case SMALLER:
            break;
          case EQUAL:
            switch (traits.compare_x(point, *next)) {
              case SMALLER: IsInside = !IsInside; break;
              case EQUAL:   return ON_BOUNDARY;
              case LARGER:  break;
            }
            break;
          case LARGER:
            if (point.x() < std::min((*current).x(), (*next).x())) {
              IsInside = !IsInside;
            }
            else if (point.x() <= std::max((*current).x(),(*next).x())) {
              switch (traits.sign(traits.determinant_2(point,
                                                       *current,
                                                       *next)))
              {
                case 0: return ON_BOUNDARY;
                case 1: IsInside = !IsInside; break;
              }
            }
            break;
        }
        break;
      case EQUAL:
        switch (CompareNext) {
          case SMALLER:
            switch (traits.compare_x(point, *current)) {
              case SMALLER: IsInside = !IsInside; break;
              case EQUAL:   return ON_BOUNDARY;
              case LARGER:  break;
            }
            break;
          case EQUAL:
            if ( (std::min((*current).x(), (*next).x()) <= point.x()) &&
                 (point.x() <= std::max((*current).x(), (*next).x()))    ) {
              return ON_BOUNDARY;
            }
            break;
          case LARGER:
            if (point.x() == (*current).x()) {
              return ON_BOUNDARY;
            }
            break;
        }
        break;
      case LARGER:
        switch (CompareNext) {
          case SMALLER:
            if (point.x() < std::min((*current).x(), (*next).x())) {
              IsInside = !IsInside;
            }
            else if (point.x() <= std::max((*current).x(),(*next).x())) {
              switch (traits.sign(traits.determinant_2(point,
                                                       *current,
                                                       *next)))
              {
                case -1: IsInside = !IsInside; break;
                case  0: return ON_BOUNDARY;
              }
            }
            break;
          case EQUAL:
            if (point.x() == (*next).x()) {
              return ON_BOUNDARY;
            }
            break;
          case LARGER:
            break;
        }
        break;
    }

    current = next;
    CompareCurrent = CompareNext;
    ++next;
    if (next == last) next = first;   
  }
  while (current != first);

  return IsInside ? ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
}

//-----------------------------------------------------------------------//
//                          orientation_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy
//      Traits::orientation

template <class ForwardIterator, class Traits>
Orientation orientation_2(ForwardIterator first,
                                    ForwardIterator last,
                                    const Traits& traits)
{
  CGAL_polygon_precondition(is_simple_2(first, last, traits));

  ForwardIterator i = left_vertex_2(first, last, traits);

  ForwardIterator prev = (i == first) ? last : i;
  --prev;

  ForwardIterator next = i;
  ++next;
  if (next == last)
    next = first;

  // if the range [first,last) contains less than three points, then some
  // of the points (prev,i,next) will coincide

  // return the orientation of the triple (prev,i,next)
  return traits.orientation(*prev, *i, *next);
}

CGAL_END_NAMESPACE

