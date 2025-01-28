// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Baruch Zukerman <baruchzu@post.tau.ac.il>
//            Ron Wein        <wein@post.tau.ac.il>
//            Efi Fogel       <efif@post.tau.ac.il>
//            Simon Giraudot  <simon.giraudot@geometryfactory.com>

#ifndef CGAL_BSO_INTERNAL_FUNCTIONS_H
#define CGAL_BSO_INTERNAL_FUNCTIONS_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <iterator>

#include <CGAL/Boolean_set_operations_2/Gps_default_traits.h>
#include <CGAL/Boolean_set_operations_2/Polygon_conversions.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/Arr_polyline_traits_2.h>

namespace CGAL {

/// \name s_do_intersect() and r_do_intersect() functions.
//@{

// Single
// With Traits
template <typename Pgn1, class Pgn2, typename Traits>
inline bool s_do_intersect(const Pgn1& pgn1, const Pgn2& pgn2, Traits& traits) {
  General_polygon_set_2<Traits> gps(pgn1, traits);
  return gps.do_intersect(pgn2);
}

// Without Traits
template <typename Pgn1, typename Pgn2>
inline bool s_do_intersect(const Pgn1& pgn1, const Pgn2& pgn2) {
  typename Gps_polyline_traits<Pgn1>::Traits traits;
  const typename Gps_polyline_traits<Pgn1>::Polyline_traits& ptraits(traits);
  return s_do_intersect(convert_polygon(pgn1, ptraits),
                        convert_polygon(pgn2, ptraits), traits);
}

// Single Range
// With Traits
template <typename InputIterator, typename Traits>
inline bool r_do_intersect(InputIterator begin, InputIterator end,
                           Traits& traits, unsigned int k=5) {
  if (begin == end) return false;
  General_polygon_set_2<Traits> gps(*begin, traits);
  return gps.do_intersect(std::next(begin), end, k);
}

// Without Traits
template <typename InputIterator>
inline bool r_do_intersect(InputIterator begin, InputIterator end,
                           unsigned int k=5) {
  typedef typename std::iterator_traits<InputIterator>::value_type Pgn;
  typename Gps_polyline_traits<Pgn>::Traits traits;
  const typename Gps_polyline_traits<Pgn>::Polyline_traits& ptraits(traits);
  return r_do_intersect(convert_polygon_iterator(begin, ptraits),
                        convert_polygon_iterator(end, ptraits), traits, k);
}

// Souble Range
// With Traits
template <typename InputIterator1, typename InputIterator2, typename Traits>
inline bool r_do_intersect(InputIterator1 begin1, InputIterator1 end1,
                           InputIterator2 begin2, InputIterator2 end2,
                           Traits& traits, unsigned int k=5) {
  if (begin1 == end1) return do_intersect(begin2, end2, traits, k);
  General_polygon_set_2<Traits> gps(*begin1, traits);
  return gps.do_intersect(std::next(begin1), end1, begin2, end2, k);
}

// Without Traits
template <typename InputIterator1, typename InputIterator2>
inline bool r_do_intersect (InputIterator1 begin1, InputIterator1 end1,
                            InputIterator2 begin2, InputIterator2 end2,
                            unsigned int k=5) {
  typedef typename std::iterator_traits<InputIterator1>::value_type Pgn;
  typename Gps_polyline_traits<Pgn>::Traits traits;
  const typename Gps_polyline_traits<Pgn>::Polyline_traits& ptraits(traits);
  return r_do_intersect(convert_polygon_iterator(begin1, ptraits),
                        convert_polygon_iterator(end1, ptraits),
                        convert_polygon_iterator(begin2, ptraits),
                        convert_polygon_iterator(end2, ptraits), traits, k);
}

//@}

/// \name _oriented_side() functions.
//@{

// With Traits
template <typename Obj, typename Pgn, typename Traits>
inline Oriented_side _oriented_side(const Obj& obj, const Pgn& pgn,
                                    Traits& traits) {
  General_polygon_set_2<Traits> gps(pgn, traits);
  return gps.oriented_side(obj);
}

// Without Traits (point, polygon)
template <typename Kernel, typename Pgn>
inline Oriented_side _oriented_side(const Point_2<Kernel>& point,
                                    const Pgn& pgn) {
  // Use the polygon to determine the (default) traits
  typename Gps_polyline_traits<Pgn>::Traits traits;
  const typename Gps_polyline_traits<Pgn>::Polyline_traits& ptraits(traits);
  return _oriented_side(point, convert_polygon(pgn, ptraits), traits);
}

// Without Traits (polygon, polygon)
template <typename Pgn1, typename Pgn2>
inline Oriented_side _oriented_side(const Pgn1& pgn1, const Pgn2& pgn2)
{
  // Use the first polygon to determine the (default) traits
  typename Gps_polyline_traits<Pgn1>::Traits traits;
  const typename Gps_polyline_traits<Pgn1>::Polyline_traits& ptraits(traits);
  return _oriented_side(convert_polygon(pgn1, ptraits),
                        convert_polygon(pgn2, ptraits), traits);
}

//@}

/// \name s_intersection() functions.
//@{

// Single
// With Traits
template <typename Pgn1, typename Pgn2, typename OutputIterator, typename Traits>
inline OutputIterator s_intersection(const Pgn1& pgn1, const Pgn2& pgn2,
                                     OutputIterator oi, Traits& traits) {
  General_polygon_set_2<Traits> gps(pgn1, traits);
  gps.intersection(pgn2);
  return gps.polygons_with_holes(oi);
}

// Without Traits
template <typename Kernel, typename Container,
          typename Pgn1, typename Pgn2, typename OutputIterator>
inline OutputIterator s_intersection(const Pgn1& pgn1, const Pgn2& pgn2,
                                     OutputIterator oi) {
  // Use the first polygon to determine the (default) traits
  typedef typename Gps_polyline_traits<Pgn1>::Polyline_traits   Polyline_traits;

  typename Gps_polyline_traits<Pgn1>::Traits traits;
  const Polyline_traits& ptraits(traits);
  s_intersection(convert_polygon(pgn1, ptraits), convert_polygon(pgn2, ptraits),
                 convert_polygon_back(oi, pgn1), traits);
  return oi;
}

// Single Range
// With Traits
template <typename InputIterator, typename OutputIterator, typename Traits>
inline OutputIterator r_intersection(InputIterator begin, InputIterator end,
                                     OutputIterator oi, Traits& traits,
                                     unsigned int k=5) {
  if (begin == end) return (oi);
  General_polygon_set_2<Traits> gps(*begin, traits);
  gps.intersection(std::next(begin), end, k);
  return gps.polygons_with_holes(oi);
}

// Without Traits
template <typename InputIterator, typename OutputIterator>
inline OutputIterator r_intersection(InputIterator begin, InputIterator end,
                                     OutputIterator oi, unsigned int k=5) {
  typedef typename std::iterator_traits<InputIterator>::value_type Pgn;
  typename Gps_polyline_traits<Pgn>::Traits traits;
  const typename Gps_polyline_traits<Pgn>::Polyline_traits& ptraits(traits);
  if (begin == end) return (oi);
  return r_intersection(convert_polygon_iterator(begin, ptraits),
                        convert_polygon_iterator(end, ptraits),
                        convert_polygon_back(oi, *begin), traits, k);
}

// Souble Range
// With Traits
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename Traits>
inline OutputIterator r_intersection(InputIterator1 begin1, InputIterator1 end1,
                                     InputIterator2 begin2, InputIterator2 end2,
                                     OutputIterator oi, Traits& traits,
                                     unsigned int k=5) {
  if (begin1 == end1) return r_intersection(begin2, end2, oi, traits, k);
  General_polygon_set_2<Traits> gps(*begin1, traits);
  gps.intersection(std::next(begin1), end1, begin2, end2, k);
  return gps.polygons_with_holes(oi);
}

// Without Traits
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline OutputIterator
r_intersection(InputIterator1 begin1, InputIterator1 end1,
               InputIterator2 begin2, InputIterator2 end2,
               OutputIterator oi, unsigned int k=5) {
  typedef typename std::iterator_traits<InputIterator1>::value_type Pgn;
  typename Gps_polyline_traits<Pgn>::Traits traits;
  const typename Gps_polyline_traits<Pgn>::Polyline_traits& ptraits(traits);
  if (begin1 == end1) {
      if (begin2 == end2) return oi;
      return r_intersection(convert_polygon_iterator(begin2, ptraits),
                            convert_polygon_iterator(end2, ptraits),
                            convert_polygon_back(oi, *begin2), traits, k);
  }
  return r_intersection(convert_polygon_iterator(begin1, ptraits),
                        convert_polygon_iterator(end1, ptraits),
                        convert_polygon_iterator(begin2, ptraits),
                        convert_polygon_iterator(end2, ptraits),
                        convert_polygon_back(oi, *begin1), traits, k);
}

//@}

/// \name is_empty() functions.
//@{

// Polygon_2
template <typename Traits>
inline bool _is_empty(const typename Traits::Polygon_2& pgn, Traits& traits) {
  typedef typename Traits::Curve_const_iterator Curve_const_iterator;
  const std::pair<Curve_const_iterator, Curve_const_iterator>& itr_pair =
    traits.construct_curves_2_object()(pgn);
  return (itr_pair.first == itr_pair.second);
}

// A polygon with holes cannot be empty.
template <typename Traits>
inline bool _is_empty(const typename Traits::Polygon_with_holes_2&, Traits&)
{ return false; }

//@}

/// \name s_join() and r_join() functions.
//@{

// Single operands
// With Traits
template <typename Pgn1, typename Pgn2, typename Traits>
inline bool s_join(const Pgn1& pgn1, const Pgn2& pgn2,
                   typename Traits::Polygon_with_holes_2& res, Traits& traits) {
  if (_is_empty(pgn1, traits) || _is_empty(pgn2, traits)) return false;

  General_polygon_set_2<Traits> gps(pgn1, traits);
  gps.join(pgn2);
  if (gps.number_of_polygons_with_holes() == 1) {
    Oneset_iterator<typename Traits::Polygon_with_holes_2> oi(res);
    gps.polygons_with_holes(oi);
    return true;
  }

  // The polygon doesn't intersect; the original pgn1, pgn2 contain the union
  return false;
}

// Without Traits
template <typename Kernel, typename Container,
          typename Pgn1, typename Pgn2, typename Pwh>
inline bool s_join(const Pgn1& pgn1, const Pgn2& pgn2, Pwh& pwh) {
  // Use the first polygon to determine the (default) traits
  typedef typename Gps_polyline_traits<Pgn1>::Polyline_traits   Polyline_traits;
  typedef General_polygon_2<Polyline_traits>                    General_pgn;
  typedef General_polygon_with_holes_2<General_pgn>             General_pwh;

  General_pwh general_pwh;
  typename Gps_polyline_traits<Pgn1>::Traits traits;
  const Polyline_traits& ptraits(traits);
  auto res = s_join(convert_polygon(pgn1, ptraits),
                    convert_polygon(pgn2, ptraits),
                    general_pwh, traits);
  pwh = convert_polygon_back<Kernel, Container>(general_pwh);
  return res;
}

// Range
// With traits
template <typename InputIterator, typename OutputIterator, typename Traits>
inline OutputIterator r_join(InputIterator begin, InputIterator end,
                             OutputIterator oi, Traits& traits,
                             unsigned int k=5) {
  if (begin == end) return oi;
  General_polygon_set_2<Traits> gps(*begin, traits);
  gps.join(std::next(begin), end, k);
  return gps.polygons_with_holes(oi);
}

// Without traits
template <typename InputIterator, typename OutputIterator>
inline OutputIterator r_join(InputIterator begin, InputIterator end,
                             OutputIterator oi, unsigned int k=5) {
  typedef typename std::iterator_traits<InputIterator>::value_type Pgn;
  typename Gps_polyline_traits<Pgn>::Traits traits;
  const typename Gps_polyline_traits<Pgn>::Polyline_traits& ptraits(traits);

  if (begin == end) return oi;

  return r_join(convert_polygon_iterator(begin, ptraits),
                convert_polygon_iterator(end, ptraits),
                convert_polygon_back(oi, *begin), traits, k);
}

// Double Range
// With traits
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename Traits>
inline OutputIterator r_join(InputIterator1 begin1, InputIterator1 end1,
                             InputIterator2 begin2, InputIterator2 end2,
                             OutputIterator oi, Traits& traits,
                             unsigned int k=5) {
  if (begin1 == end1) return r_join(begin2, end2, oi, traits, k);
  General_polygon_set_2<Traits> gps(*begin1, traits);
  gps.join(std::next(begin1), end1, begin2, end2, k);
  return gps.polygons_with_holes(oi);
}

// Without traits
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline OutputIterator r_join(InputIterator1 begin1, InputIterator1 end1,
                             InputIterator2 begin2, InputIterator2 end2,
                             OutputIterator oi, unsigned int k=5) {
  typedef typename std::iterator_traits<InputIterator1>::value_type Pgn;
  typename Gps_polyline_traits<Pgn>::Traits traits;
  const typename Gps_polyline_traits<Pgn>::Polyline_traits& ptraits(traits);
  if (begin1 == end1) {
      if (begin2 == end2) return oi;
      return r_join(convert_polygon_iterator(begin2, ptraits),
          convert_polygon_iterator(end2, ptraits),
          convert_polygon_back(oi, *begin2), traits, k);
  }
  return r_join(convert_polygon_iterator(begin1, ptraits),
                convert_polygon_iterator(end1, ptraits),
                convert_polygon_iterator(begin2, ptraits),
                convert_polygon_iterator(end2, ptraits),
                convert_polygon_back(oi, *begin1), traits, k);
}

//@}

/// \name _difference() functions.
//@{

template <typename Pgn1, typename Pgn2, typename OutputIterator, typename Traits>
inline OutputIterator _difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                  OutputIterator oi, Traits& traits) {
  General_polygon_set_2<Traits> gps(pgn1, traits);
  gps.difference(pgn2);
  return gps.polygons_with_holes(oi);
}

template <typename Kernel, typename Container,
          typename Pgn1, typename Pgn2, typename OutputIterator>
inline OutputIterator _difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                  OutputIterator oi)
{
  // Use the first polygon to determine the (default) traits
  typedef typename Gps_polyline_traits<Pgn1>::Polyline_traits   Polyline_traits;

  typename Gps_polyline_traits<Pgn1>::Traits traits;
  const Polyline_traits& ptraits(traits);
  _difference(convert_polygon(pgn1, ptraits), convert_polygon(pgn2, ptraits),
              convert_polygon_back(oi, pgn1), traits);
  return oi;
}

//@}
/// \name s_symmetric_difference() and r_symmetric_difference() functions.
//@{

// Single
// With traits
template <typename Pgn1, typename Pgn2, typename OutputIterator, typename Traits>
inline OutputIterator s_symmetric_difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                             OutputIterator oi,
                                             Traits& traits) {
  General_polygon_set_2<Traits> gps(pgn1, traits);
  gps.symmetric_difference(pgn2);
  return gps.polygons_with_holes(oi);
}

// Without traits
template <typename Kernel, typename Container,
          typename Pgn1, typename Pgn2, typename OutputIterator>
inline OutputIterator s_symmetric_difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                             OutputIterator oi) {
  // Use the first polygon to determine the (default) traits
  typedef typename Gps_polyline_traits<Pgn1>::Polyline_traits   Polyline_traits;
  typename Gps_polyline_traits<Pgn1>::Traits traits;
  const Polyline_traits& ptraits(traits);
  s_symmetric_difference(convert_polygon(pgn1, ptraits),
                         convert_polygon(pgn2, ptraits),
                         convert_polygon_back(oi, pgn1), traits);
  return oi;
}

// Single Range
// With traits
template <typename InputIterator, typename OutputIterator, typename Traits>
inline
OutputIterator r_symmetric_difference(InputIterator begin, InputIterator end,
                                      OutputIterator oi, Traits& traits,
                                      unsigned int k=5) {
  if (begin == end) return (oi);
  General_polygon_set_2<Traits> gps(*begin, traits);
  gps.symmetric_difference(std::next(begin), end, k);
  return gps.polygons_with_holes(oi);
}

// Without traits
template <typename InputIterator, typename OutputIterator>
inline OutputIterator r_symmetric_difference(InputIterator begin,
                                             InputIterator end,
                                             OutputIterator oi,
                                             unsigned int k=5)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Pgn;
  typename Gps_polyline_traits<Pgn>::Traits traits;
  const typename Gps_polyline_traits<Pgn>::Polyline_traits& ptraits(traits);
  if (begin == end) return (oi);
  return r_symmetric_difference(convert_polygon_iterator(begin, ptraits),
                                convert_polygon_iterator(end, ptraits),
                                convert_polygon_back(oi, *begin), traits, k);
}

// Double Range
// With traits
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename Traits>
inline OutputIterator r_symmetric_difference(InputIterator1 begin1,
                                             InputIterator1 end1,
                                             InputIterator2 begin2,
                                             InputIterator2 end2,
                                             OutputIterator oi, Traits& traits,
                                             unsigned int k=5)
{
  if (begin1 == end1) return r_symmetric_difference(begin2, end2, oi, traits, k);
  General_polygon_set_2<Traits> gps(*begin1, traits);
  gps.symmetric_difference(std::next(begin1), end1, begin2, end2, k);
  return gps.polygons_with_holes(oi);
}

// Without traits
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline OutputIterator r_symmetric_difference(InputIterator1 begin1,
                                             InputIterator1 end1,
                                             InputIterator2 begin2,
                                             InputIterator2 end2,
                                             OutputIterator oi,
                                             unsigned int k=5) {
  typedef typename std::iterator_traits<InputIterator1>::value_type Pgn;
  typename Gps_polyline_traits<Pgn>::Traits traits;
  const typename Gps_polyline_traits<Pgn>::Polyline_traits& ptraits(traits);
  if (begin1 == end1){
      if (begin2 == end2) return oi;
      return r_symmetric_difference(convert_polygon_iterator(begin2, ptraits),
                                    convert_polygon_iterator(end2, ptraits),
                                    convert_polygon_back(oi, *begin2), traits, k);
  }
  return r_symmetric_difference(convert_polygon_iterator(begin1, ptraits),
                                convert_polygon_iterator(end1, ptraits),
                                convert_polygon_iterator(begin2, ptraits),
                                convert_polygon_iterator(end2, ptraits),
                                convert_polygon_back(oi, *begin1), traits, k);
}

//@}

/// \name _complement() functions.
//@{

// Compute the complemenet of a polygon
template <typename Kernel, typename Container, typename Traits>
void _complement(const Polygon_2<Kernel, Container>& pgn,
                 typename Traits::Polygon_with_holes_2& res, Traits& traits) {
  General_polygon_set_2<Traits> gps(pgn, traits);
  gps.complement();
  Oneset_iterator<typename Traits::Polygon_with_holes_2> oi(res);
  gps.polygons_with_holes(oi);
}

// Compute the complemenet of a general polygon
template <typename ArrTraits, typename Traits>
void _complement(const General_polygon_2<ArrTraits>& pgn,
                 typename Traits::Polygon_with_holes_2& res, Traits& traits) {
  General_polygon_set_2<Traits> gps(pgn, traits);
  gps.complement();
  Oneset_iterator<typename Traits::Polygon_with_holes_2> oi(res);
  gps.polygons_with_holes(oi);
}

// Compute the complemenet of a polygon with holes
template <typename Kernel, typename Container, typename OutputIterator,
          typename Traits>
OutputIterator _complement(const Polygon_with_holes_2<Kernel, Container>& pgn,
                           OutputIterator oi, Traits& traits) {
  General_polygon_set_2<Traits> gps(pgn, traits);
  gps.complement();
  return gps.polygons_with_holes(oi);
}

// Compute the complemenet of a general polygon with holes
template <typename Pgn, typename OutputIterator, typename Traits>
OutputIterator _complement(const General_polygon_with_holes_2<Pgn>& pgn,
                           OutputIterator oi, Traits& traits) {
  General_polygon_set_2<Traits> gps(pgn, traits);
  gps.complement();
  return gps.polygons_with_holes(oi);
}

// Compute the complemenet of a polygon
template <typename Kernel, typename Container, typename Pwh>
void _complement(const Polygon_2<Kernel, Container>& pgn, Pwh& pwh) {
  // Use the polygon to determine the (default) traits
  typedef Polygon_2<Kernel, Container>                          Pgn;
  typedef typename Gps_polyline_traits<Pgn>::Polyline_traits    Polyline_traits;
  typedef General_polygon_2<Polyline_traits>                    General_pgn;
  typedef General_polygon_with_holes_2<General_pgn>             General_pwh;

  General_pwh general_pwh;
  typename Gps_polyline_traits<Pgn>::Traits traits;
  const Polyline_traits& ptraits(traits);
  _complement(convert_polygon(pgn, ptraits), general_pwh, traits);
  pwh = convert_polygon_back<Kernel, Container>(general_pwh);
}

// Compute the complemenet of a polygon with holes
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator _complement(const Polygon_with_holes_2<Kernel, Container>& pgn,
                           OutputIterator oi) {
  // Use the polygon with holes to determine the (default) traits
  typedef Polygon_with_holes_2<Kernel, Container>               Pgn;
  typedef typename Gps_polyline_traits<Pgn>::Polyline_traits    Polyline_traits;

  typename Gps_polyline_traits<Pgn>::Traits traits;
  const Polyline_traits& ptraits(traits);
  complement(convert_polygon(pgn, ptraits), convert_polygon_back(oi, pgn),
             traits);
  return oi;
}

//@}

} //namespace CGAL

#endif
