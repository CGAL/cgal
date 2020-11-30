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

/// \name _do_intersect() functions.
//@{

// With Traits
template <typename Pgn1, class Pgn2, typename Traits>
inline bool _do_intersect(const Pgn1& pgn1, const Pgn2& pgn2, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn1);
  return gps.do_intersect(pgn2);
}

// Without Traits
template <typename Pgn1, typename Pgn2>
inline bool _do_intersect(const Pgn1& pgn1, const Pgn2& pgn2) {
  typedef typename Gps_default_traits<Pgn1>::Arr_traits         Segment_traits;
  typedef Arr_polyline_traits_2<Segment_traits>                 Polyline_traits;
  typedef Gps_traits_2<Polyline_traits>                         Traits;
  Traits traits;
  return _do_intersect(convert_polygon<Polyline_traits>(pgn1),
                       convert_polygon<Polyline_traits>(pgn2), traits);
}

//@}

/// \name _oriented_side() functions.
//@{

// With Traits
template <typename Obj, typename Pgn, typename Traits>
inline Oriented_side _oriented_side(const Obj& obj, const Pgn& pgn, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn);
  return gps.oriented_side(obj);
}

// Without Traits
template <typename Pgn1, typename Pgn2>
inline Oriented_side _oriented_side(const Pgn1& pgn1, const Pgn2& pgn2)
{
  // Use the first polygon to determine the (default) traits
  typedef typename Gps_default_traits<Pgn1>::Arr_traits         Segment_traits;
  typedef Arr_polyline_traits_2<Segment_traits>                 Polyline_traits;
  typedef Gps_traits_2<Polyline_traits>                         Traits;
  Traits traits;
  return _oriented_side(convert_polygon<Polyline_traits>(pgn1),
                        convert_polygon<Polyline_traits>(pgn2), traits);
}

// Without Traits
template <typename Kernel, typename Pgn>
inline Oriented_side _oriented_side(const typename Kernel::Point_2& point,
                                    const Pgn& pgn)
{
  // Use the polygon to determine the (default) traits
  typedef typename Gps_default_traits<Pgn>::Arr_traits          Segment_traits;
  typedef Arr_polyline_traits_2<Segment_traits>                 Polyline_traits;
  typedef Gps_traits_2<Polyline_traits>                         Traits;
  Traits traits;
  return _oriented_side(point, convert_polygon<Polyline_traits>(pgn), traits);
}

//@}

/// \name _intersection() functions.
//@{

// With Traits
template <typename Pgn1, typename Pgn2, typename OutputIterator, typename Traits>
inline OutputIterator _intersection(const Pgn1& pgn1, const Pgn2& pgn2,
                                    OutputIterator oi, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn1);
  gps.intersection(pgn2);
  return gps.polygons_with_holes(oi);
}

// Without Traits
  template <typename Kernel, typename Container,
            typename Pgn1, typename Pgn2, typename OutputIterator>
inline OutputIterator _intersection(const Pgn1& pgn1, const Pgn2& pgn2,
                                    OutputIterator oi)
{
  // Use the first polygon to determine the (default) traits
  typedef typename Gps_default_traits<Pgn1>::Arr_traits         Segment_traits;
  typedef Arr_polyline_traits_2<Segment_traits>                 Polyline_traits;
  typedef Gps_traits_2<Polyline_traits>                         Traits;
  Traits traits;

  typedef General_polygon_2<Polyline_traits>                    General_pgn;
  typedef General_polygon_with_holes_2<General_pgn>             General_pwh;
  std::list<General_pwh> general_pwhs;

  _intersection(convert_polygon<Polyline_traits>(pgn1),
                convert_polygon<Polyline_traits>(pgn2),
               std::back_inserter(general_pwhs), traits);
  for (const auto& general_pwh : general_pwhs)
    *oi++ = convert_polygon_back<Kernel, Container>(general_pwh);
  return oi;
}

//@}

/// \name _join() functions.
//@{

// Polygon_2
template <typename Traits>
inline bool _is_empty(const typename Traits::Polygon_2& pgn, Traits& tr)
{
  typedef typename Traits::Curve_const_iterator Curve_const_iterator;
  const std::pair<Curve_const_iterator, Curve_const_iterator>& itr_pair =
    tr.construct_curves_2_object()(pgn);
  return (itr_pair.first == itr_pair.second);
}

// A polygon with holes cannot be empty.
template <typename Traits>
inline bool _is_empty(const typename Traits::Polygon_with_holes_2&, Traits&)
{ return false; }

// With Traits
template <typename Pgn1, typename Pgn2, typename Traits>
inline bool _join(const Pgn1& pgn1, const Pgn2& pgn2,
                  typename Traits::Polygon_with_holes_2& res, Traits& tr)
{
  if (_is_empty(pgn1, tr) || _is_empty(pgn2, tr)) return false;

  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn1);
  gps.join(pgn2);
  if (gps.number_of_polygons_with_holes() == 1) {
    Oneset_iterator<typename Traits::Polygon_with_holes_2> oi(res);
    gps.polygons_with_holes(oi);
    return true;
  }

  // the polygon doesnt intersect, the original pgn1, pgn2 contain the union
  return false;
}

// Without Traits
template <typename Kernel, typename Container,
          typename Pgn1, typename Pgn2, typename Pwh>
inline bool _join(const Pgn1& pgn1, const Pgn2& pgn2, Pwh& pwh)
{
  // Use the first polygon to determine the (default) traits
  typedef typename Gps_default_traits<Pgn1>::Arr_traits         Segment_traits;
  typedef Arr_polyline_traits_2<Segment_traits>                 Polyline_traits;
  typedef Gps_traits_2<Polyline_traits>                         Traits;
  Traits traits;

  typedef General_polygon_2<Polyline_traits>                    General_pgn;
  typedef General_polygon_with_holes_2<General_pgn>             General_pwh;
  General_pwh general_pwh;

  auto res = _join(convert_polygon<Polyline_traits>(pgn1),
                   convert_polygon<Polyline_traits>(pgn2),
                   general_pwh, traits);
  pwh = convert_polygon_back<Kernel, Container>(general_pwh);
  return res;
}

//@}

/// \name _difference() functions.
//@{

template <typename Pgn1, typename Pgn2, typename OutputIterator, typename Traits>
inline OutputIterator _difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                  OutputIterator oi, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn1);
  gps.difference(pgn2);
  return gps.polygons_with_holes(oi);
}

template <typename Kernel, typename Container,
          typename Pgn1, typename Pgn2, typename OutputIterator>
inline OutputIterator _difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                  OutputIterator oi)
{
  // Use the first polygon to determine the (default) traits
  typedef typename Gps_default_traits<Pgn1>::Arr_traits         Segment_traits;
  typedef Arr_polyline_traits_2<Segment_traits>                 Polyline_traits;
  typedef Gps_traits_2<Polyline_traits>                         Traits;
  Traits traits;

  typedef General_polygon_2<Polyline_traits>                    General_pgn;
  typedef General_polygon_with_holes_2<General_pgn>             General_pwh;
  std::list<General_pwh> general_pwhs;

  _difference(convert_polygon<Polyline_traits>(pgn1),
              convert_polygon<Polyline_traits>(pgn2),
              std::back_inserter(general_pwhs), traits);
  for (const auto& general_pwh : general_pwhs)
    *oi++ = convert_polygon_back<Kernel, Container>(general_pwh);
  return oi;
}

//@}
/// \name _symmetric_difference() functions.
//@{

template <typename Pgn1, typename Pgn2, typename OutputIterator, typename Traits>
inline OutputIterator _symmetric_difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                            OutputIterator oi, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn1);
  gps.symmetric_difference(pgn2);
  return gps.polygons_with_holes(oi);
}

template <typename Kernel, typename Container,
          typename Pgn1, typename Pgn2, typename OutputIterator>
inline OutputIterator _symmetric_difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                            OutputIterator oi)
{
  // Use the first polygon to determine the (default) traits
  typedef typename Gps_default_traits<Pgn1>::Arr_traits         Segment_traits;
  typedef Arr_polyline_traits_2<Segment_traits>                 Polyline_traits;
  typedef Gps_traits_2<Polyline_traits>                         Traits;
  Traits traits;

  typedef General_polygon_2<Polyline_traits>                    General_pgn;
  typedef General_polygon_with_holes_2<General_pgn>             General_pwh;
  std::list<General_pwh> general_pwhs;

  _symmetric_difference(convert_polygon<Polyline_traits>(pgn1),
                        convert_polygon<Polyline_traits>(pgn2),
              std::back_inserter(general_pwhs), traits);
  for (const auto& general_pwh : general_pwhs)
    *oi++ = convert_polygon_back<Kernel, Container>(general_pwh);
  return oi;
}

//@}

/// \name _complement() functions.
//@{

// Compute the complemenet of a (general) polygon
template <typename Pgn, typename Traits>
void _complement(const Pgn& pgn, typename Traits::Polygon_with_holes_2& res,
                 Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn);
  gps.complement();
  Oneset_iterator<typename Traits::Polygon_with_holes_2> oi(res);
  gps.polygons_with_holes(oi);
}

// Compute the complemenet of a (general) polygon
template <typename Kernel, typename Container, typename Pgn, typename Pwh>
void _complement(const Pgn& pgn, Pwh& pwh)
{
  // Use the first polygon to determine the (default) traits
  typedef typename Gps_default_traits<Pgn>::Arr_traits         Segment_traits;
  typedef Arr_polyline_traits_2<Segment_traits>                 Polyline_traits;
  typedef Gps_traits_2<Polyline_traits>                         Traits;
  Traits traits;

  typedef General_polygon_2<Polyline_traits>                    General_pgn;
  typedef General_polygon_with_holes_2<General_pgn>             General_pwh;
  General_pwh general_pwh;

  _complement(convert_polygon<Polyline_traits>(pgn), general_pwh, traits);
  pwh = convert_polygon_back<Kernel, Container>(general_pwh);
}

// Compute the complemenet of a (general) polygon with holes
template <typename Pgn, typename OutputIterator, typename Traits>
OutputIterator _complement(const Pgn& pgn, OutputIterator oi, Traits& traits)
{
  General_polygon_set_2<Traits> gps(traits);
  gps.insert(pgn, traits);
  gps.complement();
  return gps.polygons_with_holes(oi, traits);
}

// Compute the complemenet of a (general) polygon with holes
template <typename Kernel, typename Container,
          typename Pgn, typename OutputIterator>
OutputIterator _complement(const Pgn& pgn, OutputIterator oi)
{
  // Use the first polygon to determine the (default) traits
  typedef typename Gps_default_traits<Pgn>::Arr_traits          Segment_traits;
  typedef Arr_polyline_traits_2<Segment_traits>                 Polyline_traits;
  typedef Gps_traits_2<Polyline_traits>                         Traits;
  Traits traits;

  typedef General_polygon_2<Polyline_traits>                    General_pgn;
  typedef General_polygon_with_holes_2<General_pgn>             General_pwh;
  std::list<General_pwh> general_pwhs;

  complement(convert_polygon<Polyline_traits>(pgn),
             std::back_inserter(general_pwhs), traits);
  for (const auto& general_pwh : general_pwhs)
    *oi++ = convert_polygon_back<Kernel, Container>(general_pwh);
  return oi;
}

//@}

} //namespace CGAL

#endif
