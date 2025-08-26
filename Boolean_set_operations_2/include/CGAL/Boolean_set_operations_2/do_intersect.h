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

#ifndef CGAL_BOOLEAN_SET_OPERATIONS_2_DO_INTERSECT_H
#define CGAL_BOOLEAN_SET_OPERATIONS_2_DO_INTERSECT_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/iterator.h>
#include <CGAL/Boolean_set_operations_2/Bso_internal_functions.h>
#include <CGAL/Boolean_set_operations_2/Polygon_conversions.h>
#include <CGAL/type_traits/is_iterator.h>

namespace CGAL {

/// \name do_intersect() functions.
//@{

/*! We do not use polyline for do_intersect), as we relly on the overlay traits
 * to only intercept intersections between the interiors of segments that
 * comprise the boundary of polygons. Observe that The intersections between the
 * interiors of polylines that comprise the boundary of polygons may include an
 * endpoint of a segment, and we do not want that.
 */

// Polygon_2, Polygon_2 ========================================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2,
                         Traits& traits)
{ return s_do_intersect(pgn1, pgn2, traits); }

// without traits
template <typename Kernel, typename Container>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2) {
  using Polygon = Polygon_2<Kernel, Container>;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_do_intersect(pgn1, pgn2, traits);
}

// Polygon_2, Polygon_with_hole_2 ==============================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2,
                         Traits& traits)
{ return s_do_intersect(pgn1, pgn2, traits); }

// Without traits
template <typename Kernel, typename Container>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2) {
  // Use the first polygon to determine the (default) traits
  using Polygon = Polygon_2<Kernel, Container>;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_do_intersect(pgn1, pgn2, traits);
}

// Polygon_with_hole_2, Polygon_2 ==============================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2,
                         Traits& traits)
{ return s_do_intersect(pgn1, pgn2, traits); }

// Without traits
template <typename Kernel, typename Container>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2) {
  // Use the first polygon to determine the (default) traits
  using Polygon_with_holes = Polygon_with_holes_2<Kernel, Container>;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_do_intersect(pgn1, pgn2, traits);
}

// Polygon_with_hole_2, Polygon_with_hole_2 ====================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2,
                         Traits& traits)
{ return s_do_intersect(pgn1, pgn2, traits); }

// Without traits
template <typename Kernel, typename Container>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2) {
  // Use the first polygon to determine the (default) traits
  using Polygon_with_holes = Polygon_with_holes_2<Kernel, Container>;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_do_intersect(pgn1, pgn2, traits);
}

// General_polygon_2, General_polygon_2 ========================================
// With Traits
template <typename ArrTraits, typename GpsTraits>
inline bool do_intersect(const General_polygon_2<ArrTraits>& pgn1,
                         const General_polygon_2<ArrTraits>& pgn2,
                         GpsTraits& traits)
{ return s_do_intersect(pgn1, pgn2, traits); }

// Without Traits
template <typename ArrTraits>
inline bool do_intersect(const General_polygon_2<ArrTraits>& pgn1,
                         const General_polygon_2<ArrTraits>& pgn2) {
  // Use the first polygon to determine the (default) traits
  using Polygon = General_polygon_2<ArrTraits>;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_do_intersect(pgn1, pgn2, traits);
}

// General_polygon_2, General_polygon_with_holes_2 =============================
// With Traits
template <typename ArrTraits, typename GpsTraits>
inline bool do_intersect(const General_polygon_2<ArrTraits>& pgn1,
                         const General_polygon_with_holes_2
                           <General_polygon_2<ArrTraits> >& pgn2,
                         GpsTraits& traits)
{ return s_do_intersect(pgn1, pgn2, traits); }

// Without Traits
template <typename ArrTraits>
inline bool do_intersect(const General_polygon_2<ArrTraits>& pgn1,
                         const General_polygon_with_holes_2
                           <General_polygon_2<ArrTraits> >& pgn2) {
  // Use the first polygon to determine the (default) traits
  using Polygon = General_polygon_2<ArrTraits>;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_do_intersect(pgn1, pgn2, traits);
}

// General_polygon_with_holes_2, General_polygon_2 =============================
// With Traits
template <typename ArrTraits, typename GpsTraits>
inline bool do_intersect(const General_polygon_with_holes_2
                           <General_polygon_2<ArrTraits> >& pgn1,
                         const General_polygon_2<ArrTraits>& pgn2,
                         GpsTraits& traits)
{ return s_do_intersect(pgn1, pgn2, traits); }

// Without Traits
template <typename ArrTraits>
inline bool do_intersect(const General_polygon_with_holes_2
                           <General_polygon_2<ArrTraits> >& pgn1,
                         const General_polygon_2<ArrTraits>& pgn2) {
  // Use the first polygon to determine the (default) traits
  using Polygon = General_polygon_2<ArrTraits>;
  using Polygon_with_holes = General_polygon_with_holes_2<Polygon>;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_do_intersect(pgn1, pgn2, traits);
}

// General_polygon_with_holes_2, General_polygon_with_holes_2 ==================
// With Traits
template <typename Polygon_, typename Traits>
inline bool do_intersect(const General_polygon_with_holes_2<Polygon_>& pgn1,
                         const General_polygon_with_holes_2<Polygon_>& pgn2,
                         Traits& traits)
{ return s_do_intersect(pgn1, pgn2, traits); }

// Without Traits
template <typename Polygon_>
inline bool do_intersect(const General_polygon_with_holes_2<Polygon_>& pgn1,
                         const General_polygon_with_holes_2<Polygon_>& pgn2) {
  // Use the first polygon to determine the (default) traits
  using Polygon_with_holes = General_polygon_with_holes_2<Polygon_>;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_do_intersect(pgn1, pgn2, traits);
}

//@}

/// \name Aggregated do_intersect() functions.
//@{

// With Traits
template <typename InputIterator, typename Traits>
inline bool do_intersect(InputIterator begin, InputIterator end, Traits& traits,
                         unsigned int k = 5,
                         std::enable_if_t<CGAL::is_iterator<InputIterator>::value>* = 0)
{ return r_do_intersect(begin, end, traits, k); }

// Without Traits
template <typename InputIterator>
inline bool do_intersect(InputIterator begin, InputIterator end, unsigned int k = 5,
                         std::enable_if_t<CGAL::is_iterator<InputIterator>::value>* = 0,
                         Enable_if_Polygon_2_iterator<InputIterator>* = 0) {
  typename Iterator_to_gps_traits<InputIterator>::Traits traits;
  return r_do_intersect(begin, end, traits, k);
}

// General polygons or polygons with holes
template <typename InputIterator>
inline bool do_intersect(InputIterator begin, InputIterator end,
                         unsigned int k = 5,
                         std::enable_if_t<CGAL::is_iterator<InputIterator>::value>* = 0,
                         Disable_if_Polygon_2_iterator<InputIterator>* = 0) {
  typename Iterator_to_gps_traits<InputIterator>::Traits traits;
  return do_intersect(begin, end, traits, k);
}

// With Traits
template <typename InputIterator1, typename InputIterator2, typename Traits>
inline bool do_intersect(InputIterator1 begin1, InputIterator1 end1,
                         InputIterator2 begin2, InputIterator2 end2,
                         Traits& traits, unsigned int k = 5)
{ return r_do_intersect(begin1, end1, begin2, end2, traits, k); }

// Without Traits
template <typename InputIterator1, typename InputIterator2>
inline bool do_intersect(InputIterator1 begin1, InputIterator1 end1,
                         InputIterator2 begin2, InputIterator2 end2,
                         unsigned int k = 5,
                         Enable_if_Polygon_2_iterator<InputIterator1>* = 0)
{ return r_do_intersect(begin1, end1, begin2, end2, k); }

// General polygons or polygons with holes
template <typename InputIterator1, typename InputIterator2>
inline bool do_intersect(InputIterator1 begin1, InputIterator1 end1,
                         InputIterator2 begin2, InputIterator2 end2,
                         unsigned int k = 5,
                         Disable_if_Polygon_2_iterator<InputIterator1>* = 0) {
  typename Iterator_to_gps_traits<InputIterator1>::Traits traits;
  return r_do_intersect(begin1, end1, begin2, end2, traits, k);
}

//@}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
