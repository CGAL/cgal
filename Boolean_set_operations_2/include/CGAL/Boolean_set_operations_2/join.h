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

#ifndef CGAL_BOOLEAN_SET_OPERATIONS_2_JOIN_H
#define CGAL_BOOLEAN_SET_OPERATIONS_2_JOIN_H

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

namespace CGAL
{

/// \name join() functions.
//@{

// Polygon_2, Polygon_2 ========================================================
// With traits
template <typename Kernel, typename Container, typename Traits>
inline bool join(const Polygon_2<Kernel, Container>& pgn1,
                 const Polygon_2<Kernel, Container>& pgn2,
                 Polygon_with_holes_2<Kernel, Container>& res, Traits& traits)
{ return s_join(pgn1, pgn2, res, traits); }

// With Tag_true
template <typename Kernel, typename Container>
inline bool join(const Polygon_2<Kernel, Container>& pgn1,
                 const Polygon_2<Kernel, Container>& pgn2,
                 Polygon_with_holes_2<Kernel, Container>& res,
                 Tag_true = Tag_true())
{ return s_join<Kernel, Container>(pgn1, pgn2, res); }

// With Tag_false
template <typename Kernel, typename Container>
inline bool join(const Polygon_2<Kernel, Container>& pgn1,
                 const Polygon_2<Kernel, Container>& pgn2,
                 Polygon_with_holes_2<Kernel, Container>& res,
                 Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_join(pgn1, pgn2, res, traits);
}

// Polygon_2, Polygon_with_holes_2 =============================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
inline bool join(const Polygon_2<Kernel, Container>& pgn1,
                 const Polygon_with_holes_2<Kernel, Container>& pgn2,
                 Polygon_with_holes_2<Kernel, Container>& res, Traits& traits)
{ return s_join(pgn1, pgn2, res, traits); }

// With Tag_true
template <typename Kernel, typename Container>
inline bool join(const Polygon_2<Kernel, Container>& pgn1,
                 const Polygon_with_holes_2<Kernel, Container>& pgn2,
                 Polygon_with_holes_2<Kernel, Container>& res,
                 Tag_true = Tag_true())
{ return s_join<Kernel, Container>(pgn1, pgn2, res); }

// With Tag_false
template <typename Kernel, typename Container>
inline bool join(const Polygon_2<Kernel, Container>& pgn1,
                 const Polygon_with_holes_2<Kernel, Container>& pgn2,
                 Polygon_with_holes_2<Kernel, Container>& res,
                 Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_join(pgn1, pgn2, res, traits);
}

// Polygon_with_holes_2, Polygon_2 =============================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
inline bool join(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                 const Polygon_2<Kernel, Container>& pgn2,
                 Polygon_with_holes_2<Kernel, Container>& res, Traits& traits)
{ return s_join(pgn1, pgn2, res, traits); }

// With Tag_true
template <typename Kernel, typename Container>
inline bool join(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                 const Polygon_2<Kernel, Container>& pgn2,
                 Polygon_with_holes_2<Kernel, Container>& res,
                 Tag_true = Tag_true())
{ return s_join<Kernel, Container>(pgn1, pgn2, res); }

// With Tag_false
template <typename Kernel, typename Container>
inline bool join(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                 const Polygon_2<Kernel, Container>& pgn2,
                 Polygon_with_holes_2<Kernel, Container>& res,
                 Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_with_holes_2<Kernel, Container>       Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_join(pgn1, pgn2, res, traits);
}

// Polygon_with_holes_2, Polygon_with_holes_2 ==================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
inline bool join(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                 const Polygon_with_holes_2<Kernel, Container>& pgn2,
                 Polygon_with_holes_2<Kernel, Container>& res, Traits& traits)
{ return s_join(pgn1, pgn2, res, traits); }

// With Tag_true
template <typename Kernel, typename Container>
inline bool join(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                 const Polygon_with_holes_2<Kernel, Container>& pgn2,
                 Polygon_with_holes_2<Kernel, Container>& res,
                 Tag_true = Tag_true())
{ return s_join<Kernel, Container>(pgn1, pgn2, res); }

// With Tag_false
template <typename Kernel, typename Container>
inline bool join(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                 const Polygon_with_holes_2<Kernel, Container>& pgn2,
                 Polygon_with_holes_2<Kernel, Container>& res,
                 Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_with_holes_2<Kernel, Container>       Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_join(pgn1, pgn2, res, traits);
}

// General_polygon_2, General_polygon_2 ========================================
// With Traits
template <typename ArrTraits, typename Traits>
inline bool
join(const General_polygon_2<ArrTraits>& pgn1,
     const General_polygon_2<ArrTraits>& pgn2,
     General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res,
     Traits& traits)
{ return s_join(pgn1, pgn2, res, traits); }

// Without Traits
template <typename ArrTraits>
inline bool
join(const General_polygon_2<ArrTraits>& pgn1,
     const General_polygon_2<ArrTraits>& pgn2,
     General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_join(pgn1, pgn2, res, traits);
}

// General_polygon_2, General_polygon_with_holes_2 =============================
// With Traits
template <typename ArrTraits, typename Traits>
inline bool
join(const General_polygon_2<ArrTraits>& pgn1,
     const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
     General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res,
     Traits& traits)
{ return s_join(pgn1, pgn2, res, traits); }

// Without Traits
template <typename ArrTraits>
inline bool
join(const General_polygon_2<ArrTraits>& pgn1,
     const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn2,
     General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_join(pgn1, pgn2, res, traits);
}

// General_polygon_with_holes_2, General_polygon_2 =============================
// With Traits
template <typename ArrTraits, typename Traits>
inline bool
join(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
     const General_polygon_2<ArrTraits>& pgn2,
     General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res,
     Traits& traits)
{ return s_join(pgn1, pgn2, res, traits); }

// Without Traits
template <typename ArrTraits>
inline bool
join(const General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& pgn1,
     const General_polygon_2<ArrTraits>& pgn2,
     General_polygon_with_holes_2<General_polygon_2<ArrTraits> >& res)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typedef General_polygon_with_holes_2<Polygon>         Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_join(pgn1, pgn2, res, traits);
}

// General_polygon_with_holes_2, General_polygon_with_holes_2 ==================
// With Traits
template <typename Polygon_, typename Traits>
inline bool join(const General_polygon_with_holes_2<Polygon_>& pgn1,
                 const General_polygon_with_holes_2<Polygon_>& pgn2,
                 General_polygon_with_holes_2<Polygon_>& res, Traits& traits)
{ return s_join(pgn1, pgn2, res, traits); }

// Without Traits
template <typename Polygon_>
inline bool join(const General_polygon_with_holes_2<Polygon_>& pgn1,
                 const General_polygon_with_holes_2<Polygon_>& pgn2,
                 General_polygon_with_holes_2<Polygon_>& res)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_with_holes_2<Polygon_>        Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_join(pgn1, pgn2, res, traits);
}

//@}

/// \name Aggregated join() functions.
//@{

// With Traits
template <typename InputIterator, typename OutputIterator, typename Traits>
inline OutputIterator join(InputIterator begin, InputIterator end,
                           OutputIterator oi, Traits& traits, unsigned int k=5)
{ return r_join(begin, end, oi, traits, k); }

// Without Traits
// Tag_true => convert to polylines
template <typename InputIterator, typename OutputIterator>
inline OutputIterator join(InputIterator begin, InputIterator end,
                           OutputIterator oi, Tag_true = Tag_true(),
                           unsigned int k=5,
                           Enable_if_Polygon_2_iterator<InputIterator>* = 0)
{ return r_join(begin, end, oi, k); }

// Tag_false => do not convert to polylines
template <typename InputIterator, typename OutputIterator>
inline OutputIterator join(InputIterator begin, InputIterator end,
                           OutputIterator oi, Tag_false, unsigned int k=5,
                           Enable_if_Polygon_2_iterator<InputIterator>* = 0)
{
  typename Iterator_to_gps_traits<InputIterator>::Traits traits;
  return r_join(begin, end, oi, traits, k);
}

// General polygons or polygons with holes
template <typename InputIterator, typename OutputIterator>
inline OutputIterator join(InputIterator begin, InputIterator end,
                           OutputIterator oi, unsigned int k=5,
                           Disable_if_Polygon_2_iterator<InputIterator>* = 0)
{
  typename Iterator_to_gps_traits<InputIterator>::Traits traits;
  return r_join(begin, end, oi, traits, k);
}

// Join two ranges of simple polygons and polygons with holes.
// With Traits
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename Traits>
inline OutputIterator join(InputIterator1 begin1, InputIterator1 end1,
                           InputIterator2 begin2, InputIterator2 end2,
                           OutputIterator oi, Traits& traits, unsigned int k=5)
{ return r_join(begin1, end1, begin2, end2, oi, traits, k); }

// Without Traits
// Tag_true => convert to polylines
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline OutputIterator join(InputIterator1 begin1, InputIterator1 end1,
                           InputIterator2 begin2, InputIterator2 end2,
                           OutputIterator oi, Tag_true = Tag_true(),
                           unsigned int k=5,
                           Enable_if_Polygon_2_iterator<InputIterator1>* = 0)
{ return r_join(begin1, end1, begin2, end2, oi, k); }

// Tag_false => do not convert to polylines
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline OutputIterator join(InputIterator1 begin1, InputIterator1 end1,
                           InputIterator2 begin2, InputIterator2 end2,
                           OutputIterator oi, Tag_false, unsigned int k=5,
                           Enable_if_Polygon_2_iterator<InputIterator1>* = 0)
{
  typename Iterator_to_gps_traits<InputIterator1>::Traits traits;
  return r_join(begin1, end1, begin2, end2, oi, traits, k);
}

// General polygons or polygons with holes
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline OutputIterator join(InputIterator1 begin1, InputIterator1 end1,
                           InputIterator2 begin2, InputIterator2 end2,
                           OutputIterator oi, unsigned int k=5,
                           Disable_if_Polygon_2_iterator<InputIterator1>* = 0)
{
  typename Iterator_to_gps_traits<InputIterator1>::Traits traits;
  return r_join(begin1, end1, begin2, end2, oi, traits, k);
}

//@}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
