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

#ifndef CGAL_BOOLEAN_SET_OPERATIONS_H
#define CGAL_BOOLEAN_SET_OPERATIONS_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <boost/utility/enable_if.hpp>

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
#include <CGAL/is_iterator.h>

namespace CGAL {

/// \name do_intersect() functions.
//@{

// Polygon_2, Polygon_2 ========================================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2,
                         Traits& traits)
{ return s_do_intersect(pgn1, pgn2, traits); }

// With Tag_true
template <typename Kernel, typename Container>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2,
                         Tag_true = Tag_true())
{ return s_do_intersect(pgn1, pgn2); }

// With Tag_false
template <typename Kernel, typename Container>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2,
                         Tag_false)
{
  typedef Polygon_2<Kernel, Container>                          Polygon;
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

// With Tag_true
template <typename Kernel, typename Container>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2,
                         Tag_true = Tag_true())
{ return s_do_intersect(pgn1, pgn2); }

// With Tag_false
template <typename Kernel, typename Container>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2,
                         Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_2<Kernel, Container>                          Polygon;
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

// With Tag_true
template <typename Kernel, typename Container>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2,
                         Tag_true = Tag_true())
{ return s_do_intersect(pgn1, pgn2); }

// With Tag_false
template <typename Kernel, typename Container>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2,
                         Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_with_holes_2<Kernel, Container>       Polygon_with_holes;
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

// With Tag_true
template <typename Kernel, typename Container>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2,
                         Tag_true = Tag_true())
{ return s_do_intersect(pgn1, pgn2); }

// With Tag_false
template <typename Kernel, typename Container>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2,
                         Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_with_holes_2<Kernel, Container>       Polygon_with_holes;
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
                         const General_polygon_2<ArrTraits>& pgn2)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
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
                           <General_polygon_2<ArrTraits> >& pgn2)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
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
                         const General_polygon_2<ArrTraits>& pgn2)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typedef General_polygon_with_holes_2<Polygon>         Polygon_with_holes;
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
                         const General_polygon_with_holes_2<Polygon_>& pgn2)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_with_holes_2<Polygon_>        Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_do_intersect(pgn1, pgn2, traits);
}

//@}

/// \name intersection() functions.
//@{

// Polygon_2, Polygon_2 ========================================================
// With Traits
template <typename Kernel, typename Container, typename OutputIterator,
          typename Traits>
inline OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1,
                                   const Polygon_2<Kernel, Container>& pgn2,
                                   OutputIterator out, Traits& traits)
{ return s_intersection(pgn1, pgn2, out, traits); }

// With Tag_true
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1,
                                   const Polygon_2<Kernel, Container>& pgn2,
                                   OutputIterator out, Tag_true = Tag_true())
{ return s_intersection<Kernel, Container>(pgn1, pgn2, out); }

// With Tag_false
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1,
                                   const Polygon_2<Kernel, Container>& pgn2,
                                   OutputIterator out, Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_intersection(pgn1, pgn2, out, traits);
}

// Polygon_2, Polygon_with_holes_2 =============================================
// With Traits
template <typename Kernel, typename Container, typename OutputIterator,
          typename Traits>
inline OutputIterator
intersection(const Polygon_2<Kernel, Container>& pgn1,
             const Polygon_with_holes_2<Kernel, Container>& pgn2,
             OutputIterator out, Traits& traits)
{ return s_intersection(pgn1, pgn2, out, traits); }

// With Tag_true
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
intersection(const Polygon_2<Kernel, Container>& pgn1,
             const Polygon_with_holes_2<Kernel, Container>& pgn2,
             OutputIterator out, Tag_true = Tag_true())
{ return s_intersection<Kernel, Container>(pgn1, pgn2, out); }

// With Tag_false
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
intersection(const Polygon_2<Kernel, Container>& pgn1,
             const Polygon_with_holes_2<Kernel, Container>& pgn2,
             OutputIterator out, Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_intersection(pgn1, pgn2, out, traits);
}

// Polygon_with_holes_2, Polygon_2 =============================================
// With Traits
template <typename Kernel, typename Container, typename OutputIterator,
          typename Traits>
inline OutputIterator
intersection(const Polygon_with_holes_2<Kernel, Container>& pgn1,
             const Polygon_2<Kernel, Container>& pgn2,
             OutputIterator out, Traits& traits)
{ return s_intersection(pgn1, pgn2, out, traits); }

// With Tag_true
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
intersection(const Polygon_with_holes_2<Kernel, Container>& pgn1,
             const Polygon_2<Kernel, Container>& pgn2,
             OutputIterator out, Tag_true = Tag_true())
{ return s_intersection<Kernel, Container>(pgn1, pgn2, out); }

// With Tag_false
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
intersection(const Polygon_with_holes_2<Kernel, Container>& pgn1,
             const Polygon_2<Kernel, Container>& pgn2,
             OutputIterator out, Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_with_holes_2<Kernel, Container>       Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_intersection(pgn1, pgn2, out, traits);
}

// Polygon_with_holes_2, Polygon_with_holes_2 ==================================
// With Traits
template <typename Kernel, typename Container, typename OutputIterator,
          typename Traits>
inline OutputIterator
intersection(const Polygon_with_holes_2<Kernel, Container>& pgn1,
             const Polygon_with_holes_2<Kernel, Container>& pgn2,
             OutputIterator out, Traits& traits)
{ return s_intersection(pgn1, pgn2, out, traits); }

// With Tag_true
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
intersection(const Polygon_with_holes_2<Kernel, Container>& pgn1,
             const Polygon_with_holes_2<Kernel, Container>& pgn2,
             OutputIterator out, Tag_true = Tag_true())
{ return s_intersection<Kernel, Container>(pgn1, pgn2, out); }

// With Tag_false
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
intersection(const Polygon_with_holes_2<Kernel, Container>& pgn1,
             const Polygon_with_holes_2<Kernel, Container>& pgn2,
             OutputIterator out, Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_with_holes_2<Kernel, Container>       Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_intersection(pgn1, pgn2, out, traits);
}

// General_polygon_2, General_polygon_2 ========================================
// With Traits
template <typename ArrTraits, typename OutputIterator, typename Traits>
inline OutputIterator intersection(const General_polygon_2<ArrTraits>& pgn1,
                                   const General_polygon_2<ArrTraits>& pgn2,
                                   OutputIterator out, Traits& traits)
{ return s_intersection(pgn1, pgn2, out, traits); }

// Without Traits
template <typename ArrTraits, typename OutputIterator>
inline OutputIterator intersection(const General_polygon_2<ArrTraits>& pgn1,
                                   const General_polygon_2<ArrTraits>& pgn2,
                                   OutputIterator out)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_intersection(pgn1, pgn2, out, traits);
}

// General_polygon_2, General_polygon_with_holes_2 =============================
// With traits
template <typename ArrTraits, typename OutputIterator, typename Traits>
inline OutputIterator intersection(const General_polygon_2<ArrTraits>& pgn1,
                                   const General_polygon_with_holes_2
                                     <General_polygon_2<ArrTraits> >& pgn2,
                                   OutputIterator out, Traits& traits)
{ return s_intersection(pgn1, pgn2, out, traits); }

// Without Traits
template <typename ArrTraits, typename OutputIterator>
inline OutputIterator intersection(const General_polygon_2<ArrTraits>& pgn1,
                                   const General_polygon_with_holes_2
                                     <General_polygon_2<ArrTraits> >& pgn2,
                                   OutputIterator out)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_intersection(pgn1, pgn2, out, traits);
}

// General_polygon_with_holes_2, General_polygon_2 =============================
// With traits
template <typename ArrTraits, typename OutputIterator, typename Traits>
inline OutputIterator intersection(const General_polygon_with_holes_2
                                     <General_polygon_2<ArrTraits> >& pgn1,
                                   const General_polygon_2<ArrTraits>& pgn2,
                                   OutputIterator out, Traits& traits)
{ return s_intersection(pgn1, pgn2, out, traits); }

// Without Traits
template <typename ArrTraits, typename OutputIterator>
inline OutputIterator intersection(const General_polygon_with_holes_2
                                     <General_polygon_2<ArrTraits> >& pgn1,
                                   const General_polygon_2<ArrTraits>& pgn2,
                                   OutputIterator out)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typedef General_polygon_with_holes_2<Polygon>         Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_intersection(pgn1, pgn2, out, traits);
}

// General_polygon_with_holes_2, General_polygon_with_holes_2 =================
// With traits
template <typename Polygon_, typename OutputIterator, typename Traits>
inline OutputIterator
intersection(const General_polygon_with_holes_2<Polygon_>& pgn1,
             const General_polygon_with_holes_2<Polygon_>& pgn2,
             OutputIterator out, Traits& traits)
{ return s_intersection(pgn1, pgn2, out, traits); }

// Without Traits
template <typename Polygon_, typename OutputIterator>
inline OutputIterator
intersection(const General_polygon_with_holes_2<Polygon_>& pgn1,
             const General_polygon_with_holes_2<Polygon_>& pgn2,
             OutputIterator out)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_with_holes_2<Polygon_>        Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_intersection(pgn1, pgn2, out, traits);
}

//@}

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

/// \name difference() functions.
//@{

// Polygon_2, Polygon_2 ========================================================
// With Traits
template <typename Kernel, typename Container, typename OutputIterator,
          typename Traits>
inline OutputIterator difference(const Polygon_2<Kernel, Container>& pgn1,
                                 const Polygon_2<Kernel, Container>& pgn2,
                                 OutputIterator oi, Traits& traits)
{ return _difference(pgn1, pgn2, oi, traits); }

// With Tag_true
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator difference(const Polygon_2<Kernel, Container>& pgn1,
                                 const Polygon_2<Kernel, Container>& pgn2,
                                 OutputIterator oi, Tag_true = Tag_true())
{ return _difference<Kernel, Container>(pgn1, pgn2, oi); }

// With Tag_false
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator difference(const Polygon_2<Kernel, Container>& pgn1,
                                 const Polygon_2<Kernel, Container>& pgn2,
                                 OutputIterator oi, Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return _difference(pgn1, pgn2, oi, traits);
}

// Polygon_2, Polygon_with_holes_2 =============================================
// With Traits
template <typename Kernel, typename Container, typename OutputIterator,
          typename Traits>
inline OutputIterator
difference(const Polygon_2<Kernel, Container>& pgn1,
           const Polygon_with_holes_2<Kernel, Container>& pgn2,
           OutputIterator oi, Traits& traits)
{ return _difference(pgn1, pgn2, oi, traits); }

// With Tag_true
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
difference(const Polygon_2<Kernel, Container>& pgn1,
           const Polygon_with_holes_2<Kernel, Container>& pgn2,
           OutputIterator oi, Tag_true = Tag_true())
{ return _difference<Kernel, Container>(pgn1, pgn2, oi); }

// With Tag_false
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
difference(const Polygon_2<Kernel, Container>& pgn1,
           const Polygon_with_holes_2<Kernel, Container>& pgn2,
           OutputIterator oi, Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return _difference(pgn1, pgn2, oi, traits);
}

// Polygon_with_holes_2, Polygon_2 =============================================
// With Traits
template <typename Kernel, typename Container, typename OutputIterator,
          typename Traits>
inline OutputIterator
difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
           const Polygon_2<Kernel, Container>& pgn2,
           OutputIterator oi, Traits& traits)
{ return _difference(pgn1, pgn2, oi, traits); }

// With Tag_true
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
           const Polygon_2<Kernel, Container>& pgn2,
           OutputIterator oi, Tag_true = Tag_true())
{ return _difference<Kernel, Container>(pgn1, pgn2, oi); }

// With Tag_false
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
           const Polygon_2<Kernel, Container>& pgn2,
           OutputIterator oi, Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_with_holes_2<Kernel, Container>       Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return _difference(pgn1, pgn2, oi, traits);
}

// Polygon_with_holes_2, Polygon_with_holes_2 ==================================
// With Traits
template <typename Kernel, typename Container, typename OutputIterator,
          typename Traits>
inline OutputIterator
difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
           const Polygon_with_holes_2<Kernel, Container>& pgn2,
           OutputIterator oi, Traits& traits)
{ return _difference(pgn1, pgn2, oi, traits); }

// With Tag_true
template <typename Kernel, typename Container, typename OutputIterator>
 OutputIterator
difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
           const Polygon_with_holes_2<Kernel, Container>& pgn2,
           OutputIterator oi, Tag_true = Tag_true())
{ return _difference<Kernel, Container>(pgn1, pgn2, oi); }

// With Tag_false
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
           const Polygon_with_holes_2<Kernel, Container>& pgn2,
           OutputIterator oi, Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_with_holes_2<Kernel, Container>       Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return _difference(pgn1, pgn2, oi, traits);
}

// General_polygon_2, General_polygon_2 ========================================
// With Traits
template <typename ArrTraits, typename OutputIterator, typename Traits>
inline OutputIterator difference(const General_polygon_2<ArrTraits>& pgn1,
                                 const General_polygon_2<ArrTraits>& pgn2,
                                 OutputIterator oi,
                                 Traits& traits)
{ return _difference(pgn1, pgn2, oi, traits); }

// Without Traits
template <typename ArrTraits, typename OutputIterator>
inline OutputIterator difference(const General_polygon_2<ArrTraits>& pgn1,
                                 const General_polygon_2<ArrTraits>& pgn2,
                                 OutputIterator oi)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return difference(pgn1, pgn2, oi, traits);
}

// General_polygon_2, General_polygon_with_holes_2 =============================
// With Traits
template <typename ArrTraits, typename OutputIterator, typename Traits>
inline OutputIterator difference(const General_polygon_2<ArrTraits>& pgn1,
                                 const General_polygon_with_holes_2
                                   <General_polygon_2<ArrTraits> >& pgn2,
                                 OutputIterator oi, Traits& traits)
{ return _difference(pgn1, pgn2, oi, traits); }

// Without Traits
template <typename ArrTraits, typename OutputIterator>
inline OutputIterator difference(const General_polygon_2<ArrTraits>& pgn1,
                                 const General_polygon_with_holes_2
                                   <General_polygon_2<ArrTraits> >& pgn2,
                                 OutputIterator oi)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return difference(pgn1, pgn2, oi, traits);
}

// General_polygon_with_holes_2, General_polygon_2 =============================
// With Traits
template <typename ArrTraits, typename OutputIterator, typename Traits>
inline OutputIterator difference(const General_polygon_with_holes_2
                                   <General_polygon_2<ArrTraits> >& pgn1,
                                 const General_polygon_2<ArrTraits>& pgn2,
                                 OutputIterator oi, Traits& traits)
{ return _difference(pgn1, pgn2, oi, traits); }

// Without Traits
template <typename ArrTraits, typename OutputIterator>
inline OutputIterator difference(const General_polygon_with_holes_2
                                   <General_polygon_2<ArrTraits> >& pgn1,
                                 const General_polygon_2<ArrTraits>& pgn2,
                                 OutputIterator oi)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typedef General_polygon_with_holes_2<Polygon>         Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return difference(pgn1, pgn2, oi, traits);
}

// General_polygon_with_holes_2, General_polygon_with_holes_2 ==================
// With Traits
template <typename Polygon_, typename OutputIterator, typename Traits>
inline OutputIterator
difference(const General_polygon_with_holes_2<Polygon_>& pgn1,
           const General_polygon_with_holes_2<Polygon_>& pgn2,
           OutputIterator oi, Traits& traits)
{ return _difference(pgn1, pgn2, oi, traits); }

// Without Traits
template <typename Polygon_, typename OutputIterator>
inline OutputIterator
difference(const General_polygon_with_holes_2<Polygon_>& pgn1,
           const General_polygon_with_holes_2<Polygon_>& pgn2,
           OutputIterator oi)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_with_holes_2<Polygon_>        Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return difference(pgn1, pgn2, oi, traits);
}

//@}

/// \name symmetric_difference() functions.
//@{

// Polygon_2, Polygon_2 ========================================================
// With Traits
template <typename Kernel, typename Container, typename OutputIterator,
          typename Traits>
inline OutputIterator
symmetric_difference(const Polygon_2<Kernel, Container>& pgn1,
                     const Polygon_2<Kernel, Container>& pgn2,
                     OutputIterator oi, Traits& traits)
{ return s_symmetric_difference(pgn1, pgn2, oi, traits); }

// With Tag_true
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
symmetric_difference(const Polygon_2<Kernel, Container>& pgn1,
                     const Polygon_2<Kernel, Container>& pgn2,
                     OutputIterator oi, Tag_true = Tag_true())
{ return s_symmetric_difference<Kernel, Container>(pgn1, pgn2, oi); }

// With Tag_false
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
symmetric_difference(const Polygon_2<Kernel, Container>& pgn1,
                     const Polygon_2<Kernel, Container>& pgn2,
                     OutputIterator oi, Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_symmetric_difference(pgn1, pgn2, oi, traits);
}

// Polygon_2, Polygon_with_holes_2 =============================================
// With Traits
template <typename Kernel, typename Container, typename OutputIterator,
          typename Traits>
inline OutputIterator
symmetric_difference(const Polygon_2<Kernel, Container>& pgn1,
                     const Polygon_with_holes_2<Kernel, Container>& pgn2,
                     OutputIterator oi, Traits& traits)
{ return s_symmetric_difference(pgn1, pgn2, oi, traits); }

// With Tag_true
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
symmetric_difference(const Polygon_2<Kernel, Container>& pgn1,
                     const Polygon_with_holes_2<Kernel, Container>& pgn2,
                     OutputIterator oi, Tag_true = Tag_true())
{ return s_symmetric_difference<Kernel, Container>(pgn1, pgn2, oi); }

// With Tag_false
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
symmetric_difference(const Polygon_2<Kernel, Container>& pgn1,
                     const Polygon_with_holes_2<Kernel, Container>& pgn2,
                     OutputIterator oi, Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_symmetric_difference(pgn1, pgn2, oi, traits);
}

// Polygon_with_holes_2, Polygon_2 =============================================
// With Traits
template <typename Kernel, typename Container, typename OutputIterator,
          typename Traits>
inline OutputIterator
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                     const Polygon_2<Kernel, Container>& pgn2,
                     OutputIterator oi, Traits& traits)
{ return s_symmetric_difference(pgn1, pgn2, oi, traits); }

// With Tag_true
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                     const Polygon_2<Kernel, Container>& pgn2,
                     OutputIterator oi, Tag_true = Tag_true())
{ return s_symmetric_difference<Kernel, Container>(pgn1, pgn2, oi); }

// With Tag_false
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                     const Polygon_2<Kernel, Container>& pgn2,
                     OutputIterator oi, Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_with_holes_2<Kernel, Container>       Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_symmetric_difference(pgn1, pgn2, oi, traits);
}

// Polygon_with_holes_2, Polygon_with_holes_2 ==================================
// With Traits
template <typename Kernel, typename Container, typename OutputIterator,
          typename Traits>
inline OutputIterator
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                     const Polygon_with_holes_2<Kernel, Container>& pgn2,
                     OutputIterator oi, Traits& traits)
{ return s_symmetric_difference(pgn1, pgn2, oi, traits); }

// With Tag_true
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                     const Polygon_with_holes_2<Kernel, Container>& pgn2,
                     OutputIterator oi, Tag_true = Tag_true())
{ return s_symmetric_difference<Kernel, Container>(pgn1, pgn2, oi); }

// With Tag_false
template <typename Kernel, typename Container, typename OutputIterator>
inline OutputIterator
symmetric_difference(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                     const Polygon_with_holes_2<Kernel, Container>& pgn2,
                     OutputIterator oi, Tag_false)
{
  // Use the first polygon to determine the (default) traits
  typedef Polygon_with_holes_2<Kernel, Container>       Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_symmetric_difference(pgn1, pgn2, oi, traits);
}

// General_polygon_2, General_polygon_2 ========================================
// With Traits
template <typename ArrTraits, typename OutputIterator, typename Traits>
inline OutputIterator
symmetric_difference(const General_polygon_2<ArrTraits>& pgn1,
                     const General_polygon_2<ArrTraits>& pgn2,
                     OutputIterator oi, Traits& traits)
{ return s_symmetric_difference(pgn1, pgn2, oi, traits); }

template <typename ArrTraits, typename OutputIterator>
inline OutputIterator
symmetric_difference(const General_polygon_2<ArrTraits>& pgn1,
                     const General_polygon_2<ArrTraits>& pgn2,
                     OutputIterator oi)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_symmetric_difference(pgn1, pgn2, oi, traits);
}

// General_polygon_2, General_polygon_with_holes_2 =============================
// With Traits
template <typename ArrTraits, typename OutputIterator, typename Traits>
inline OutputIterator
symmetric_difference(const General_polygon_2<ArrTraits>& pgn1,
                     const General_polygon_with_holes_2
                       <General_polygon_2<ArrTraits> >& pgn2,
                     OutputIterator oi, Traits& traits)
{ return s_symmetric_difference(pgn1, pgn2, oi, traits); }

// Without Traits
template <typename ArrTraits, typename OutputIterator>
inline OutputIterator
symmetric_difference(const General_polygon_2<ArrTraits>& pgn1,
                     const General_polygon_with_holes_2
                           <General_polygon_2<ArrTraits> >& pgn2,
                     OutputIterator oi)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return s_symmetric_difference(pgn1, pgn2, oi, traits);
}

// General_polygon_with_holes_2, General_polygon_2 =============================
// With Traits
template <typename ArrTraits, typename OutputIterator, typename Traits>
inline OutputIterator
symmetric_difference(const General_polygon_with_holes_2
                       <General_polygon_2<ArrTraits> >& pgn1,
                     const General_polygon_2<ArrTraits>& pgn2,
                     OutputIterator oi, Traits& traits)
{ return s_symmetric_difference(pgn1, pgn2, oi, traits); }

// Without Traits
template <typename ArrTraits, typename OutputIterator>
inline OutputIterator
symmetric_difference(const General_polygon_with_holes_2
                       <General_polygon_2<ArrTraits> >& pgn1,
                     const General_polygon_2<ArrTraits>& pgn2,
                     OutputIterator oi)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typedef General_polygon_with_holes_2<Polygon>         Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_symmetric_difference(pgn1, pgn2, oi, traits);
}

// General_polygon_with_holes_2, General_polygon_with_holes_2 ==================
// With Traits
template <typename Polygon_, typename OutputIterator, typename Traits>
inline OutputIterator
symmetric_difference(const General_polygon_with_holes_2<Polygon_>& pgn1,
                     const General_polygon_with_holes_2<Polygon_>& pgn2,
                     OutputIterator oi, Traits& traits)
{ return s_symmetric_difference(pgn1, pgn2, oi, traits); }

// Without Traits
template <typename Polygon_, typename OutputIterator>
inline OutputIterator
symmetric_difference(const General_polygon_with_holes_2<Polygon_>& pgn1,
                     const General_polygon_with_holes_2<Polygon_>& pgn2,
                     OutputIterator oi)
{
  typedef General_polygon_with_holes_2<Polygon_>        Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return s_symmetric_difference(pgn1, pgn2, oi, traits);
}

//@}

/// \name complement() functions.
//@{

// Polygon_2 ===================================================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
void complement(const Polygon_2<Kernel, Container>& pgn,
                Polygon_with_holes_2<Kernel, Container>& res,
                Traits& traits)
{ _complement(pgn, res, traits); }

// With Tag_true
template <typename Kernel, typename Container>
void complement(const Polygon_2<Kernel, Container>& pgn,
                Polygon_with_holes_2<Kernel, Container>& res,
                Tag_true = Tag_true())
{ _complement(pgn, res); }

// With Tag_false
template <typename Kernel, typename Container>
void complement(const Polygon_2<Kernel, Container>& pgn,
                Polygon_with_holes_2<Kernel, Container>& res,
                Tag_false)
{
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  complement(pgn, res, traits);
}

// Polygon_with_holes_2 ========================================================
// With Traits
template <typename Kernel, typename Container, typename OutputIterator,
          typename Traits>
OutputIterator complement(const Polygon_with_holes_2<Kernel, Container>& pgn,
                          OutputIterator oi, Traits& traits)
{ return _complement(pgn, oi, traits); }

// With Tag_true
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator complement(const Polygon_with_holes_2<Kernel, Container>& pgn,
                          OutputIterator oi, Tag_true = Tag_true())
{ return _complement(pgn, oi); }

// With Tag_false
template <typename Kernel, typename Container, typename OutputIterator>
OutputIterator complement(const Polygon_with_holes_2<Kernel, Container>& pgn,
                          OutputIterator oi, Tag_false)
{
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return complement(pgn, oi, traits);
}

// General_polygon_2 ===========================================================
// With Traits
template <typename ArrTraits, typename Traits>
void complement(const General_polygon_2<ArrTraits>& pgn,
                General_polygon_with_holes_2
                  <General_polygon_2<ArrTraits> >& res,
                Traits& traits)
{ _complement(pgn, res, traits); }

// Without Traits
template <typename ArrTraits>
void complement(const General_polygon_2<ArrTraits>& pgn,
                General_polygon_with_holes_2
                  <General_polygon_2<ArrTraits> >& res)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  complement(pgn, res, traits);
}

// General_polygon_with_holes_2 ================================================
// With Traits
template <typename Polygon_, typename OutputIterator, typename Traits>
OutputIterator complement(const General_polygon_with_holes_2<Polygon_>& pgn,
                          OutputIterator oi, Traits& traits)
{ return _complement(pgn, oi, traits); }

// Without Traits
template <typename Polygon_, typename OutputIterator>
OutputIterator complement(General_polygon_with_holes_2<Polygon_>& pgn,
                          OutputIterator oi)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_with_holes_2<Polygon_>        Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return complement(pgn, oi, traits);
}

//@}

/// \name Aggregated join() functions.
//@{

template <typename InputIterator>
struct map_iterator_to_traits {
  typedef typename std::iterator_traits<InputIterator>::value_type InputPolygon;
  typedef typename Gps_default_traits<InputPolygon>::Traits    Traits;
};

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
  typename map_iterator_to_traits<InputIterator>::Traits traits;
  return r_join(begin, end, oi, traits, k);
}

// General polygons or polygons with holes
template <typename InputIterator, typename OutputIterator>
inline OutputIterator join(InputIterator begin, InputIterator end,
                           OutputIterator oi, unsigned int k=5,
                           Disable_if_Polygon_2_iterator<InputIterator>* = 0)
{
  typename map_iterator_to_traits<InputIterator>::Traits traits;
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
  typename map_iterator_to_traits<InputIterator1>::Traits traits;
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
  typename map_iterator_to_traits<InputIterator1>::Traits traits;
  return r_join(begin1, end1, begin2, end2, oi, traits, k);
}

//@}

/// \name Aggregated intersection() functions.
//@{

// With Traits
template <typename InputIterator, typename OutputIterator, typename Traits>
inline OutputIterator intersection(InputIterator begin, InputIterator end,
                                   OutputIterator oi, Traits& traits,
                                   unsigned int k=5)
{ return r_intersection(begin, end, oi, traits, k); }

// Without Traits
// Tag_true => convert to polylines
template <typename InputIterator, typename OutputIterator>
inline OutputIterator
intersection(InputIterator begin, InputIterator end,
             OutputIterator oi, Tag_true = Tag_true(),
             unsigned int k=5,
             Enable_if_Polygon_2_iterator<InputIterator>* = 0)
{ return r_intersection(begin, end, oi, k); }

// Tag_false => do not convert to polylines
template <typename InputIterator, typename OutputIterator>
inline OutputIterator
intersection(InputIterator begin, InputIterator end,
             OutputIterator oi, Tag_false, unsigned int k=5,
             Enable_if_Polygon_2_iterator<InputIterator>* = 0)
{
  typename map_iterator_to_traits<InputIterator>::Traits traits;
  return r_intersection(begin, end, oi, traits, k);
}

// General polygons or polygons with holes
template <typename InputIterator, typename OutputIterator>
inline OutputIterator
intersection(InputIterator begin, InputIterator end,
             OutputIterator oi, unsigned int k=5,
             // workaround to avoid ambiguous calls with kernel functions
             typename boost::enable_if
               <typename CGAL::is_iterator<InputIterator>>::type* = 0,
             Disable_if_Polygon_2_iterator<InputIterator>* = 0)
{
  typename map_iterator_to_traits<InputIterator>::Traits traits;
  return r_intersection(begin, end, oi, traits, k);
}

// Inersect two ranges of simple polygons and polygons with holes.
// With Traits
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename Traits>
inline OutputIterator intersection(InputIterator1 begin1, InputIterator1 end1,
                                   InputIterator2 begin2, InputIterator2 end2,
                                   OutputIterator oi, Traits& traits,
                                   unsigned int k=5)
{ return r_intersection(begin1, end1, begin2, end2, oi, traits, k); }

// Without Traits
// Tag_true => convert to polylines
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline OutputIterator
intersection(InputIterator1 begin1, InputIterator1 end1,
             InputIterator2 begin2, InputIterator2 end2,
             OutputIterator oi, Tag_true = Tag_true(), unsigned int k=5,
             Enable_if_Polygon_2_iterator<InputIterator1>* = 0)
{ return r_intersection(begin1, end1, begin2, end2, oi, k); }

// Tag_false => do not convert to polylines
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline OutputIterator
intersection(InputIterator1 begin1, InputIterator1 end1,
             InputIterator2 begin2, InputIterator2 end2,
             OutputIterator oi, Tag_false, unsigned int k=5,
             Enable_if_Polygon_2_iterator<InputIterator1>* = 0)
{
  typename map_iterator_to_traits<InputIterator1>::Traits traits;
  return r_intersection(begin1, end1, begin2, end2, oi, traits, k);
}

// General polygons or polygons with holes
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline OutputIterator
intersection(InputIterator1 begin1, InputIterator1 end1,
             InputIterator2 begin2, InputIterator2 end2,
             OutputIterator oi, unsigned int k=5,
             Disable_if_Polygon_2_iterator<InputIterator1>* = 0)
{
  typename map_iterator_to_traits<InputIterator1>::Traits traits;
  return r_intersection(begin1, end1, begin2, end2, oi, traits, k);
}

//@}

/// \name Aggregated symmetric_difference() functions.
//@{

// With Traits
template <typename InputIterator, typename OutputIterator, typename Traits>
inline
OutputIterator symmetric_difference(InputIterator begin, InputIterator end,
                                    OutputIterator oi, Traits& traits,
                                    unsigned int k=5)
{ return r_symmetric_difference(begin, end, oi, traits, k); }

// Without Traits
// Tag_true => convert to polylines
template <typename InputIterator, typename OutputIterator>
inline OutputIterator
symmetric_difference(InputIterator begin, InputIterator end,
                     OutputIterator oi, Tag_true = Tag_true(), unsigned int k=5,
                     Enable_if_Polygon_2_iterator<InputIterator>* = 0)
{ return r_symmetric_difference(begin, end, oi, k); }

// Tag_false => do not convert to polylines
template <typename InputIterator, typename OutputIterator>
inline OutputIterator
symmetric_difference(InputIterator begin, InputIterator end,
                     OutputIterator oi, Tag_false, unsigned int k=5,
                     Enable_if_Polygon_2_iterator<InputIterator>* = 0)
{
  typename map_iterator_to_traits<InputIterator>::Traits traits;
  return r_symmetric_difference(begin, end, oi, traits, k);
}

// General polygons or polygons with holes
template <typename InputIterator, typename OutputIterator>
inline OutputIterator
symmetric_difference(InputIterator begin, InputIterator end,
                     OutputIterator oi, unsigned int k=5,
                     Disable_if_Polygon_2_iterator<InputIterator>* = 0)
{
  typename map_iterator_to_traits<InputIterator>::Traits traits;
  return r_symmetric_difference(begin, end, oi, traits, k);
}

// Xor two ranges of simple polygons and polygons with holes.
// With Traits
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, typename Traits>
inline
OutputIterator symmetric_difference(InputIterator1 begin1, InputIterator1 end1,
                                    InputIterator2 begin2, InputIterator2 end2,
                                    OutputIterator oi, Traits& traits,
                                    unsigned int k=5)
{ return r_symmetric_difference(begin1, end1, begin2, end2, oi, traits, k); }

// Without Traits
// Tag_true => convert to polylines
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline OutputIterator
symmetric_difference(InputIterator1 begin1, InputIterator1 end1,
                     InputIterator2 begin2, InputIterator2 end2,
                     OutputIterator oi, Tag_true = Tag_true(), unsigned int k=5,
                     Enable_if_Polygon_2_iterator<InputIterator1>* = 0)
{ return r_symmetric_difference(begin1, end1, begin2, end2, oi, k); }

// Tag_false => do not convert to polylines
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline OutputIterator
symmetric_difference(InputIterator1 begin1, InputIterator1 end1,
                     InputIterator2 begin2, InputIterator2 end2,
                     OutputIterator oi, Tag_false, unsigned int k=5,
                     Enable_if_Polygon_2_iterator<InputIterator1>* = 0)
{
  typename map_iterator_to_traits<InputIterator1>::Traits traits;
  return r_symmetric_difference(begin1, end1, begin2, end2, oi, traits, k);
}

// General polygons or polygons with holes
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline OutputIterator
symmetric_difference(InputIterator1 begin1, InputIterator1 end1,
                     InputIterator2 begin2, InputIterator2 end2,
                     OutputIterator oi, unsigned int k=5,
                     Disable_if_Polygon_2_iterator<InputIterator1>* = 0)
{
  typename map_iterator_to_traits<InputIterator1>::Traits traits;
  return r_symmetric_difference(begin1, end1, begin2, end2, oi, traits, k);
}

//@}

/// \name Aggregated do_intersect() functions.
//@{

// With Traits
template <typename InputIterator, typename Traits>
inline bool do_intersect(InputIterator begin, InputIterator end, Traits& traits,
                         unsigned int k=5)
{ return r_do_intersect(begin, end, traits, k); }

// Without Traits
// Tag_true => convert to polylines
template <typename InputIterator>
inline bool do_intersect(InputIterator begin, InputIterator end,
                         Tag_true = Tag_true(), unsigned int k=5,
                         Enable_if_Polygon_2_iterator<InputIterator>* = 0)
{ return r_do_intersect(begin, end, k); }

// Tag_false => do not convert to polylines
template <typename InputIterator>
inline bool do_intersect(InputIterator begin, InputIterator end,
                         Tag_false, unsigned int k=5,
                         Enable_if_Polygon_2_iterator<InputIterator>* = 0)
{
  typename map_iterator_to_traits<InputIterator>::Traits traits;
  return r_do_intersect(begin, end, traits, k);
}

// General polygons or polygons with holes
template <typename InputIterator>
inline bool do_intersect(InputIterator begin, InputIterator end,
                         unsigned int k=5,
                         Disable_if_Polygon_2_iterator<InputIterator>* = 0)
{
  typename map_iterator_to_traits<InputIterator>::Traits traits;
  return do_intersect(begin, end, traits, k);
}

// With Traits
template <typename InputIterator1, typename InputIterator2, typename Traits>
inline bool do_intersect(InputIterator1 begin1, InputIterator1 end1,
                         InputIterator2 begin2, InputIterator2 end2,
                         Traits& traits, unsigned int k=5)
{ return r_do_intersect(begin1, end1, begin2, end2, traits, k); }

// Without Traits
// Tag_true => convert to polylines
template <typename InputIterator1, typename InputIterator2>
inline bool do_intersect (InputIterator1 begin1, InputIterator1 end1,
                          InputIterator2 begin2, InputIterator2 end2,
                          Tag_true = Tag_true(), unsigned int k=5,
                          Enable_if_Polygon_2_iterator<InputIterator1>* = 0)
{ return r_do_intersect(begin1, end1, begin2, end2, k); }

// Tag_false => do not convert to polylines
template <typename InputIterator1, typename InputIterator2>
inline bool do_intersect (InputIterator1 begin1, InputIterator1 end1,
                          InputIterator2 begin2, InputIterator2 end2,
                          Tag_false, unsigned int k=5,
                          Enable_if_Polygon_2_iterator<InputIterator1>* = 0)
{ return r_do_intersect(begin1, end1, begin2, end2, k); }

// General polygons or polygons with holes
template <typename InputIterator1, typename InputIterator2>
inline bool do_intersect (InputIterator1 begin1, InputIterator1 end1,
                          InputIterator2 begin2, InputIterator2 end2,
                          unsigned int k=5,
                          Disable_if_Polygon_2_iterator<InputIterator1>* = 0)
{
  typename map_iterator_to_traits<InputIterator1>::Traits traits;
  return r_do_intersect(begin1, end1, begin2, end2, traits, k);
}

//@}

/// \name oriented_side() functions.
//@{

// Polygon_2, Polygon_2 ========================================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
inline Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
                                   const Polygon_2<Kernel, Container>& pgn2,
                                   Traits& traits)
{ return _oriented_side(pgn1, pgn2, traits); }

// With Tag_true
template <typename Kernel, typename Container>
inline Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
                                   const Polygon_2<Kernel, Container>& pgn2,
                                   Tag_true = Tag_true())
{ return _oriented_side(pgn1, pgn2); }

// With Tag_false
template <typename Kernel, typename Container>
inline Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
                                   const Polygon_2<Kernel, Container>& pgn2,
                                   Tag_false)
{
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return oriented_side(pgn1, pgn2, traits);
}

// Polygon_2, Polygon_with_holes_2 =============================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
inline
Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            Traits& traits)
{ return _oriented_side(pgn1, pgn2, traits); }

// With Tag_true
template <typename Kernel, typename Container>
inline
Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            Tag_true = Tag_true())
{ return _oriented_side(pgn1, pgn2); }

// With Tag_false
template <typename Kernel, typename Container>
inline
Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            Tag_false)
{
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return oriented_side(pgn1, pgn2, traits);
}

// Polygon_with_holes_2, Polygon_2 =============================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
inline
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            Traits& traits)
{ return _oriented_side(pgn1, pgn2, traits); }

// With Tag_true
template <typename Kernel, typename Container>
inline
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            Tag_true = Tag_true())
{ return _oriented_side(pgn1, pgn2); }

// With Tag_false
template <typename Kernel, typename Container>
inline
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            Tag_false)
{
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return oriented_side(pgn1, pgn2, traits);
}

// Polygon_with_holes_2, Polygon_with_holes_2 ==================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
inline
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            Traits& traits)
{ return _oriented_side(pgn1, pgn2, traits); }

// With Tag_true
template <typename Kernel, typename Container>
inline
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            Tag_true = Tag_true())
{ return _oriented_side(pgn1, pgn2); }

// With Tag_false
template <typename Kernel, typename Container>
inline
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            Tag_false)
{
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return oriented_side(pgn1, pgn2, traits);
}

// General_polygon_2, General_polygon_2 ========================================
// With Traits
template <typename ArrTraits, typename GpsTraits>
inline Oriented_side oriented_side(const General_polygon_2<ArrTraits>& pgn1,
                                   const General_polygon_2<ArrTraits>& pgn2,
                                   GpsTraits& traits)
{ return _oriented_side(pgn1, pgn2, traits); }

// Without Traits
template <typename ArrTraits>
inline Oriented_side oriented_side(const General_polygon_2<ArrTraits>& pgn1,
                                   const General_polygon_2<ArrTraits>& pgn2)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return oriented_side(pgn1, pgn2, traits);
}

// General_polygon_2, General_polygon_with_holes_2 =============================
// With Traits
template <typename ArrTraits, typename GpsTraits>
inline Oriented_side oriented_side(const General_polygon_2<ArrTraits>& pgn1,
                                   const General_polygon_with_holes_2
                                     <General_polygon_2<ArrTraits> >& pgn2,
                                   GpsTraits& traits)
{ return _oriented_side(pgn1, pgn2, traits); }

// Without Traits
template <typename ArrTraits>
inline Oriented_side oriented_side(const General_polygon_2<ArrTraits>& pgn1,
                                   const General_polygon_with_holes_2
                                     <General_polygon_2<ArrTraits> >& pgn2)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return oriented_side(pgn1, pgn2, traits);
}

// General_polygon_with_holes_2, General_polygon_2 =============================
// With Traits
template <typename ArrTraits, typename GpsTraits>
inline Oriented_side oriented_side(const General_polygon_with_holes_2
                                     <General_polygon_2<ArrTraits> >& pgn1,
                                   const General_polygon_2<ArrTraits>& pgn2,
                                   GpsTraits& traits)
{ return _oriented_side(pgn1, pgn2, traits); }

// Without Traits
template <typename ArrTraits>
inline Oriented_side oriented_side(const General_polygon_with_holes_2
                                     <General_polygon_2<ArrTraits> >& pgn1,
                                   const General_polygon_2<ArrTraits>& pgn2)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typedef General_polygon_with_holes_2<Polygon>         Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return oriented_side(pgn1, pgn2, traits);
}

// General_polygon_with_holes_2, General_polygon_with_holes_2 ==================
// With Traits
template <typename Polygon_, typename Traits>
inline
Oriented_side oriented_side(const General_polygon_with_holes_2<Polygon_>& pgn1,
                            const General_polygon_with_holes_2<Polygon_>& pgn2,
                            Traits& traits)
{ return _oriented_side(pgn1, pgn2, traits); }

// Without Traits
template <typename Polygon_>
inline
Oriented_side oriented_side(const General_polygon_with_holes_2<Polygon_>& pgn1,
                            const General_polygon_with_holes_2<Polygon_>& pgn2)
{
  typedef General_polygon_with_holes_2<Polygon_>        Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return oriented_side(pgn1, pgn2, traits);
}

// Point Query:
// Polygon_2 ===================================================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
inline Oriented_side oriented_side(const typename Kernel::Point_2& p,
                                   const Polygon_2<Kernel, Container>& pgn,
                                   Traits& traits)
{ return _oriented_side(p, pgn, traits); }

// With Tag_true
template <typename Kernel, typename Container>
inline Oriented_side oriented_side(const typename Kernel::Point_2& p,
                                   const Polygon_2<Kernel, Container>& pgn,
                                   Tag_true = Tag_true())
{ return _oriented_side(p, pgn); }

// With Tag_false
template <typename Kernel, typename Container>
inline Oriented_side oriented_side(const typename Kernel::Point_2& p,
                                   const Polygon_2<Kernel, Container>& pgn,
                                   Tag_false)
{
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return oriented_side(p, pgn, traits);
}

// Polygon_with_holes_2 ========================================================
// With Traits
template <typename Kernel, typename Container, typename Traits>
inline
Oriented_side oriented_side(const typename Kernel::Point_2& p,
                            const Polygon_with_holes_2<Kernel, Container>& pgn,
                            Traits& traits)
{ return _oriented_side(p, pgn, traits); }

// With Tag_true
template <typename Kernel, typename Container>
inline
Oriented_side oriented_side(const typename Kernel::Point_2& p,
                            const Polygon_with_holes_2<Kernel, Container>& pgn,
                            Tag_true = Tag_true())
{ return _oriented_side(p, pgn); }

// With Tag_false
template <typename Kernel, typename Container>
inline
Oriented_side oriented_side(const typename Kernel::Point_2& p,
                            const Polygon_with_holes_2<Kernel, Container>& pgn,
                            Tag_false)
{
  typedef Polygon_2<Kernel, Container>                          Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return oriented_side(p, pgn, traits);
}

// General_polygon_2 ===========================================================
// With Traits
template <typename ArrTraits, typename GpsTraits>
inline Oriented_side oriented_side(const typename ArrTraits::Point_2& p,
                                   const General_polygon_2<ArrTraits>& pgn,
                                   GpsTraits& traits)
{ return _oriented_side(p, pgn, traits); }

// Without Traits
template <typename ArrTraits>
inline Oriented_side oriented_side(const typename ArrTraits::Point_2& p,
                                   const General_polygon_2<ArrTraits>& pgn)
{
  // Use the first polygon to determine the (default) traits
  typedef General_polygon_2<ArrTraits>                  Polygon;
  typename Gps_default_traits<Polygon>::Traits traits;
  return oriented_side(p, pgn, traits);
}

// General_polygon_with_holes_2 ================================================
// With Traits
template <typename Polygon_, typename Traits>
inline
Oriented_side oriented_side(const typename Polygon_::Point_2& p,
                            const General_polygon_with_holes_2<Polygon_>& pgn,
                            Traits& traits)
{ return _oriented_side(p, pgn, traits); }

// Without Traits
template <typename Polygon_>
inline
Oriented_side oriented_side(const typename Polygon_::Point_2& p,
                            const General_polygon_with_holes_2<Polygon_>& pgn)
{
  typedef General_polygon_with_holes_2<Polygon_>        Polygon_with_holes;
  typename Gps_default_traits<Polygon_with_holes>::Traits traits;
  return oriented_side(p, pgn, traits);
}

//@}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
