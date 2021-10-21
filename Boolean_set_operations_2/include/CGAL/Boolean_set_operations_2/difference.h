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

#ifndef CGAL_BOOLEAN_SET_OPERATIONS_2_DIFFERENCE_H
#define CGAL_BOOLEAN_SET_OPERATIONS_2_DIFFERENCE_H

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

namespace CGAL
{

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

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
