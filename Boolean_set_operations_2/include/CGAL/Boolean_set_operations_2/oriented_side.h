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

#ifndef CGAL_BOOLEAN_SET_OPERATIONS_2_ORIENTED_SIDE_H
#define CGAL_BOOLEAN_SET_OPERATIONS_2_ORIENTED_SIDE_H

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
