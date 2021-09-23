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

#ifndef CGAL_BOOLEAN_SET_OPERATIONS_2_COMPLEMENT_H
#define CGAL_BOOLEAN_SET_OPERATIONS_2_COMPLEMENT_H

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

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
