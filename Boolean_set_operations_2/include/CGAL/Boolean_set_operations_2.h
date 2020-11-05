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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein        <wein@post.tau.ac.il>
//                 Efi Fogel       <efif@post.tau.ac.il>

#ifndef CGAL_BOOLEAN_SET_OPERATIONS_H
#define CGAL_BOOLEAN_SET_OPERATIONS_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/iterator.h>
#include <CGAL/Boolean_set_operations_2/Bso_internal_functions.h>
#include <CGAL/is_iterator.h>
#include <boost/utility/enable_if.hpp>

namespace CGAL {

/// \name do_intersect() functions.
//@{

template <class Kernel, class Container>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2)
{
  return (_do_intersect(pgn1, pgn2));
}

template <class Kernel, class Container, class Traits>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2,
                         Traits& tr)
{
  return (_do_intersect(pgn1, pgn2, tr));
}

template <class Kernel, class Container>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2)
{
  return (_do_intersect(pgn1, pgn2));
}

template <class Kernel, class Container, class Traits>
inline bool do_intersect(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2,
                         Traits& tr)
{
  return (_do_intersect(pgn1, pgn2, tr));
}

template <class Kernel, class Container>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2)
{
  return (_do_intersect(pgn1, pgn2));
}

template <class Kernel, class Container, class Traits>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_2<Kernel, Container>& pgn2,
                         Traits& tr)
{
  return (_do_intersect(pgn1, pgn2, tr));
}

template <class Kernel, class Container>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2)
{
  return (_do_intersect(pgn1, pgn2));
}

template <class Kernel, class Container, class Traits>
inline bool do_intersect(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2,
                         Traits& tr)
{
  return (_do_intersect(pgn1, pgn2, tr));
}

template <class Arr_traits>
inline bool do_intersect(const General_polygon_2<Arr_traits>& pgn1,
                         const General_polygon_2<Arr_traits>& pgn2)
{
  return (_do_intersect(pgn1, pgn2));
}

template <class Arr_traits, class Gps_traits>
inline bool do_intersect(const General_polygon_2<Arr_traits>& pgn1,
                         const General_polygon_2<Arr_traits>& pgn2,
                         Gps_traits& tr)
{
  return (_do_intersect(pgn1, pgn2, tr));
}

template <class Arr_traits>
inline bool do_intersect(const General_polygon_2<Arr_traits>& pgn1,
                         const General_polygon_with_holes_2
                               <General_polygon_2<Arr_traits> >& pgn2)
{
  return (_do_intersect(pgn1, pgn2));
}

template <class Arr_traits, class Gps_traits>
inline bool do_intersect(const General_polygon_2<Arr_traits>& pgn1,
                         const General_polygon_with_holes_2
                               <General_polygon_2<Arr_traits> >& pgn2,
                         Gps_traits& tr)
{
  return (_do_intersect(pgn1, pgn2, tr));
}



template <class Arr_traits>
inline bool do_intersect(const General_polygon_with_holes_2
                               <General_polygon_2<Arr_traits> >& pgn1,
                         const General_polygon_2<Arr_traits>& pgn2)
{
  return (_do_intersect(pgn1, pgn2));
}

template <class Arr_traits, class Gps_traits>
inline bool do_intersect(const General_polygon_with_holes_2
                               <General_polygon_2<Arr_traits> >& pgn1,
                         const General_polygon_2<Arr_traits>& pgn2,
                         Gps_traits& tr)
{
  return (_do_intersect(pgn1, pgn2, tr));
}

template <class Polygon_>
inline bool do_intersect(const General_polygon_with_holes_2<Polygon_>& pgn1,
                         const General_polygon_with_holes_2<Polygon_>& pgn2)
{
  return (_do_intersect(pgn1, pgn2));
}

template <class Polygon_, class Traits>
inline bool do_intersect(const General_polygon_with_holes_2<Polygon_>& pgn1,
                         const General_polygon_with_holes_2<Polygon_>& pgn2,
                         Traits& tr)
{
  return (_do_intersect(pgn1, pgn2, tr));
}

//@}
/// \name intersection() functions.
//@{

template <class Kernel, class Container, typename OutputIterator>
inline OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1,
                                   const Polygon_2<Kernel, Container>& pgn2,
                                   OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Kernel, class Container, typename OutputIterator, class Traits>
inline OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1,
                                   const Polygon_2<Kernel, Container>& pgn2,
                                   OutputIterator out,
                                   Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

template <class Kernel, class Container, typename OutputIterator>
inline OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1,
                         const Polygon_with_holes_2<Kernel, Container>& pgn2,
                           OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Kernel, class Container, typename OutputIterator, class Traits>
inline OutputIterator
intersection (const Polygon_2<Kernel, Container>& pgn1,
              const Polygon_with_holes_2<Kernel, Container>& pgn2,
              OutputIterator out, Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

template <class Kernel, class Container, typename OutputIterator>
inline OutputIterator
intersection (const Polygon_with_holes_2<Kernel, Container>& pgn1,
              const Polygon_2<Kernel, Container>& pgn2,
              OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Kernel, class Container, typename OutputIterator, class Traits>
inline OutputIterator
intersection (const Polygon_with_holes_2<Kernel, Container>& pgn1,
              const Polygon_2<Kernel, Container>& pgn2,
              OutputIterator out, Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

template <class Kernel, class Container, typename OutputIterator>
inline OutputIterator
intersection (const Polygon_with_holes_2<Kernel, Container>& pgn1,
              const Polygon_with_holes_2<Kernel, Container>& pgn2,
              OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Kernel, class Container, typename OutputIterator, class Traits>
inline OutputIterator
intersection (const Polygon_with_holes_2<Kernel, Container>& pgn1,
              const Polygon_with_holes_2<Kernel, Container>& pgn2,
              OutputIterator out,
              Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

template <class Arr_traits, typename OutputIterator>
inline OutputIterator intersection (const General_polygon_2<Arr_traits>& pgn1,
                                    const General_polygon_2<Arr_traits>& pgn2,
                                    OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Arr_traits, typename OutputIterator, class Traits>
inline OutputIterator intersection (const General_polygon_2<Arr_traits>& pgn1,
                                    const General_polygon_2<Arr_traits>& pgn2,
                                    OutputIterator out,
                                    Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

template <class Arr_traits, typename OutputIterator>
inline OutputIterator
intersection (const General_polygon_2<Arr_traits>& pgn1,
              const General_polygon_with_holes_2
                    <General_polygon_2<Arr_traits> >& pgn2,
              OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Arr_traits, typename OutputIterator, class Traits>
inline OutputIterator
intersection (const General_polygon_2<Arr_traits>& pgn1,
              const General_polygon_with_holes_2
                    <General_polygon_2<Arr_traits> >& pgn2,
              OutputIterator out, Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

template <class Arr_traits, typename OutputIterator>
inline OutputIterator
intersection (const General_polygon_with_holes_2
                    <General_polygon_2<Arr_traits> >& pgn1,
              const General_polygon_2<Arr_traits>& pgn2,
              OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Arr_traits, typename OutputIterator, class Traits>
inline OutputIterator
intersection (const General_polygon_with_holes_2
                    <General_polygon_2<Arr_traits> >& pgn1,
              const General_polygon_2<Arr_traits>& pgn2,
              OutputIterator out, Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

template <class Polygon_, typename OutputIterator>
inline OutputIterator
intersection (const General_polygon_with_holes_2<Polygon_>& pgn1,
              const General_polygon_with_holes_2<Polygon_>& pgn2,
              OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Polygon_, typename OutputIterator, class Traits>
inline OutputIterator
intersection (const General_polygon_with_holes_2<Polygon_>& pgn1,
              const General_polygon_with_holes_2<Polygon_>& pgn2,
              OutputIterator out, Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

//@}
/// \name _join() functions.
//@{

template <class Kernel, class Container>
inline bool join (const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2,
                  Polygon_with_holes_2<Kernel, Container>& res)
{
  return (_join(pgn1, pgn2, res));
}

template <class Kernel, class Container, class Traits>
inline bool join (const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2,
                  Polygon_with_holes_2<Kernel, Container>& res, Traits& tr)
{
  return (_join(pgn1, pgn2, res, tr));
}

template <class Kernel, class Container>
inline bool join (const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2,
                  Polygon_with_holes_2<Kernel, Container>& res)
{
  return (_join(pgn1, pgn2, res));
}

template <class Kernel, class Container, class Traits>
inline bool join (const Polygon_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2,
                  Polygon_with_holes_2<Kernel, Container>& res, Traits& tr)
{
  return (_join(pgn1, pgn2, res, tr));
}

template <class Kernel, class Container>
inline bool join (const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2,
                  Polygon_with_holes_2<Kernel, Container>& res )
{
  return (_join(pgn1, pgn2, res));
}

template <class Kernel, class Container, class Traits>
inline bool join (const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_2<Kernel, Container>& pgn2,
                  Polygon_with_holes_2<Kernel, Container>& res, Traits& tr)
{
  return (_join(pgn1, pgn2, res, tr));
}

template <class Kernel, class Container>
inline bool join (const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2,
                  Polygon_with_holes_2<Kernel, Container>& res )
{
  return (_join(pgn1, pgn2, res));
}

template <class Kernel, class Container, class Traits>
inline bool join (const Polygon_with_holes_2<Kernel, Container>& pgn1,
                  const Polygon_with_holes_2<Kernel, Container>& pgn2,
                  Polygon_with_holes_2<Kernel, Container>& res, Traits& tr)
{
  return (_join(pgn1, pgn2, res, tr));
}

template <class Arr_traits>
inline bool
join (const General_polygon_2<Arr_traits>& pgn1,
      const General_polygon_2<Arr_traits>& pgn2,
      General_polygon_with_holes_2<General_polygon_2<Arr_traits> >& res)
{
  return (_join(pgn1, pgn2, res));
}

template <class Arr_traits, class Traits>
inline bool
join (const General_polygon_2<Arr_traits>& pgn1,
      const General_polygon_2<Arr_traits>& pgn2,
      General_polygon_with_holes_2<General_polygon_2<Arr_traits> >& res,
      Traits& tr)
{
  return (_join(pgn1, pgn2, res, tr));
}

template <class Arr_traits>
inline bool
join (const General_polygon_2<Arr_traits>& pgn1,
      const General_polygon_with_holes_2<General_polygon_2<Arr_traits> >& pgn2,
      General_polygon_with_holes_2<General_polygon_2<Arr_traits> >& res)
{
  return (_join(pgn1, pgn2, res));
}

template <class Arr_traits, class Traits>
inline bool
join (const General_polygon_2<Arr_traits>& pgn1,
      const General_polygon_with_holes_2<General_polygon_2<Arr_traits> >& pgn2,
      General_polygon_with_holes_2<General_polygon_2<Arr_traits> >& res,
      Traits& tr)
{
  return (_join(pgn1, pgn2, res, tr));
}

template <class Arr_traits>
inline bool
join (const General_polygon_with_holes_2<General_polygon_2<Arr_traits> >& pgn1,
      const General_polygon_2<Arr_traits>& pgn2,
      General_polygon_with_holes_2<General_polygon_2<Arr_traits> >& res)
{
  return (_join(pgn1, pgn2, res));
}

template <class Arr_traits, class Traits>
inline bool
join (const General_polygon_with_holes_2<General_polygon_2<Arr_traits> >& pgn1,
      const General_polygon_2<Arr_traits>& pgn2,
      General_polygon_with_holes_2<General_polygon_2<Arr_traits> >& res,
      Traits& tr)
{
  return (_join(pgn1, pgn2, res, tr));
}

template <class Polygon_>
inline bool join (const General_polygon_with_holes_2<Polygon_>& pgn1,
                  const General_polygon_with_holes_2<Polygon_>& pgn2,
                  General_polygon_with_holes_2<Polygon_>& res)
{
  return (_join(pgn1, pgn2, res));
}

template <class Polygon_, class Traits>
inline bool join (const General_polygon_with_holes_2<Polygon_>& pgn1,
                  const General_polygon_with_holes_2<Polygon_>& pgn2,
                  General_polygon_with_holes_2<Polygon_>& res, Traits& tr)
{
  return (_join(pgn1, pgn2, res, tr));
}

//@}
/// \name difference() functions.
//@{

template <class Kernel, class Container, typename OutputIterator>
inline OutputIterator difference (const Polygon_2<Kernel, Container>& pgn1,
                                  const Polygon_2<Kernel, Container>& pgn2,
                                  OutputIterator oi)
{
  return(_difference(pgn1, pgn2, oi));
}

template <class Kernel, class Container, typename OutputIterator, class Traits>
inline OutputIterator difference (const Polygon_2<Kernel, Container>& pgn1,
                                  const Polygon_2<Kernel, Container>& pgn2,
                                  OutputIterator oi, Traits& tr)
{
  return(_difference(pgn1, pgn2, oi, tr));
}

template <class Kernel, class Container, typename OutputIterator>
inline OutputIterator
difference (const Polygon_2<Kernel, Container>& pgn1,
            const Polygon_with_holes_2<Kernel, Container>& pgn2,
            OutputIterator oi)
{
  return(_difference(pgn1, pgn2, oi));
}

template <class Kernel, class Container, typename OutputIterator, class Traits>
inline OutputIterator
difference (const Polygon_2<Kernel, Container>& pgn1,
            const Polygon_with_holes_2<Kernel, Container>& pgn2,
            OutputIterator oi, Traits& tr)
{
  return(_difference(pgn1, pgn2, oi, tr));
}

template <class Kernel, class Container, typename OutputIterator>
inline OutputIterator
difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
            const Polygon_2<Kernel, Container>& pgn2,
            OutputIterator oi)
{
  return (_difference(pgn1, pgn2, oi));
}

template <class Kernel, class Container, typename OutputIterator, class Traits>
inline OutputIterator
difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
            const Polygon_2<Kernel, Container>& pgn2,
            OutputIterator oi,
            Traits& tr)
{
  return (_difference(pgn1, pgn2, oi, tr));
}

template <class Kernel, class Container, typename OutputIterator>
inline OutputIterator
difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
            const Polygon_with_holes_2<Kernel, Container>& pgn2,
            OutputIterator oi)
{
  return (_difference(pgn1, pgn2, oi));
}

template <class Kernel, class Container, typename OutputIterator, class Traits>
inline OutputIterator
difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
            const Polygon_with_holes_2<Kernel, Container>& pgn2,
            OutputIterator oi, Traits& tr)
{
  return (_difference(pgn1, pgn2, oi, tr));
}

template <class Arr_traits, typename OutputIterator>
inline OutputIterator difference (const General_polygon_2<Arr_traits>& pgn1,
                                  const General_polygon_2<Arr_traits>& pgn2,
                                  OutputIterator oi)
{
  return (_difference(pgn1, pgn2, oi));
}

template <class Arr_traits, typename OutputIterator, class Traits>
inline OutputIterator difference (const General_polygon_2<Arr_traits>& pgn1,
                                  const General_polygon_2<Arr_traits>& pgn2,
                                  OutputIterator oi,
                                  Traits& tr)
{
  return (_difference(pgn1, pgn2, oi, tr));
}

template <class Arr_traits, typename OutputIterator>
inline OutputIterator difference (const General_polygon_2<Arr_traits>& pgn1,
                                  const General_polygon_with_holes_2
                                        <General_polygon_2<Arr_traits> >& pgn2,
                                  OutputIterator oi)
{
  return (_difference(pgn1, pgn2, oi));
}

template <class Arr_traits, typename OutputIterator, class Traits>
inline OutputIterator difference (const General_polygon_2<Arr_traits>& pgn1,
                                  const General_polygon_with_holes_2
                                        <General_polygon_2<Arr_traits> >& pgn2,
                                  OutputIterator oi, Traits& tr)
{
  return (_difference(pgn1, pgn2, oi, tr));
}

template <class Arr_traits, typename OutputIterator>
inline OutputIterator difference (const General_polygon_with_holes_2
                                        <General_polygon_2<Arr_traits> >& pgn1,
                                  const General_polygon_2<Arr_traits>& pgn2,
                                  OutputIterator oi)
{
  return (_difference(pgn1, pgn2, oi));
}

template <class Arr_traits, typename OutputIterator, class Traits>
inline OutputIterator difference (const General_polygon_with_holes_2
                                        <General_polygon_2<Arr_traits> >& pgn1,
                                  const General_polygon_2<Arr_traits>& pgn2,
                                  OutputIterator oi, Traits& tr)
{
  return (_difference(pgn1, pgn2, oi, tr));
}

template <class Polygon_, typename OutputIterator>
inline OutputIterator
difference (const General_polygon_with_holes_2<Polygon_>& pgn1,
            const General_polygon_with_holes_2<Polygon_>& pgn2,
            OutputIterator oi)
{
  return (_difference(pgn1, pgn2, oi));
}

template <class Polygon_, typename OutputIterator, class Traits>
inline OutputIterator
difference (const General_polygon_with_holes_2<Polygon_>& pgn1,
            const General_polygon_with_holes_2<Polygon_>& pgn2,
            OutputIterator oi, Traits& tr)
{
  return (_difference(pgn1, pgn2, oi, tr));
}

//@}
/// \name symmetric_difference() functions.
//@{

template <class Kernel, class Container, typename OutputIterator, class Traits>
inline OutputIterator
symmetric_difference (const Polygon_2<Kernel, Container>& pgn1,
                      const Polygon_2<Kernel, Container>& pgn2,
                      OutputIterator oi, Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Kernel, class Container, typename OutputIterator>
inline OutputIterator
symmetric_difference (const Polygon_2<Kernel, Container>& pgn1,
                      const Polygon_2<Kernel, Container>& pgn2,
                      OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

template <class Kernel, class Container, typename OutputIterator, class Traits>
inline OutputIterator
symmetric_difference (const Polygon_2<Kernel, Container>& pgn1,
                      const Polygon_with_holes_2<Kernel, Container>& pgn2,
                      OutputIterator oi, Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}


template <class Kernel, class Container, typename OutputIterator>
inline OutputIterator
symmetric_difference (const Polygon_2<Kernel, Container>& pgn1,
                      const Polygon_with_holes_2<Kernel, Container>& pgn2,
                      OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

template <class Kernel, class Container, typename OutputIterator, class Traits>
inline OutputIterator
symmetric_difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
                      const Polygon_2<Kernel, Container>& pgn2,
                      OutputIterator oi, Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Kernel, class Container, typename OutputIterator>
inline OutputIterator
symmetric_difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
                      const Polygon_2<Kernel, Container>& pgn2,
                      OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

template <class Kernel, class Container, typename OutputIterator, class Traits>
inline OutputIterator
symmetric_difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
                      const Polygon_with_holes_2<Kernel, Container>& pgn2,
                      OutputIterator oi, Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Kernel, class Container, typename OutputIterator>
inline OutputIterator
symmetric_difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
                      const Polygon_with_holes_2<Kernel, Container>& pgn2,
                      OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

template <class Arr_traits, typename OutputIterator, class Traits>
inline OutputIterator
symmetric_difference (const General_polygon_2<Arr_traits>& pgn1,
                      const General_polygon_2<Arr_traits>& pgn2,
                      OutputIterator oi,
                      Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Arr_traits, typename OutputIterator>
inline OutputIterator
symmetric_difference (const General_polygon_2<Arr_traits>& pgn1,
                      const General_polygon_2<Arr_traits>& pgn2,
                      OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

template <class Arr_traits, typename OutputIterator, class Traits>
inline OutputIterator
symmetric_difference (const General_polygon_2<Arr_traits>& pgn1,
                      const General_polygon_with_holes_2
                            <General_polygon_2<Arr_traits> >& pgn2,
                      OutputIterator oi, Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Arr_traits, typename OutputIterator>
inline OutputIterator
symmetric_difference(const General_polygon_2<Arr_traits>& pgn1,
                     const General_polygon_with_holes_2
                           <General_polygon_2<Arr_traits> >& pgn2,
                     OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

template <class Arr_traits, typename OutputIterator, class Traits>
inline OutputIterator
symmetric_difference (const General_polygon_with_holes_2
                        <General_polygon_2<Arr_traits> >& pgn1,
                      const General_polygon_2<Arr_traits>& pgn2,
                      OutputIterator oi, Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Arr_traits, typename OutputIterator>
inline OutputIterator
symmetric_difference (const General_polygon_with_holes_2
                            <General_polygon_2<Arr_traits> >& pgn1,
                      const General_polygon_2<Arr_traits>& pgn2,
                      OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

template <class Polygon_, typename OutputIterator, class Traits>
inline OutputIterator
symmetric_difference (const General_polygon_with_holes_2<Polygon_>& pgn1,
                      const General_polygon_with_holes_2<Polygon_>& pgn2,
                      OutputIterator oi,
                      Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Polygon_, typename OutputIterator>
inline OutputIterator
symmetric_difference (const General_polygon_with_holes_2<Polygon_>& pgn1,
                      const General_polygon_with_holes_2<Polygon_>& pgn2,
                      OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

//@}
/// \name complement() functions.
//@{

template <class Kernel, class Container>
void complement (const Polygon_2<Kernel, Container>& pgn,
                 Polygon_with_holes_2<Kernel, Container>& res)
{
  _complement(pgn, res);
}

template <class Kernel, class Container, class Traits>
void complement (const Polygon_2<Kernel, Container>& pgn,
                 Polygon_with_holes_2<Kernel, Container>& res,
                 Traits& tr)
{
  _complement(pgn, res, tr);
}

template <class Arr_traits>
void complement (const General_polygon_2<Arr_traits>& pgn,
                 General_polygon_with_holes_2
                   <General_polygon_2<Arr_traits> >& res)
{
  _complement(pgn, res);
}

template <class Arr_traits, class Traits>
void complement (const General_polygon_2<Arr_traits>& pgn,
                 General_polygon_with_holes_2
                   <General_polygon_2<Arr_traits> >& res,
                 Traits& tr)
{
  _complement(pgn, res, tr);
}

template <class Kernel, class Container, typename OutputIterator, class Traits>
OutputIterator complement (const Polygon_with_holes_2<Kernel, Container>& pgn,
                           OutputIterator oi, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn);
  gps.complement();
  return (gps.polygons_with_holes(oi));
}

template <class Kernel, class Container, typename OutputIterator>
OutputIterator complement (const Polygon_with_holes_2<Kernel, Container>& pgn,
                           OutputIterator oi)
{
  typename Gps_default_traits<Polygon_2<Kernel, Container> >::Traits  tr;
  return (complement(pgn, oi, tr));
}

template <class Arr_traits, typename OutputIterator, class Traits>
OutputIterator complement (const General_polygon_with_holes_2<Arr_traits>& pgn,
                           OutputIterator oi, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn);
  gps.complement();
  return (gps.polygons_with_holes(oi));
}

template <class Arr_traits, typename OutputIterator>
OutputIterator complement (General_polygon_with_holes_2
                             <General_polygon_2<Arr_traits> >& pgn,
                           OutputIterator oi)
{
  typename Gps_default_traits<General_polygon_2<Arr_traits> >::Traits  tr;
  return (complement(pgn, oi, tr));
}

//@}
/// \name Aggregated join() functions.
//@{

template <typename InputIterator>
struct map_iterator_to_traits
{
  typedef typename std::iterator_traits<InputIterator>::value_type InputPolygon;
  typedef typename Gps_default_traits<InputPolygon>::Traits    Traits;
};

template <typename InputIterator, typename OutputIterator, class Traits>
inline OutputIterator join(InputIterator begin, InputIterator end,
                           OutputIterator oi, Traits&, unsigned int k=5)
{
  if (begin == end)
    return (oi);

  General_polygon_set_2<Traits> gps(*begin);
  gps.join(++begin, end, k);
  return (gps.polygons_with_holes(oi));
}

template <typename InputIterator, typename OutputIterator>
inline OutputIterator join(InputIterator begin, InputIterator end,
                           OutputIterator oi, unsigned int k=5)
{
  typename map_iterator_to_traits<InputIterator>::Traits          tr;
  return join(begin, end, oi, tr, k);
}

// Join two ranges of simple polygons and polygons with holes.
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, class Traits>
inline OutputIterator join(InputIterator1 begin1, InputIterator1 end1,
                           InputIterator2 begin2, InputIterator2 end2,
                           OutputIterator oi, Traits& tr, unsigned int k=5)
{
  if (begin1 == end1)
    return (join(begin2, end2, oi, tr, k));

  General_polygon_set_2<Traits> gps(*begin1);
  gps.join(++begin1, end1, begin2, end2, k);
  return (gps.polygons_with_holes(oi));
}

template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline OutputIterator join(InputIterator1 begin1, InputIterator1 end1,
                           InputIterator2 begin2, InputIterator2 end2,
                           OutputIterator oi, unsigned int k=5)
{
  typename map_iterator_to_traits<InputIterator1>::Traits  tr;
  return join(begin1, end1, begin2, end2, oi, tr, k);
}

//@}
/// \name Aggregated intersection() functions.
//@{
template <typename InputIterator, typename OutputIterator, class Traits>
inline OutputIterator intersection (InputIterator begin, InputIterator end,
                                    OutputIterator oi, Traits&, unsigned int k=5)
{
  if (begin == end)
    return (oi);

  General_polygon_set_2<Traits> gps(*begin);
  gps.intersection(++begin, end, k);
  return (gps.polygons_with_holes(oi));
}

template <typename InputIterator, typename OutputIterator>
inline OutputIterator intersection (InputIterator begin, InputIterator end,
                                    OutputIterator oi, unsigned int k=5,
                                    typename boost::enable_if<
                                      typename CGAL::is_iterator<InputIterator>
                                    >::type* =0 // workaround to avoid ambiguous calls with kernel functions
)
{
  typename map_iterator_to_traits<InputIterator>::Traits          tr;
  return intersection(begin, end, oi, tr, k);
}

// Inersect two ranges of simple polygons and polygons with holes.
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, class Traits>
inline OutputIterator intersection(InputIterator1 begin1, InputIterator1 end1,
                                   InputIterator2 begin2, InputIterator2 end2,
                                   OutputIterator oi, Traits& tr,
                                   unsigned int k=5)
{
  if (begin1 == end1)
    return (intersection(begin2, end2, oi, tr, k));

  General_polygon_set_2<Traits> gps(*begin1);
  gps.intersection(++begin1, end1, begin2, end2, k);
  return (gps.polygons_with_holes(oi));

}

template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline OutputIterator intersection(InputIterator1 begin1, InputIterator1 end1,
                                   InputIterator2 begin2, InputIterator2 end2,
                                   OutputIterator oi, unsigned int k=5)
{
  typename map_iterator_to_traits<InputIterator1>::Traits  tr;
  return intersection(begin1, end1, begin2, end2, oi, tr, k);
}

//@}
/// \name Aggregated symmetric_difference() functions.
//@{
template <typename InputIterator, typename OutputIterator, class Traits>
inline
OutputIterator symmetric_difference(InputIterator begin, InputIterator end,
                                    OutputIterator oi, Traits& tr,
                                    unsigned int k=5)
{
  if (begin == end)
    return (oi);

  General_polygon_set_2<Traits> gps(tr);
  gps.insert(*begin);
  gps.symmetric_difference(++begin, end, k);
  return (gps.polygons_with_holes(oi));
}

template <typename InputIterator, typename OutputIterator>
inline
OutputIterator symmetric_difference (InputIterator begin, InputIterator end,
                                     OutputIterator oi, unsigned int k=5)
{
  typename map_iterator_to_traits<InputIterator>::Traits          tr;
  return symmetric_difference(begin, end, oi, tr, k);
}

// Xor two ranges of simple polygons and polygons with holes.
template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator, class Traits>
inline
OutputIterator symmetric_difference (InputIterator1 begin1, InputIterator1 end1,
                                     InputIterator2 begin2, InputIterator2 end2,
                                     OutputIterator oi, Traits& tr,
                                     unsigned int k=5)
{
  if (begin1 == end1)
    return (symmetric_difference(begin2, end2, oi, tr, k));

  General_polygon_set_2<Traits> gps(tr);
  gps.insert(*begin1);
  gps.symmetric_difference(++begin1, end1, begin2, end2, k);
  return (gps.polygons_with_holes(oi));

}

template <typename InputIterator1, typename InputIterator2,
          typename OutputIterator>
inline
OutputIterator symmetric_difference (InputIterator1 begin1, InputIterator1 end1,
                                     InputIterator2 begin2, InputIterator2 end2,
                                     OutputIterator oi, unsigned int k=5)
{
  typename map_iterator_to_traits<InputIterator1>::Traits  tr;
  return symmetric_difference(begin1, end1, begin2, end2, oi, tr, k);
}

//@}
/// \name Aggregated do_intersect() functions.
//@{

template <typename InputIterator, class Traits>
inline bool do_intersect(InputIterator begin, InputIterator end, Traits& tr,
                         unsigned int k=5)
{
  if (begin == end)
    return false;

  General_polygon_set_2<Traits> gps(tr);
  gps.insert(*begin);
  return gps.do_intersect(++begin, end, k);
}

template <typename InputIterator>
inline bool do_intersect(InputIterator begin, InputIterator end,
                         unsigned int k=5)
{
  typename map_iterator_to_traits<InputIterator>::Traits  tr;
  return do_intersect(begin, end, tr, k);
}

template <typename InputIterator1, typename InputIterator2, class Traits>
inline bool do_intersect (InputIterator1 begin1, InputIterator1 end1,
                          InputIterator2 begin2, InputIterator2 end2,
                          Traits& tr, unsigned int k=5)
{
  if (begin1 == end1)
    return (do_intersect(begin2, end2, tr, k));

  General_polygon_set_2<Traits> gps(tr);
  gps.insert(*begin1);
  return gps.do_intersect(++begin1, end1, begin2, end2, k);
}

template <typename InputIterator1, typename InputIterator2>
inline bool do_intersect (InputIterator1 begin1, InputIterator1 end1,
                          InputIterator2 begin2, InputIterator2 end2,
                          unsigned int k=5)
{
  typename map_iterator_to_traits<InputIterator1>::Traits  tr;
  return do_intersect(begin1, end1, begin2, end2, tr, k);
}

//@}

/// \name oriented_side() functions.
//@{

template <class Kernel, class Container>
inline Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
                                   const Polygon_2<Kernel, Container>& pgn2)
{
  return (_oriented_side(pgn1, pgn2));
}

template <class Kernel, class Container, class Traits>
inline Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
                                   const Polygon_2<Kernel, Container>& pgn2,
                                   Traits& tr)
{
  return (_oriented_side(pgn1, pgn2, tr));
}

template <class Kernel, class Container>
inline
Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2)
{
  return (_oriented_side(pgn1, pgn2));
}

template <class Kernel, class Container, class Traits>
inline
Oriented_side oriented_side(const Polygon_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            Traits& tr)
{
  return (_oriented_side(pgn1, pgn2, tr));
}

template <class Kernel, class Container>
inline
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2)
{
  return (_oriented_side(pgn1, pgn2));
}

template <class Kernel, class Container, class Traits>
inline
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_2<Kernel, Container>& pgn2,
                            Traits& tr)
{
  return (_oriented_side(pgn1, pgn2, tr));
}

template <class Kernel, class Container>
inline
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2)
{
  return (_oriented_side(pgn1, pgn2));
}

template <class Kernel, class Container, class Traits>
inline
Oriented_side oriented_side(const Polygon_with_holes_2<Kernel, Container>& pgn1,
                            const Polygon_with_holes_2<Kernel, Container>& pgn2,
                            Traits& tr)
{
  return (_oriented_side(pgn1, pgn2, tr));
}

template <class Arr_traits>
inline Oriented_side oriented_side(const General_polygon_2<Arr_traits>& pgn1,
                                   const General_polygon_2<Arr_traits>& pgn2)
{
  return (_oriented_side(pgn1, pgn2));
}

template <class Arr_traits, class Gps_traits>
inline Oriented_side oriented_side(const General_polygon_2<Arr_traits>& pgn1,
                                   const General_polygon_2<Arr_traits>& pgn2,
                                   Gps_traits& tr)
{
  return (_oriented_side(pgn1, pgn2, tr));
}

template <class Arr_traits>
inline Oriented_side oriented_side(const General_polygon_2<Arr_traits>& pgn1,
                                   const General_polygon_with_holes_2
                                   <General_polygon_2<Arr_traits> >& pgn2)
{
  return (_oriented_side(pgn1, pgn2));
}

template <class Arr_traits, class Gps_traits>
inline Oriented_side oriented_side(const General_polygon_2<Arr_traits>& pgn1,
                                   const General_polygon_with_holes_2
                                   <General_polygon_2<Arr_traits> >& pgn2,
                                   Gps_traits& tr)
{
  return (_oriented_side(pgn1, pgn2, tr));
}



template <class Arr_traits>
inline Oriented_side oriented_side(const General_polygon_with_holes_2
                                   <General_polygon_2<Arr_traits> >& pgn1,
                                   const General_polygon_2<Arr_traits>& pgn2)
{
  return (_oriented_side(pgn1, pgn2));
}

template <class Arr_traits, class Gps_traits>
inline Oriented_side oriented_side(const General_polygon_with_holes_2
                                   <General_polygon_2<Arr_traits> >& pgn1,
                                   const General_polygon_2<Arr_traits>& pgn2,
                                   Gps_traits& tr)
{
  return (_oriented_side(pgn1, pgn2, tr));
}

template <class Polygon_>
inline
Oriented_side oriented_side(const General_polygon_with_holes_2<Polygon_>& pgn1,
                            const General_polygon_with_holes_2<Polygon_>& pgn2)
{
  return (_oriented_side(pgn1, pgn2));
}

template <class Polygon_, class Traits>
inline
Oriented_side oriented_side(const General_polygon_with_holes_2<Polygon_>& pgn1,
                            const General_polygon_with_holes_2<Polygon_>& pgn2,
                            Traits& tr)
{
  return (_oriented_side(pgn1, pgn2, tr));
}

// Point Query:

template <class Kernel, class Container>
inline Oriented_side oriented_side(const typename Kernel::Point_2& p,
                                   const Polygon_2<Kernel, Container>& pgn)
{
  return (_oriented_side(p, pgn));
}

template <class Kernel, class Container, class Traits>
inline Oriented_side oriented_side(const typename Kernel::Point_2& p,
                                   const Polygon_2<Kernel, Container>& pgn,
                                   Traits& tr)
{
  return (_oriented_side(p, pgn, tr));
}

template <class Kernel, class Container>
inline
Oriented_side oriented_side(const typename Kernel::Point_2& p,
                            const Polygon_with_holes_2<Kernel, Container>& pgn)
{
  return (_oriented_side(p, pgn));
}

template <class Kernel, class Container, class Traits>
inline
Oriented_side oriented_side(const typename Kernel::Point_2& p,
                            const Polygon_with_holes_2<Kernel, Container>& pgn,
                            Traits& tr)
{
  return (_oriented_side(p, pgn, tr));
}

template <class Arr_traits>
inline Oriented_side oriented_side(const typename Arr_traits::Point_2& p,
                                   const General_polygon_2<Arr_traits>& pgn)
{
  return (_oriented_side(p, pgn));
}

template <class Arr_traits, class Gps_traits>
inline Oriented_side oriented_side(const typename Arr_traits::Point_2& p,
                                   const General_polygon_2<Arr_traits>& pgn,
                                   Gps_traits& tr)
{
  return (_oriented_side(p, pgn, tr));
}

template <class Polygon_>
inline
Oriented_side oriented_side(const typename Polygon_::Point_2& p,
                            const General_polygon_with_holes_2<Polygon_>& pgn)
{
  return (_oriented_side(p, pgn));
}

template <class Polygon_, class Traits>
inline
Oriented_side oriented_side(const typename Polygon_::Point_2& p,
                            const General_polygon_with_holes_2<Polygon_>& pgn,
                            Traits& tr)
{
  return (_oriented_side(p, pgn, tr));
}

//@}


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
