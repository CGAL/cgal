// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein        <wein@post.tau.ac.il>
//                 Efi Fogel       <efif@post.tau.ac.il>

#ifndef CGAL_BOOLEAN_SET_OPERATIONS_H
#define CGAL_BOOLEAN_SET_OPERATIONS_H

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/iterator.h> 
#include <CGAL/Boolean_set_operations_2/Bso_internal_functions.h>

CGAL_BEGIN_NAMESPACE

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

template <class Kernel, class Container, class OutputIterator>
inline OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1, 
                                   const Polygon_2<Kernel, Container>& pgn2,
                                   OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Kernel, class Container, class OutputIterator, class Traits>
inline OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1, 
                                   const Polygon_2<Kernel, Container>& pgn2,
                                   OutputIterator out,
                                   Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

template <class Kernel, class Container, class OutputIterator>
inline OutputIterator intersection(const Polygon_2<Kernel, Container>& pgn1, 
                         const Polygon_with_holes_2<Kernel, Container>& pgn2,
                           OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Kernel, class Container, class OutputIterator, class Traits>
inline OutputIterator
intersection (const Polygon_2<Kernel, Container>& pgn1, 
              const Polygon_with_holes_2<Kernel, Container>& pgn2,
              OutputIterator out, Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

template <class Kernel, class Container, class OutputIterator>
inline OutputIterator
intersection (const Polygon_with_holes_2<Kernel, Container>& pgn1,
              const Polygon_2<Kernel, Container>& pgn2,
              OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Kernel, class Container, class OutputIterator, class Traits>
inline OutputIterator
intersection (const Polygon_with_holes_2<Kernel, Container>& pgn1,
              const Polygon_2<Kernel, Container>& pgn2,
              OutputIterator out, Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

template <class Kernel, class Container, class OutputIterator>
inline OutputIterator 
intersection (const Polygon_with_holes_2<Kernel, Container>& pgn1,
              const Polygon_with_holes_2<Kernel, Container>& pgn2,
              OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Kernel, class Container, class OutputIterator, class Traits>
inline OutputIterator
intersection (const Polygon_with_holes_2<Kernel, Container>& pgn1,
              const Polygon_with_holes_2<Kernel, Container>& pgn2,
              OutputIterator out,
              Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

template <class Arr_traits, class OutputIterator>
inline OutputIterator intersection (const General_polygon_2<Arr_traits>& pgn1, 
                                    const General_polygon_2<Arr_traits>& pgn2,
                                    OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Arr_traits, class OutputIterator, class Traits>
inline OutputIterator intersection (const General_polygon_2<Arr_traits>& pgn1, 
                                    const General_polygon_2<Arr_traits>& pgn2,
                                    OutputIterator out,
                                    Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

template <class Arr_traits, class OutputIterator>
inline OutputIterator
intersection (const General_polygon_2<Arr_traits>& pgn1, 
              const General_polygon_with_holes_2
                    <General_polygon_2<Arr_traits> >& pgn2,
              OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Arr_traits, class OutputIterator, class Traits>
inline OutputIterator
intersection (const General_polygon_2<Arr_traits>& pgn1, 
              const General_polygon_with_holes_2
                    <General_polygon_2<Arr_traits> >& pgn2,
              OutputIterator out, Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

template <class Arr_traits, class OutputIterator>
inline OutputIterator
intersection (const General_polygon_with_holes_2
                    <General_polygon_2<Arr_traits> >& pgn1,
              const General_polygon_2<Arr_traits>& pgn2,
              OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Arr_traits, class OutputIterator, class Traits>
inline OutputIterator
intersection (const General_polygon_with_holes_2
                    <General_polygon_2<Arr_traits> >& pgn1,
              const General_polygon_2<Arr_traits>& pgn2,
              OutputIterator out, Traits& tr)
{
  return (_intersection(pgn1, pgn2, out, tr));
}

template <class Polygon_, class OutputIterator>
inline OutputIterator
intersection (const General_polygon_with_holes_2<Polygon_>& pgn1,
              const General_polygon_with_holes_2<Polygon_>& pgn2,
              OutputIterator out)
{
  return (_intersection(pgn1, pgn2, out));
}

template <class Polygon_, class OutputIterator, class Traits>
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

template <class Kernel, class Container, class OutputIterator>
inline OutputIterator difference (const Polygon_2<Kernel, Container>& pgn1, 
                                  const Polygon_2<Kernel, Container>& pgn2,
                                  OutputIterator oi)
{
  return(_difference(pgn1, pgn2, oi));
}

template <class Kernel, class Container, class OutputIterator, class Traits>
inline OutputIterator difference (const Polygon_2<Kernel, Container>& pgn1, 
                                  const Polygon_2<Kernel, Container>& pgn2,
                                  OutputIterator oi, Traits& tr)
{
  return(_difference(pgn1, pgn2, oi, tr));
}

template <class Kernel, class Container, class OutputIterator>
inline OutputIterator
difference (const Polygon_2<Kernel, Container>& pgn1, 
            const Polygon_with_holes_2<Kernel, Container>& pgn2,
            OutputIterator oi)
{
  return(_difference(pgn1, pgn2, oi));
}

template <class Kernel, class Container, class OutputIterator, class Traits>
inline OutputIterator
difference (const Polygon_2<Kernel, Container>& pgn1, 
            const Polygon_with_holes_2<Kernel, Container>& pgn2,
            OutputIterator oi, Traits& tr)
{
  return(_difference(pgn1, pgn2, oi, tr));
}

template <class Kernel, class Container, class OutputIterator>
inline OutputIterator
difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
            const Polygon_2<Kernel, Container>& pgn2,
            OutputIterator oi)
{
  return (_difference(pgn1, pgn2, oi));
}

template <class Kernel, class Container, class OutputIterator, class Traits>
inline OutputIterator
difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
            const Polygon_2<Kernel, Container>& pgn2,
            OutputIterator oi,
            Traits& tr)
{
  return (_difference(pgn1, pgn2, oi, tr));
}

template <class Kernel, class Container, class OutputIterator>
inline OutputIterator
difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
            const Polygon_with_holes_2<Kernel, Container>& pgn2,
            OutputIterator oi)
{
  return (_difference(pgn1, pgn2, oi));
}

template <class Kernel, class Container, class OutputIterator, class Traits>
inline OutputIterator
difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
            const Polygon_with_holes_2<Kernel, Container>& pgn2,
            OutputIterator oi, Traits& tr)
{
  return (_difference(pgn1, pgn2, oi, tr));
}

template <class Arr_traits, class OutputIterator>
inline OutputIterator difference (const General_polygon_2<Arr_traits>& pgn1, 
                                  const General_polygon_2<Arr_traits>& pgn2,
                                  OutputIterator oi)
{
  return (_difference(pgn1, pgn2, oi));
}

template <class Arr_traits, class OutputIterator, class Traits>
inline OutputIterator difference (const General_polygon_2<Arr_traits>& pgn1, 
                                  const General_polygon_2<Arr_traits>& pgn2,
                                  OutputIterator oi,
                                  Traits& tr)
{
  return (_difference(pgn1, pgn2, oi, tr));
}

template <class Arr_traits, class OutputIterator>
inline OutputIterator difference (const General_polygon_2<Arr_traits>& pgn1, 
                                  const General_polygon_with_holes_2
                                        <General_polygon_2<Arr_traits> >& pgn2,
                                  OutputIterator oi)
{
  return (_difference(pgn1, pgn2, oi));
}

template <class Arr_traits, class OutputIterator, class Traits>
inline OutputIterator difference (const General_polygon_2<Arr_traits>& pgn1, 
                                  const General_polygon_with_holes_2
                                        <General_polygon_2<Arr_traits> >& pgn2,
                                  OutputIterator oi, Traits& tr)
{
  return (_difference(pgn1, pgn2, oi, tr));
}

template <class Arr_traits, class OutputIterator>
inline OutputIterator difference (const General_polygon_with_holes_2
                                        <General_polygon_2<Arr_traits> >& pgn1,
                                  const General_polygon_2<Arr_traits>& pgn2,
                                  OutputIterator oi)
{
  return (_difference(pgn1, pgn2, oi));
}

template <class Arr_traits, class OutputIterator, class Traits>
inline OutputIterator difference (const General_polygon_with_holes_2
                                        <General_polygon_2<Arr_traits> >& pgn1,
                                  const General_polygon_2<Arr_traits>& pgn2,
                                  OutputIterator oi, Traits& tr)
{
  return (_difference(pgn1, pgn2, oi, tr));
}

template <class Polygon_, class OutputIterator>
inline OutputIterator
difference (const General_polygon_with_holes_2<Polygon_>& pgn1,
            const General_polygon_with_holes_2<Polygon_>& pgn2,
            OutputIterator oi)
{
  return (_difference(pgn1, pgn2, oi));
}

template <class Polygon_, class OutputIterator, class Traits>
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

template <class Kernel, class Container, class OutputIterator, class Traits>
inline OutputIterator
symmetric_difference (const Polygon_2<Kernel, Container>& pgn1, 
                      const Polygon_2<Kernel, Container>& pgn2,
                      OutputIterator oi, Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Kernel, class Container, class OutputIterator>
inline OutputIterator
symmetric_difference (const Polygon_2<Kernel, Container>& pgn1, 
                      const Polygon_2<Kernel, Container>& pgn2,
                      OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

template <class Kernel, class Container, class OutputIterator, class Traits>
inline OutputIterator
symmetric_difference (const Polygon_2<Kernel, Container>& pgn1, 
                      const Polygon_with_holes_2<Kernel, Container>& pgn2,
                      OutputIterator oi, Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}


template <class Kernel, class Container, class OutputIterator>
inline OutputIterator
symmetric_difference (const Polygon_2<Kernel, Container>& pgn1, 
                      const Polygon_with_holes_2<Kernel, Container>& pgn2,
                      OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

template <class Kernel, class Container, class OutputIterator, class Traits>
inline OutputIterator
symmetric_difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
                      const Polygon_2<Kernel, Container>& pgn2,
                      OutputIterator oi, Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Kernel, class Container, class OutputIterator>
inline OutputIterator
symmetric_difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
                      const Polygon_2<Kernel, Container>& pgn2,
                      OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

template <class Kernel, class Container, class OutputIterator, class Traits>
inline OutputIterator
symmetric_difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
                      const Polygon_with_holes_2<Kernel, Container>& pgn2,
                      OutputIterator oi, Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Kernel, class Container, class OutputIterator>
inline OutputIterator
symmetric_difference (const Polygon_with_holes_2<Kernel, Container>& pgn1,
                      const Polygon_with_holes_2<Kernel, Container>& pgn2,
                      OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

template <class Arr_traits, class OutputIterator, class Traits>
inline OutputIterator 
symmetric_difference (const General_polygon_2<Arr_traits>& pgn1, 
                      const General_polygon_2<Arr_traits>& pgn2,
                      OutputIterator oi,
                      Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Arr_traits, class OutputIterator>
inline OutputIterator 
symmetric_difference (const General_polygon_2<Arr_traits>& pgn1, 
                      const General_polygon_2<Arr_traits>& pgn2,
                      OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

template <class Arr_traits, class OutputIterator, class Traits>
inline OutputIterator 
symmetric_difference (const General_polygon_2<Arr_traits>& pgn1, 
                      const General_polygon_with_holes_2
                            <General_polygon_2<Arr_traits> >& pgn2,
                      OutputIterator oi, Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Arr_traits, class OutputIterator>
inline OutputIterator 
symmetric_difference(const General_polygon_2<Arr_traits>& pgn1, 
                     const General_polygon_with_holes_2
                           <General_polygon_2<Arr_traits> >& pgn2,
                     OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

template <class Arr_traits, class OutputIterator, class Traits>
inline OutputIterator
symmetric_difference (const General_polygon_with_holes_2
                        <General_polygon_2<Arr_traits> >& pgn1,
                      const General_polygon_2<Arr_traits>& pgn2,
                      OutputIterator oi, Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Arr_traits, class OutputIterator>
inline OutputIterator
symmetric_difference (const General_polygon_with_holes_2
                            <General_polygon_2<Arr_traits> >& pgn1,
                      const General_polygon_2<Arr_traits>& pgn2,
                      OutputIterator oi)
{
  return (_symmetric_difference(pgn1, pgn2, oi));
}

template <class Polygon_, class OutputIterator, class Traits>
inline OutputIterator
symmetric_difference (const General_polygon_with_holes_2<Polygon_>& pgn1,
                      const General_polygon_with_holes_2<Polygon_>& pgn2,
                      OutputIterator oi,
                      Traits& tr)
{
  return (_symmetric_difference(pgn1, pgn2, oi, tr));
}

template <class Polygon_, class OutputIterator>
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

template <class Kernel, class Container, class OutputIterator, class Traits>
OutputIterator complement (const Polygon_with_holes_2<Kernel, Container>& pgn,
                           OutputIterator oi, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn);
  gps.complement();
  return (gps.polygons_with_holes(oi));
}

template <class Kernel, class Container, class OutputIterator>
OutputIterator complement (const Polygon_with_holes_2<Kernel, Container>& pgn,
                           OutputIterator oi)
{
  typename Gps_default_traits<Polygon_2<Kernel, Container> >::Traits  tr;
  return (complement(pgn, oi, tr));
}

template <class Arr_traits, class OutputIterator, class Traits>
OutputIterator complement (const General_polygon_with_holes_2<Arr_traits>& pgn,
                           OutputIterator oi, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn);
  gps.complement();
  return (gps.polygons_with_holes(oi));
}

template <class Arr_traits, class OutputIterator>
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

template <class InputIterator>
struct map_iterator_to_traits
{
  typedef typename std::iterator_traits<InputIterator>::value_type InputPolygon;
  typedef typename Gps_default_traits<InputPolygon>::Traits    Traits;
};

template <class InputIterator, class OutputIterator, class Traits>
inline OutputIterator join(InputIterator begin, InputIterator end,
                           OutputIterator oi, Traits& )
{
  if (begin == end)
    return (oi);

  General_polygon_set_2<Traits> gps(*begin);
  gps.join(++begin, end);
  return (gps.polygons_with_holes(oi));
}

template <class InputIterator, class OutputIterator>
inline OutputIterator join(InputIterator begin, InputIterator end,
                           OutputIterator oi)
{
  typename map_iterator_to_traits<InputIterator>::Traits          tr;
  return join(begin, end, oi, tr);
}

// Join two ranges of simple polygons and polygons with holes.
template <class InputIterator1, class InputIterator2,
          class OutputIterator, class Traits>
inline OutputIterator join(InputIterator1 begin1, InputIterator1 end1,
                           InputIterator2 begin2, InputIterator2 end2,
                           OutputIterator oi, Traits& tr)
{
  if (begin1 == end1)
    return (join(begin2, end2, oi, tr));

  General_polygon_set_2<Traits> gps(*begin1);
  gps.join(++begin1, end1, begin2, end2);
  return (gps.polygons_with_holes(oi));
}

template <class InputIterator1, class InputIterator2, class OutputIterator>
inline OutputIterator join(InputIterator1 begin1, InputIterator1 end1,
                           InputIterator2 begin2, InputIterator2 end2,
                           OutputIterator oi)
{
  typename map_iterator_to_traits<InputIterator1>::Traits  tr;
  return join(begin1, end1, begin2, end2, oi, tr);
}

//@}
/// \name Aggregated intersection() functions.
//@{
template <class InputIterator, class OutputIterator, class Traits>
inline OutputIterator intersection (InputIterator begin, InputIterator end,
                                    OutputIterator oi, Traits&)
{
  if (begin == end)
    return (oi);

  General_polygon_set_2<Traits> gps(*begin);
  gps.intersection(++begin, end);
  return (gps.polygons_with_holes(oi));
}

template <class InputIterator, class OutputIterator>
inline OutputIterator intersection (InputIterator begin, InputIterator end,
                                    OutputIterator oi)
{
  typename map_iterator_to_traits<InputIterator>::Traits          tr;
  return intersection(begin, end, oi, tr);
}

// Inersect two ranges of simple polygons and polygons with holes.
template <class InputIterator1, class InputIterator2,
          class OutputIterator, class Traits>
inline OutputIterator intersection(InputIterator1 begin1, InputIterator1 end1,
                                   InputIterator2 begin2, InputIterator2 end2,
                                   OutputIterator oi, Traits& tr)
{
  if (begin1 == end1)
    return (intersection(begin2, end2, oi, tr));

  General_polygon_set_2<Traits> gps(*begin1);
  gps.intersection(++begin1, end1, begin2, end2);
  return (gps.polygons_with_holes(oi));
 
}

template <class InputIterator1, class InputIterator2, class OutputIterator>
inline OutputIterator intersection(InputIterator1 begin1, InputIterator1 end1,
                                   InputIterator2 begin2, InputIterator2 end2,
                                   OutputIterator oi)
{
  typename map_iterator_to_traits<InputIterator1>::Traits  tr;
  return intersection(begin1, end1, begin2, end2, oi, tr);
}

//@}
/// \name Aggregated symmetric_difference() functions.
//@{
template <class InputIterator, class OutputIterator, class Traits>
inline
OutputIterator symmetric_difference(InputIterator begin, InputIterator end,
                                    OutputIterator oi, Traits& tr)
{
  if (begin == end)
    return (oi);

  General_polygon_set_2<Traits> gps(tr);
  gps.insert(*begin);
  gps.symmetric_difference(++begin, end);
  return (gps.polygons_with_holes(oi));
}

template <class InputIterator, class OutputIterator>
inline
OutputIterator symmetric_difference (InputIterator begin, InputIterator end,
                                     OutputIterator oi)
{
  typename map_iterator_to_traits<InputIterator>::Traits          tr;
  return symmetric_difference(begin, end, oi, tr);
}

// Xor two ranges of simple polygons and polygons with holes.
template <class InputIterator1, class InputIterator2,
          class OutputIterator, class Traits>
inline
OutputIterator symmetric_difference (InputIterator1 begin1, InputIterator1 end1,
                                     InputIterator2 begin2, InputIterator2 end2,
                                     OutputIterator oi, Traits& tr)
{
  if (begin1 == end1)
    return (symmetric_difference(begin2, end2, oi, tr));

  General_polygon_set_2<Traits> gps(tr);
  gps.insert(*begin1);
  gps.symmetric_difference(++begin1, end1, begin2, end2);
  return (gps.polygons_with_holes(oi));
 
}

template <class InputIterator1, class InputIterator2, class OutputIterator>
inline
OutputIterator symmetric_difference (InputIterator1 begin1, InputIterator1 end1,
                                     InputIterator2 begin2, InputIterator2 end2,
                                     OutputIterator oi)
{
  typename map_iterator_to_traits<InputIterator1>::Traits  tr;
  return symmetric_difference(begin1, end1, begin2, end2, oi, tr);
}

//@}
/// \name Aggregated do_intersect() functions.
//@{

template <class InputIterator, class Traits>
inline bool do_intersect(InputIterator begin, InputIterator end, Traits& tr)
{
  if (begin == end)
    return false;

  General_polygon_set_2<Traits> gps(tr);
  gps.insert(*begin);

  return (gps.do_intersect(++begin, end));
}

template <class InputIterator>
inline bool do_intersect(InputIterator begin, InputIterator end)
{
  typename map_iterator_to_traits<InputIterator>::Traits  tr;
  return do_intersect(begin, end, tr);
}

template <class InputIterator1, class InputIterator2, class Traits>
inline bool do_intersect (InputIterator1 begin1, InputIterator1 end1,
                          InputIterator2 begin2, InputIterator2 end2,
                          Traits& tr)
{
  if (begin1 == end1)
    return (do_intersect(begin2, end2, tr));

  General_polygon_set_2<Traits> gps(tr);
  gps.insert(*begin1);
  return (gps.do_intersect(++begin1, end1, begin2, end2));
}

template <class InputIterator1, class InputIterator2>
inline bool do_intersect (InputIterator1 begin1, InputIterator1 end1,
                          InputIterator2 begin2, InputIterator2 end2)
{
  typename map_iterator_to_traits<InputIterator1>::Traits  tr;
  return do_intersect(begin1, end1, begin2, end2, tr);
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


CGAL_END_NAMESPACE

#endif
