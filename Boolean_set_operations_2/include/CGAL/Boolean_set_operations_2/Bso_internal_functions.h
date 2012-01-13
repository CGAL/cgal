// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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

#ifndef CGAL_BSO_INTERNAL_FUNCTIONS_H
#define CGAL_BSO_INTERNAL_FUNCTIONS_H

#include <CGAL/Boolean_set_operations_2/Gps_default_traits.h>
#include <iterator>

namespace CGAL {

/// \name _do_intersect() functions.
//@{

template <class Pgn1, class Pgn2, class Traits>
inline bool _do_intersect(const Pgn1& pgn1, const Pgn2& pgn2, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn1);
  return (gps.do_intersect(pgn2));
}

template <class Pgn1, class Pgn2>
inline bool _do_intersect(const Pgn1& pgn1, const Pgn2& pgn2)
{
  typename Gps_default_traits<Pgn1>::Traits    tr;
  return _do_intersect(pgn1, pgn2, tr);
}

//@}
/// \name _oriented_side() functions.
//@{

template <class Obj, class Pgn, class Traits>
inline
Oriented_side _oriented_side(const Obj& obj, const Pgn& pgn, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn);
  return (gps.oriented_side(obj));
}

template <class Obj, class Pgn>
inline Oriented_side _oriented_side(const Obj& obj, const Pgn& pgn)
{
  typename Gps_default_traits<Pgn>::Traits    tr;
  return _oriented_side(obj, pgn, tr);
}

//@}
/// \name _intersection() functions.
//@{

template <class Pgn1, class Pgn2, class OutputIterator, class Traits>
inline OutputIterator _intersection(const Pgn1& pgn1, const Pgn2& pgn2,
                                    OutputIterator out, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn1);
  gps.intersection(pgn2);
  return (gps.polygons_with_holes(out));
}

template <class Pgn1, class Pgn2, class OutputIterator>
inline OutputIterator _intersection(const Pgn1& pgn1, const Pgn2& pgn2,
                                    OutputIterator out)
{
  typename Gps_default_traits<Pgn1>::Traits    tr;
  return (_intersection(pgn1, pgn2, out, tr));
}

//@}
/// \name _join() functions.
//@{

template <class Traits>
inline bool _is_empty (const typename Traits:: Polygon_2& pgn, Traits& tr) 
{
  typedef typename Traits::Curve_const_iterator Curve_const_iterator;
  const std::pair<Curve_const_iterator, Curve_const_iterator>& itr_pair =
    tr.construct_curves_2_object()(pgn);
  return (itr_pair.first == itr_pair.second);
}

template <class Traits>
inline bool _is_empty (const typename Traits::Polygon_with_holes_2&, Traits&)
{
  return false;
}

template <class Pgn1, class Pgn2, class Traits>
inline bool _join(const Pgn1& pgn1, const Pgn2& pgn2,
                  typename Traits::Polygon_with_holes_2& res, Traits& tr)
{
  if (_is_empty(pgn1, tr) || _is_empty(pgn2, tr))
    return false;

  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn1);
  gps.join(pgn2);
  if (gps.number_of_polygons_with_holes() == 1)
  {
    Oneset_iterator<typename Traits::Polygon_with_holes_2> oi (res);
    gps.polygons_with_holes(oi);
    return true;
  }

  // the polygon doesnt intersect, the original pgn1, pgn2 contain the union
  return false;
}

template <class Pgn1, class Pgn2, class Pwh>
inline bool _join(const Pgn1& pgn1, const Pgn2& pgn2, Pwh& res)
{
  typename Gps_default_traits<Pgn1>::Traits  tr;
  return _join(pgn1, pgn2, res, tr);
}

//@}
/// \name _difference() functions.
//@{

template <class Pgn1, class Pgn2, class OutputIterator, class Traits>
inline OutputIterator _difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                  OutputIterator oi, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn1);
  gps.difference(pgn2);
  return gps.polygons_with_holes(oi);
}

template <class Pgn1, class Pgn2, class OutputIterator>
inline OutputIterator _difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                  OutputIterator oi)
{
  typename Gps_default_traits<Pgn1>::Traits  tr;
  return _difference(pgn1, pgn2, oi, tr);
}

//@}
/// \name _symmetric_difference() functions.
//@{

template <class Pgn1, class Pgn2, class OutputIterator, class Traits>
inline OutputIterator _symmetric_difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                            OutputIterator oi, Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn1);
  gps.symmetric_difference(pgn2);
  return gps.polygons_with_holes(oi);
}

template <class Pgn1, class Pgn2, class OutputIterator>
inline OutputIterator _symmetric_difference(const Pgn1& pgn1, const Pgn2& pgn2,
                                            OutputIterator oi)
{
  typename Gps_default_traits<Pgn1>::Traits    tr;
  return _symmetric_difference(pgn1, pgn2, oi, tr);
}

//@}
/// \name _complement() functions.
//@{

template <class Pgn, class Traits>
void _complement(const Pgn& pgn, typename Traits::Polygon_with_holes_2& res,
                 Traits& tr)
{
  General_polygon_set_2<Traits> gps(tr);
  gps.insert(pgn);
  gps.complement();
  Oneset_iterator<typename Traits::Polygon_with_holes_2> oi(res);
  gps.polygons_with_holes(oi);
}

template <class Pgn, class Pwh>
void _complement(const Pgn& pgn, Pwh& res)
{
  typename Gps_default_traits<Pgn>::Traits    tr;
  _complement(pgn, res, tr);
}

//@}

} //namespace CGAL

#endif
