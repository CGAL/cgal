// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)         : Oren Nechushtan <theoren@math.tau.ac.il>
//                   Iddo Hanniel <hanniel@math.tau.ac.il>
#ifndef CGAL_TD_PREDICATES_H
#define CGAL_TD_PREDICATES_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <CGAL/basic.h>

#include <functional>

namespace CGAL {

template < class Td_traits> class Trapezoidal_decomposition_2;

////MICHAL: not in use
//template <class map_item>
//struct Td_active_map_item : public CGAL::cpp98::unary_function<map_item,bool>
//{
//  bool operator()(const map_item& item) const
//  {
//    return item.is_active();
//  }
//};

////MICHAL: not in use
//template <class X_trapezoid,class Traits>
//struct Td_active_non_degenerate_trapezoid :
//public CGAL::cpp98::unary_function<X_trapezoid,bool>
//{
//  Td_active_non_degenerate_trapezoid(Traits& t) : traits(t) {}
//  bool operator()(const X_trapezoid& tr) const
//  {
//    return tr.is_active() && !traits.is_degenerate(tr);
//  }
//protected:
//  const Traits& traits;
//};

template <class map_item,class Traits>
struct Td_active_edge_item:
  public CGAL::cpp98::unary_function<map_item,bool>
{
  Td_active_edge_item(const Traits& t) : traits(t) {}
  bool operator()(const map_item& item) const
  {
    return traits.is_active(item) && traits.is_td_edge(item);
  }
  protected:
  const Traits& traits;
};

template <class _Tp>
struct Td_map_item_handle_less : public CGAL::cpp98::binary_function<_Tp, _Tp, bool>
{
  bool operator()(const _Tp& __x, const _Tp& __y) const {
    return __x->id() < __y->id(); }
};
/* Return if two trapezoids are the same */

} //namespace CGAL

#endif //CGAL_TD_PREDICATES_H
