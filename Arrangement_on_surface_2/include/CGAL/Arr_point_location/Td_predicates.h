// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)	 : Oren Nechushtan <theoren@math.tau.ac.il>
//		   Iddo Hanniel <hanniel@math.tau.ac.il>
#ifndef CGAL_TD_PREDICATES_H
#define CGAL_TD_PREDICATES_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <CGAL/basic.h>

#include <functional>

namespace CGAL {

template < class Td_traits> class Trapezoidal_decomposition_2;

////MICHAL: not in use
//template <class map_item>
//struct Td_active_map_item : public std::unary_function<map_item,bool>
//{
//  bool operator()(const map_item& item) const
//  {
//    return item.is_active();
//  }
//};

////MICHAL: not in use
//template <class X_trapezoid,class Traits>
//struct Td_active_non_degenerate_trapezoid : 
//public std::unary_function<X_trapezoid,bool>
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
  public std::unary_function<map_item,bool>
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
struct Td_map_item_handle_less : public std::binary_function<_Tp, _Tp, bool>
{
  bool operator()(const _Tp& __x, const _Tp& __y) const { 
    return __x->id() < __y->id(); }
};
/* Return if two trapezoids are the same */

} //namespace CGAL

#endif //CGAL_TD_PREDICATES_H
