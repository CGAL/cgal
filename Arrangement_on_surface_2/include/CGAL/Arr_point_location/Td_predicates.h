// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// Author(s)	 : Oren Nechushtan <theoren@math.tau.ac.il>
//		   Iddo Hanniel <hanniel@math.tau.ac.il>
#ifndef CGAL_TD_PREDICATES_H
#define CGAL_TD_PREDICATES_H

#include <CGAL/basic.h>

#include <functional>

namespace CGAL {

template < class Td_traits> class Trapezoidal_decomposition_2;

template <class X_trapezoid>
struct Td_active_trapezoid : public std::unary_function<X_trapezoid,bool>
{
  bool operator()(const X_trapezoid& tr) const
  {
    return tr.is_active();
  }
};

template <class X_trapezoid,class Traits>
struct Td_active_non_degenerate_trapezoid : 
public std::unary_function<X_trapezoid,bool>
{
  Td_active_non_degenerate_trapezoid(Traits& t) : traits(t) {}
  bool operator()(const X_trapezoid& tr) const
  {
    return tr.is_active() && !traits.is_degenerate(tr);
  }
protected:
  const Traits& traits;
};

template <class X_trapezoid,class Traits>
struct Td_active_right_degenerate_curve_trapezoid:
  public std::unary_function<X_trapezoid,bool>
{
  typedef const Traits& const_Traits_ref;
  Td_active_right_degenerate_curve_trapezoid(const_Traits_ref t) : traits(t) {}
  bool operator()(const X_trapezoid& tr) const
  {
    return tr.is_active() && traits.is_degenerate_curve(tr) && 
      !tr.right_bottom_neighbour();
  }
  protected:
  const Traits& traits;
};

template <class _Tp>
struct Trapezoid_handle_less : public std::binary_function<_Tp, _Tp, bool>
{
  bool operator()(const _Tp& __x, const _Tp& __y) const { 
    return __x->id() < __y->id(); }
};
/* Return if two trapezoids are the same */

} //namespace CGAL

#endif //CGAL_TD_PREDICATES_H
