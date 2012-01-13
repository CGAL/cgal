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
// Author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
#ifndef CGAL_TD_TRAITS_H
#include <CGAL/Arr_point_location/Td_traits.h>
#endif

namespace CGAL {

template <class Traits,class X_curve_plus>
const typename Td_traits<Traits,X_curve_plus>::Point&
Td_traits<Traits,X_curve_plus>::point_at_left_top_infinity(){
  /*
  static typename Td_traits<Traits,X_curve_plus>::Point
    Td_traits<Traits,X_curve_plus>::POINT_AT_LEFT_TOP_INFINITY;
  return Td_traits<Traits,X_curve_plus>::POINT_AT_LEFT_TOP_INFINITY;
  */
  //  static Point POINT_AT_LEFT_TOP_INFINITY;
  if (!POINT_AT_LEFT_TOP_INFINITY)
    POINT_AT_LEFT_TOP_INFINITY = new Point();
  return *POINT_AT_LEFT_TOP_INFINITY;
}

template <class Traits,class X_curve_plus>
const typename Td_traits<Traits,X_curve_plus>::Point&
Td_traits<Traits,X_curve_plus>::point_at_right_bottom_infinity(){
  /*
  static typename Td_traits<Traits,X_curve_plus>::Point
    Td_traits<Traits,X_curve_plus>::POINT_AT_RIGHT_BOTTOM_INFINITY;
  return Td_traits<Traits,X_curve_plus>::POINT_AT_RIGHT_BOTTOM_INFINITY;
  */
  //  static Point POINT_AT_RIGHT_BOTTOM_INFINITY;
  if (!POINT_AT_RIGHT_BOTTOM_INFINITY)
    POINT_AT_RIGHT_BOTTOM_INFINITY = new Point();
  return *POINT_AT_RIGHT_BOTTOM_INFINITY;
}

template <class Traits,class X_curve_plus>
const typename Td_traits<Traits,X_curve_plus>::X_curve&
Td_traits<Traits,X_curve_plus>::curve_at_infinity(){
  /*
  static typename typename Traits::X_curveTraits::X_curve 
    Td_traits<Traits,X_curve_plus>::CURVE_AT_INFINITY;
  return Td_traits<Traits,X_curve_plus>::CURVE_AT_INFINITY;
  */
  //  static X_curve CURVE_AT_INFINITY;
  if (!CURVE_AT_INFINITY)
    CURVE_AT_INFINITY = new X_curve();
  return *CURVE_AT_INFINITY;
}

template <class Traits,class X_curve_plus>
typename Td_traits<Traits,X_curve_plus>::Point *
Td_traits<Traits,X_curve_plus>::POINT_AT_LEFT_TOP_INFINITY = 0;

template <class Traits,class X_curve_plus>
typename Td_traits<Traits,X_curve_plus>::Point *
Td_traits<Traits,X_curve_plus>::POINT_AT_RIGHT_BOTTOM_INFINITY = 0;

template <class Traits,class X_curve_plus>
typename Td_traits<Traits,X_curve_plus>::X_curve *
Td_traits<Traits,X_curve_plus>::CURVE_AT_INFINITY = 0;

} //namespace CGAL
