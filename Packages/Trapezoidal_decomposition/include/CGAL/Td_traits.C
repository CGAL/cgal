// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/Td_traits.C
// package       : Trapezoidal decomposition 2
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// maintainer(s) : Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_TD_TRAITS_H
#include <CGAL/Td_traits.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Traits,class X_curve_plus>
const typename Td_traits<Traits,X_curve_plus>::Point&
Td_traits<Traits,X_curve_plus>::get_point_at_left_top_infinity(){
  /*
  static typename Td_traits<Traits,X_curve_plus>::Point
    Td_traits<Traits,X_curve_plus>::POINT_AT_LEFT_TOP_INFINITY;
  return Td_traits<Traits,X_curve_plus>::POINT_AT_LEFT_TOP_INFINITY;
  */
  //  static Point POINT_AT_LEFT_TOP_INFINITY;
  return POINT_AT_LEFT_TOP_INFINITY;
}

template <class Traits,class X_curve_plus>
const typename Td_traits<Traits,X_curve_plus>::Point&
Td_traits<Traits,X_curve_plus>::get_point_at_right_bottom_infinity(){
  /*
  static typename Td_traits<Traits,X_curve_plus>::Point
    Td_traits<Traits,X_curve_plus>::POINT_AT_RIGHT_BOTTOM_INFINITY;
  return Td_traits<Traits,X_curve_plus>::POINT_AT_RIGHT_BOTTOM_INFINITY;
  */
  //  static Point POINT_AT_RIGHT_BOTTOM_INFINITY;
  return POINT_AT_RIGHT_BOTTOM_INFINITY;
}

template <class Traits,class X_curve_plus>
const typename Td_traits<Traits,X_curve_plus>::X_curve&
Td_traits<Traits,X_curve_plus>::get_curve_at_infinity(){
  /*
  static typename typename Traits::X_curveTraits::X_curve 
    Td_traits<Traits,X_curve_plus>::CURVE_AT_INFINITY;
  return Td_traits<Traits,X_curve_plus>::CURVE_AT_INFINITY;
  */
  //  static X_curve CURVE_AT_INFINITY;
  return CURVE_AT_INFINITY;
}

template <class Traits,class X_curve_plus>
typename Td_traits<Traits,X_curve_plus>::Point
Td_traits<Traits,X_curve_plus>::POINT_AT_LEFT_TOP_INFINITY = 
Point(0,0);

template <class Traits,class X_curve_plus>
typename Td_traits<Traits,X_curve_plus>::Point
Td_traits<Traits,X_curve_plus>::POINT_AT_RIGHT_BOTTOM_INFINITY = 
Point(0,0);


template <class Traits,class X_curve_plus>
typename Td_traits<Traits,X_curve_plus>::X_curve
Td_traits<Traits,X_curve_plus>::CURVE_AT_INFINITY = X_curve();

CGAL_END_NAMESPACE
