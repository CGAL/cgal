// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Sebastien Loriot, Sylvain Pion

#ifndef CGAL_TUPLE_H
#define CGAL_TUPLE_H

// A wrapper around C++0x, TR1 or Boost tuple<>.
// Together with the Is_in_tuple<> tool.

#include <CGAL/config.h>

// workaround issue in the tuple implementation of libc++
// by falling back onto boost's
#ifndef CGAL_CFG_NO_CPP0X_TUPLE
  #include <cstddef>
  #if defined( _LIBCPP_VERSION ) // check if libc++ is used
    #define CGAL_CFG_NO_CPP0X_TUPLE
  #endif
  #if !defined(__clang__) && defined(__GNUC__) && defined(__GNUC_MINOR__) && (__GNUC__ <= 4) && (__GNUC_MINOR__ <= 4)
    #define CGAL_CFG_NO_CPP0X_TUPLE
    #define CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  #endif
#endif

#ifndef CGAL_CFG_NO_CPP0X_TUPLE
#  include <tuple>
#else
#  include <boost/tuple/tuple.hpp>
#  include <boost/tuple/tuple_comparison.hpp>
#  include <utility>
#endif

namespace CGAL {

namespace cpp11 {

#ifndef CGAL_CFG_NO_CPP0X_TUPLE
using std::tuple;
using std::make_tuple;
using std::tie;
using std::get;
using std::tuple_size;
using std::tuple_element;
#else
using boost::tuple;
using boost::make_tuple;
using boost::tie;
using boost::get;
  
//tuple_size
template <class T>
struct tuple_size:public boost::tuples::length<T> {};
  
//tuple_element
template <int N,class T>
struct tuple_element: public boost::tuples::element<N,T>{};
  
#endif


#if defined(CGAL_CFG_NO_CPP0X_TUPLE)
// If not C++11 tuple, we need to add get<N>(std::pair).

////////////////////////////////////////////////////////////
//                                                        //
// Allow CGAL::cpp11::get<N>(std::pair), if N==0 or N==1. //
//                                                        //
// That is already in TR1 and C++11, but not in Boost.    //
//                                                        //
////////////////////////////////////////////////////////////
template <std::size_t N, typename T1, typename T2>
struct pair_get;

template <typename T1, typename T2>
struct pair_get<0, T1, T2> {
  static T1& get(std::pair<T1, T2>& pair) { return pair.first; }
  static const T1& get(const std::pair<T1, T2>& pair) { return pair.first; }
}; // end specialization struct pair_get<0, T2, T2>

template <typename T1, typename T2>
struct pair_get<1, T1, T2> {
  static T2& get(std::pair<T1, T2>& pair) { return pair.second; }
  static const T2& get(const std::pair<T1, T2>& pair) { return pair.second; }
}; // end specialization struct pair_get<0, T2, T2>

template <std::size_t N, typename T1, typename T2>
inline typename boost::tuples::element<N, boost::tuple<T1, T2> >::type&
get(std::pair<T1, T2>& pair) {
  return pair_get<N, T1, T2>::get(pair);
}

template <std::size_t N, typename T1, typename T2>
inline const typename boost::tuples::element<N, boost::tuple<T1, T2> >::type&
get(const std::pair<T1, T2>& pair) { 
  return pair_get<N, T1, T2>::get(pair);
}

#endif // end if not C++11 tuple


} // cpp11

namespace cpp0x = cpp11;

#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

// Tool to test whether a type V is among the types of a tuple<...> = T.
template <typename V, typename T>
struct Is_in_tuple;

template <typename V, typename T0, typename... T>
struct Is_in_tuple <V, cpp11::tuple<T0, T...> >
{
  static const bool value = Is_in_tuple<V, cpp11::tuple<T...> >::value;
};

template <typename V, typename... T>
struct Is_in_tuple <V, cpp11::tuple<V, T...> >
{
  static const bool value = true;
};

template <typename V>
struct Is_in_tuple <V, cpp11::tuple<> >
{
  static const bool value = false;
};

#else

// Non-variadic version

template <typename V,typename T>
struct Is_in_tuple;

template <typename V,typename T0,typename T1>
struct Is_in_tuple <V,cpp11::tuple<T0,T1> >
{
  static const bool value = Is_in_tuple<V,cpp11::tuple<T1> >::value;
};

template <typename V, typename T0,typename T1,typename T2>
struct Is_in_tuple <V, cpp11::tuple<T0,T1,T2> >
{
  static const bool value = Is_in_tuple<V,cpp11::tuple<T1,T2> >::value;
};

template <typename V, typename T0,typename T1,typename T2,typename T3>
struct Is_in_tuple <V, cpp11::tuple<T0,T1,T2,T3> >
{
  static const bool value = Is_in_tuple<V,cpp11::tuple<T1,T2,T3> >::value;
};

template <typename V, typename T0,typename T1,typename T2,typename T3,typename T4>
struct Is_in_tuple <V, cpp11::tuple<T0,T1,T2,T3,T4> >
{
  static const bool value = Is_in_tuple<V,cpp11::tuple<T1,T2,T3,T4> >::value;
};

template <typename V, typename T0,typename T1,typename T2,typename T3,typename T4,typename T5>
struct Is_in_tuple <V, cpp11::tuple<T0,T1,T2,T3,T4,T5> >
{
  static const bool value = Is_in_tuple<V,cpp11::tuple<T1,T2,T3,T4,T5> >::value;
};

template <typename V, typename T0,typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
struct Is_in_tuple <V, cpp11::tuple<T0,T1,T2,T3,T4,T5,T6> >
{
  static const bool value = Is_in_tuple<V,cpp11::tuple<T1,T2,T3,T4,T5,T6> >::value;
};


//Conclusions

template <typename V,typename T1>
struct Is_in_tuple <V,cpp11::tuple<T1> >
{
  static const bool value = false;
};

template <typename V>
struct Is_in_tuple <V,cpp11::tuple<V> >
{
  static const bool value = true;
};

template <typename V,typename T1>
struct Is_in_tuple <V,cpp11::tuple<V,T1> >
{
  static const bool value = true;
};

template <typename V,typename T1,typename T2>
struct Is_in_tuple <V,cpp11::tuple<V,T1,T2> >
{
  static const bool value = true;
};

template <typename V,typename T1,typename T2,typename T3>
struct Is_in_tuple <V,cpp11::tuple<V,T1,T2,T3> >
{
  static const bool value = true;
};

template <typename V,typename T1,typename T2,typename T3,typename T4>
struct Is_in_tuple <V,cpp11::tuple<V,T1,T2,T3,T4> >
{
  static const bool value = true;
};

template <typename V,typename T1,typename T2,typename T3,typename T4,typename T5>
struct Is_in_tuple <V,cpp11::tuple<V,T1,T2,T3,T4,T5> >
{
  static const bool value = true;
};

template <typename V,typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
struct Is_in_tuple <V,cpp11::tuple<V,T1,T2,T3,T4,T5,T6> >
{
  static const bool value = true;
};


#endif 

} //namespace CGAL

#endif // CGAL_TUPLE_H
