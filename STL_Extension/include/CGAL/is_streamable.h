// Copyright (c) 2012  GeometryFactory Sarl (France)
// All rights reserved.
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
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
// Author(s)     : Laurent Rineau, Philipp Moeller

#ifndef CGAL_IS_STREAMABLE_H
#define CGAL_IS_STREAMABLE_H

#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_cv.hpp>
#include <boost/static_assert.hpp>
#include <iostream>

namespace CGAL {
namespace internal {
namespace is_streamable 
{
  // A tag type returned by operator == for the any struct in this namespace
  // when T does not support ==.
  struct tag {};

  // This type soaks up any implicit conversions and makes the following operator ==
  // less preferred than any other such operator found via ADL.
  struct any
  {
      // Conversion constructor for any type.
      template <class T>
      any(T const&);
  };

  tag operator<<(any const&, any const&);
  tag operator>>(any const&, any const&);

  // Two overloads to distinguish whether T supports a certain operator expression.
  // The first overload returns a reference to a two-element character array and is chosen if
  // T does not support the expression, such as ==, whereas the second overload returns a char 
  // directly and is chosen if T supports the expression. So using sizeof(check(<expression>))
  // returns 2 for the first overload and 1 for the second overload.
  typedef char yes;
  typedef char (&no)[2];

  no check(tag);

  template <class T>
  yes check(T const&);

  template <class T>
  struct is_streamable_impl
  {
    static typename boost::remove_cv<typename boost::remove_reference<T>::type>::type const & x;
    static typename boost::remove_cv<typename boost::remove_reference<T>::type>::type  & y;
    
    static const bool value = 
      sizeof(is_streamable::check(std::cout << x)) == sizeof(is_streamable::yes) &&
      sizeof(is_streamable::check(std::cin >> y)) == sizeof(is_streamable::yes);
  };

} // end namespace internal::is_streamable
} // end namespace internal


/// is_streamable is a meta-function that checks if a type is streamable
/// 
/// is_streamable<T>::value is true iff the type T has stream operators <<
/// and >>. Otherwise it is false.
template <class T>
struct is_streamable
  : internal::is_streamable::is_streamable_impl<T> {};

} // end namespace CGAL

#endif // CGAL_IS_STREAMABLE_H
