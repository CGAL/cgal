// Copyright (c) 2007  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Meyer

#ifndef CGAL_TYPE_TRAITS_H
#define CGAL_TYPE_TRAITS_H

#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/mpl/or.hpp>

#include <type_traits>

namespace CGAL {

template< class Base, class Derived >
struct is_same_or_derived :
  public std::bool_constant<
    ::std::is_same_v< Base, Derived > ||
    ::boost::is_base_and_derived< Base, Derived >::value
  >
{};

template<typename T>
struct is_nothrow_movable :
  public std::bool_constant<
    std::is_nothrow_move_constructible_v<T> &&
    std::is_nothrow_move_assignable_v<T>
  >
{};

template<typename T>
inline constexpr bool is_nothrow_movable_v = is_nothrow_movable<T>::value;

namespace cpp20 {

  template<class T>
  struct type_identity { using type = T; };

  template<class T>
  using type_identity_t = typename type_identity<T>::type;

  template< class T >
  struct remove_cvref {
      typedef std::remove_cv_t<std::remove_reference_t<T>> type;
  };

  template< class T >
  using remove_cvref_t = typename remove_cvref<T>::type;

} // end namespace cpp20

namespace details {
  template <typename From, typename To, typename = void>
  struct is_convertible_without_narrowing : std::false_type
  {};

  template <typename From, typename To>
  struct is_convertible_without_narrowing<From,
                                          To,
                                          std::void_t<decltype(cpp20::type_identity_t<To[]>{std::declval<From>()})>>
      : std::is_convertible<From, To>
  {};
}

template <typename From, typename To>
struct is_convertible_without_narrowing : details::is_convertible_without_narrowing<From, To>
{};

template <typename From, typename To>
inline constexpr bool is_convertible_without_narrowing_v = is_convertible_without_narrowing<From, To>::value;

namespace is_complete_internals
{
    template<class T>
    std::enable_if_t<sizeof(T) != 0, std::true_type>
    check(T(*)());
    
    std::false_type check(...);
};

template<class T, class Base = decltype(is_complete_internals::check(typename std::enable_if<true, T(*)()>::type()))>
struct is_complete : Base { };

namespace is_complete_testsuite
{

  struct S;
  static_assert(!is_complete<S>::value, "error");
  struct S
  {
    static_assert(!is_complete<S>::value, "error");
  };
  static_assert(is_complete<S>::value, "error");
}

} // end namespace CGAL

#endif // CGAL_TYPE_TRAITS_H
