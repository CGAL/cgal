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
#include <CGAL/Dimension.h>
#include <CGAL/tags.h>

namespace CGAL {

template< class Base, class Derived >
struct is_same_or_derived :
  public std::bool_constant<
    ::std::is_same_v< Base, Derived > ||
    ::boost::is_base_and_derived< Base, Derived >::value
  >
{};

template <typename... Ts> using void_t = void;

template <int N> struct Priority_tag : Priority_tag<N+1> {};
template <> struct Priority_tag<8> {};

template <typename...>
constexpr bool always_false_v = false;

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

  template <typename T, typename Dimension> // high priority Priority_tag<0>
  auto
  detect_Bare_point_type(T*, void_t<typename T::Bare_point_3>*, Dimension, Priority_tag<0> = {}) {
    return typename T::Bare_point_3{};
  }

  template <typename T, typename Dimension>
  auto
  detect_Bare_point_type(T*, void_t<typename T::Bare_point>*, Dimension, Priority_tag<1> = {}) {
    return typename T::Bare_point{};
  }

  template <typename T> // low priority Priority_tag<2>: use Tr::Point_3
  auto
  detect_Bare_point_type(T*, void_t<typename T::Point_3>*, Dimension_tag<3>, Priority_tag<2> = {}) {
    return typename T::Point_3{};
  }

  template <typename T> // low priority Priority_tag<2>: use Tr::Point_2 (2D case)
  auto
  detect_Bare_point_type(T*, void_t<typename T::Point_2>*, Dimension_tag<2>, Priority_tag<2> = {}) {
    return typename T::Point_2{};
  }

  template <typename T, typename Dimension> // lowest priority Priority_tag<3>: use Tr::Point, and
                                            // assume it is bare
  auto detect_Bare_point_type(T*, void_t<typename T::Point>*, Dimension, Priority_tag<3> = {}) {
    return typename T::Point{};
  }
} // namespace details

template <typename T, int dim = 3>
struct Bare_point_type {
  using type = decltype(details::detect_Bare_point_type<T>(
      nullptr, nullptr, Dimension_tag<dim>{}, Priority_tag<0>()));
};

template <typename T, int dim = 3>
using Bare_point_type_t = typename Bare_point_type<T, dim>::type;

template <typename Tr>
constexpr bool is_regular_triangulation_v = (!std::is_same_v<Bare_point_type_t<Tr>, typename Tr::Point>);

template <typename Tr, typename = void>
struct Is_periodic_triangulation : Tag_false {};

template <typename Tr>
struct Is_periodic_triangulation<Tr, void_t<typename Tr::Periodic_tag>> : Tr::Periodic_tag {};

template <typename Tr>
constexpr bool is_periodic_triangulation_v = Is_periodic_triangulation<Tr>::value;

template <typename From, typename To>
struct is_convertible_without_narrowing : details::is_convertible_without_narrowing<From, To>
{};

template <typename From, typename To>
inline constexpr bool is_convertible_without_narrowing_v = is_convertible_without_narrowing<From, To>::value;

} // end namespace CGAL

#endif // CGAL_TYPE_TRAITS_H
