// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_ARRAY_H
#define CGAL_ARRAY_H

#include <CGAL/config.h>
#include <array>
#include <utility>

namespace CGAL {

// The make_array() function simply constructs an std::array.
// It is needed for cases where a std::array is used as a class data
// member and you want to initialize it in the member initializers list.
// It is also optimized: no spurious copies of the objects are made,
// provided the compiler does the NRVO.  So this is better than
// default construction followed by assignment.

// I proposed it for Boost.Array, but it has not been integrated so far.
// See the thread starting at
// https://lists.boost.org/Archives/boost/2006/08/109003.php
//
// C++0x has it under discussion here :
// https://www.open-std.org/jtc1/sc22/wg21/docs/lwg-active.html#851

// Hopefully C++0x will fix this properly with initializer_lists.
// So, it's temporary, therefore I do not document it and keep it internal.

// NOTE : The above is actually untrue !  It is possible to do :
//     struct S2 {
//       typedef std::array<M,2> Ar;
//       Ar m;
//       S2 (const M&a) : m ((Ar) { { a, a } }) {}
//     };
// without spurious copies...  Except VC++ does not eat it :-(

// It's also untrue that this is not documented...  It is !

template <typename T, typename ...Args>
struct Make_array_element_type {
  using type = T;
};

template <typename ...Args>
struct Make_array_element_type<void, Args...> {
  using type = typename std::common_type_t<Args...>;
};

template <typename T, typename ...Args>
using Make_array_element_type_t = typename Make_array_element_type<T, Args...>::type;

template<typename T = void, typename... Args>
constexpr
std::array<Make_array_element_type_t<T, Args...>, sizeof...(Args) >
make_array(Args&& ... args)
{
  return {{ std::forward<Args>(args)... }};
}

// Functor version
struct Construct_array
{
  template <typename... Args>
  constexpr
  decltype(auto)
  operator()(Args&& ... args) const
  {
    return make_array( std::forward<Args>(args)... );
  }
};

template <std::size_t...Is, typename T>
constexpr std::array<T, sizeof...(Is)>
make_filled_array_aux(const T& value, std::index_sequence<Is...>)
{
  return {(static_cast<void>(Is), value)...};
}

template <std::size_t N, typename T>
constexpr std::array<T, N> make_filled_array(const T& value)
{
  return make_filled_array_aux(value, std::make_index_sequence<N>());
}

} //namespace CGAL

#endif // CGAL_ARRAY_H
