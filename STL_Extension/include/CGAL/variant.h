// Copyright (c) 2023  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sébastien Loriot
//

#ifndef CGAL_VARIANT_H
#define CGAL_VARIANT_H

#include <variant>

namespace CGAL
{

template <class T, class Variant>
struct Is_in_variant;

template <class T, class V1, class ... Vn>
struct Is_in_variant<T, std::variant<V1,Vn...>>
{
  inline static constexpr bool value =
    std::is_same_v<T, V1> || Is_in_variant<T, std::variant<Vn...>>::value;
};

template <class T, class V>
struct Is_in_variant<T, std::variant<V>>
{
  inline static constexpr bool value = std::is_same_v<T, V>;
};

/// equals true iif `T` is a possible type in `Variant` with `Variant` being a `std::variant`
template <class T, class Variant>
inline constexpr bool Is_in_variant_v = Is_in_variant<T, Variant>::value;

// --
template <class T, class Variant>
struct Add_to_variant;

template <class T, class ... Vn>
struct Add_to_variant<T, std::variant<Vn...>>
{
  using type = std::variant<Vn..., T>;
};

/// a `std::variant` with `T` appended to the types of the `std::variant` `Variant`
template< class T, class Variant >
using Add_to_variant_t = typename Add_to_variant<T,Variant>::type;

namespace internal{
template <class Variant, class T1, class ... Tn>
struct Get_variant_impl
{
  using type = typename Get_variant_impl<
    std::conditional_t<Is_in_variant_v<T1, Variant>,
                       Variant,
                       Add_to_variant_t<T1,Variant>>,
    Tn...>::type;
};

template <class Variant, class T1>
struct Get_variant_impl<Variant, T1>
{
  using type = std::conditional_t<Is_in_variant_v<T1, Variant>,
                                  Variant,
                                  Add_to_variant_t<T1,Variant>>;
};
} // end of internal namespace

template <class T1, class ... Tn>
struct Variant_with_no_duplicate
{
  using type = typename internal::Get_variant_impl<std::variant<T1>, Tn ...>::type;
};

/// a `std::variant` with types being all different
template< class ... Tn >
using Variant_with_no_duplicate_t = typename Variant_with_no_duplicate<Tn ... >::type;

} //end of CGAL namespace

#endif
