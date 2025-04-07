// Copyright (c) 2024  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_STL_EXTENSION_INTERNAL_TUPLE_LIKE_HELPERS_H
#define CGAL_STL_EXTENSION_INTERNAL_TUPLE_LIKE_HELPERS_H

#include <type_traits>
#include <utility>

namespace CGAL::STL_Extension::internal {

  template <typename, typename = void>
  constexpr bool has_tuple_size_v = false;

  template <typename T>
  constexpr bool has_tuple_size_v<T, std::void_t<decltype(std::tuple_size<const T>::value)>> = true;

  template <typename T, bool = has_tuple_size_v<T>>
  constexpr bool tuple_like_of_size_2 = false;

  template <typename T>
  constexpr bool tuple_like_of_size_2<T, true> = (std::tuple_size_v<T> == 2);

} // end namespace CGAL::STL_Extension::internal

#endif // CGAL_STL_EXTENSION_INTERNAL_TUPLE_LIKE_HELPERS_H
