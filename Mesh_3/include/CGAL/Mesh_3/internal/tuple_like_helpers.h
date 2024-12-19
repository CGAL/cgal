// Copyright (c) 2024  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_3_INTERNAL_TUPLE_LIKE_HELPERS_H
#define CGAL_MESH_3_INTERNAL_TUPLE_LIKE_HELPERS_H

#include <CGAL/license/Mesh_3.h>

#include <type_traits>
#include <tuple>

namespace CGAL::Mesh_3::internal {

  template <typename, typename = void>
  constexpr bool has_tuple_size_v = false;

  template <typename T>
  constexpr bool has_tuple_size_v<T, std::void_t<decltype(std::tuple_size<const T>::value)>> = true;

  template <typename T, bool = has_tuple_size_v<T>>
  constexpr bool tuple_like_of_size_2 = false;

  template <typename T>
  constexpr bool tuple_like_of_size_2<T, true> = (std::tuple_size_v<T> == 2);

} // end namespace CGAL::Mesh_3::internal

#endif // CGAL_MESH_3_INTERNAL_TUPLE_LIKE_HELPERS_H
