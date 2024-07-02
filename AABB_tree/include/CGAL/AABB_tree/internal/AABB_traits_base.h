// Copyright (c) 2024 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Stéphane Tayeb, Pierre Alliez, Camille Wormser
//

#ifndef CGAL_AABB_TRAITS_BASE_H
#define CGAL_AABB_TRAITS_BASE_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/AABB_tree/internal/Has_nested_type_Shared_data.h>

namespace CGAL {
namespace internal {
namespace AABB_tree {

//helper controlling whether extra data should be stored in the AABB_tree traits class
template <class Primitive, bool has_shared_data = Has_nested_type_Shared_data<Primitive>::value>
struct AABB_traits_base;

template <class Primitive>
struct AABB_traits_base<Primitive, false> {};

template <class Primitive>
struct AABB_traits_base<Primitive, true> {
  typename  Primitive::Shared_data m_primitive_data;

  template <typename ... T>
  void set_shared_data(T&& ... t) {
    m_primitive_data = Primitive::construct_shared_data(std::forward<T>(t)...);
  }
  const typename Primitive::Shared_data& shared_data() const { return m_primitive_data; }
};

}
}
}

#endif
