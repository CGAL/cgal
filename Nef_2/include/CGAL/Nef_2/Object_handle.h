// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra
//                 Andreas Fabri
//                 Geert-Jan Giezeman
//                 Michael Seel
//                 Sylvain Pion
// ======================================================================

#ifndef CGAL_OBJECT_HANDLE_H
#define CGAL_OBJECT_HANDLE_H

#include <CGAL/license/Nef_2.h>

#ifdef CGAL_NEF_USE_ANY_OBJECT
#include <CGAL/Object.h>
#include <type_traits>
#include <utility>
#else
#include <variant>
#endif

namespace CGAL {

template <typename... Args>
struct Type_pack {};

#ifdef CGAL_NEF_USE_ANY_OBJECT

template <typename U>
struct Object_handle : Object
{
  using Object::Object;

  template <typename T, typename = std::void_t
    <typename std::iterator_traits<std::decay_t<T>>::value_type>>
  Object_handle(T&& t)
    : Object(std::forward<T>(t), Object::private_tag{})
  {}
};

template <typename T, typename U>
inline bool assign(T& t, const Object_handle<U>& o) { return o.assign(t); }

#else

template <typename U>
struct Object_handle;

template <typename... Args>
struct Object_handle<Type_pack<Args...>>
  : std::variant<std::monostate, Args...> {
  using std::variant<std::monostate, Args...>::variant;
  // needed for compatabiity with CGAL::Object API
  bool empty() const
  { return std::holds_alternative<std::monostate>(*this); }
};

template <typename T, typename U>
inline bool assign(T& t, const Object_handle<U>& oh) {
  const auto* p = std::get_if<T>(&oh);
  if (!p) return false;
  t = *p;
  return true;
}

#endif

} //namespace CGAL

#endif // CGAL_OBJECT_HANDLE_H
