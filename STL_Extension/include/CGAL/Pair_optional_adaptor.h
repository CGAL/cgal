// Copyright (c) 2024 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau

#ifndef CGAL_PAIR_OPTIONAL_ADAPTOR
#define CGAL_PAIR_OPTIONAL_ADAPTOR

namespace CGAL {

// T is supposed to be a handle
template<typename T>
class Pair_optional_adaptor : public std::optional<T> {
public:
  Pair_optional_adaptor(std::optional<T>& obj) : std::optional<T>(obj), second(obj.has_value()), first(t_storage.t) {
    if (obj.has_value())
      t_storage.t = *obj;
  }

  Pair_optional_adaptor(const std::nullopt_t& obj) : std::optional<T>(std::nullopt), second(false), first(t_storage.t) {}

  Pair_optional_adaptor(std::pair<T, bool>& p) : std::optional<T>(p.second ? p.first : std::optional<T>()), first(t_storage.t), second(p.second), t_storage(p.first) {}

  operator std::pair<T, bool>() {
    return std::pair<T, bool>(first, second);
  }

  T &first;
  bool second;

private:
  union T_value {
    T t;
    int i;
    T_value() : i(0) {}
    T_value(T t) : t(t) {}
  } t_storage;
};

} // CGAL

#endif