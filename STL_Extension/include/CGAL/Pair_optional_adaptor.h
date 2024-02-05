// Copyright (c) 2024 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
  Pair_optional_adaptor(std::optional<T>& obj) : std::optional<T>(obj), second(obj.has_value()), first(u.t) {
    std::cout << "optional constructor" << std::endl;
    if (obj.has_value())
      u.t = *obj;
  }

  Pair_optional_adaptor(const std::nullopt_t& obj) : std::optional<T>(std::nullopt), second(false), first(u.t) {
    std::cout << "nullopt constructor" << std::endl;
  }

  Pair_optional_adaptor(std::pair<T, bool>& p) : std::optional<T>(b ? p.first : std::nullopt), first(p.first), second(b) {
    std::cout << "pair constructor" << std::endl;
  }

  operator std::pair<T, bool>() {
    return std::pair<T, bool>(first, second);
  }

  operator std::optional<T>() {
    if (second)
      return std::optional<T>(first);
    else return std::nullopt;
  }

  T &first;
  bool second;

private:
  union U {
    T t;
    int i;
    U() : i(0) {}
    U(T t) : t(t) {}
  } u;
};

} // CGAL

#endif