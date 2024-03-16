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
class Pair_optional_adaptor {
public:
  Pair_optional_adaptor(const std::optional<T>& obj) : first(), second(obj.has_value()) {
    if (obj.has_value())
      first = *obj;
  }

  Pair_optional_adaptor(const std::nullopt_t& obj) : first(), second(false) {}

  Pair_optional_adaptor(std::pair<T, bool>& p) : first(p.first), second(p.second) {}

  operator std::optional<T>() {
    if (second)
      return std::optional<T>(first);
    else
      return std::nullopt;
  }

#ifndef CGAL_NO_DEPRECATED_CODE
  CGAL_DEPRECATED_MSG("you are using the deprecated API, please update your code")
  operator std::pair<T, bool>() {
    return std::pair<T, bool>(first, second);
  }

  CGAL_DEPRECATED_MSG("you are using the deprecated API, please update your code")
  operator boost::tuple<T&, bool&>() {
    return boost::tuple<T&, bool&>(first, second);
  }

  CGAL_DEPRECATED_MSG("you are using the deprecated API, please update your code")
  operator boost::tuple<T&, boost::tuples::detail::swallow_assign>() {
    return boost::tuple<T&, boost::tuples::detail::swallow_assign>(first, boost::tuples::ignore);
  }

  CGAL_DEPRECATED_MSG("you are using the deprecated API, please update your code")
  operator boost::tuple<boost::tuples::detail::swallow_assign, bool&>() {
    return boost::tuple<boost::tuples::detail::swallow_assign, bool&>(boost::tuples::ignore, second);
  }

  CGAL_DEPRECATED_MSG("you are using the deprecated API, please update your code")
  operator boost::tuple<boost::tuples::detail::swallow_assign, boost::tuples::detail::swallow_assign>() {
    return boost::tuple<boost::tuples::detail::swallow_assign, boost::tuples::detail::swallow_assign>(boost::tuples::ignore, boost::tuples::ignore);
  }

  CGAL_DEPRECATED_MSG("you are using the deprecated API, please update your code")
  operator std::tuple<T&, bool&>() {
    return std::tuple<T&, bool&>(first, second);
  }

  CGAL_DEPRECATED_MSG("you are using the deprecated API, please update your code")
  operator std::tuple<T&, decltype(std::ignore)&>() {
    return std::tuple<T&, decltype(std::ignore)&>(first, std::ignore);
  }

  CGAL_DEPRECATED_MSG("you are using the deprecated API, please update your code")
  operator std::tuple<decltype(std::ignore)&, bool&>() {
    return std::tuple<decltype(std::ignore)&, bool&>(std::ignore, second);
  }

  CGAL_DEPRECATED_MSG("you are using the deprecated API, please update your code")
  operator std::tuple<decltype(std::ignore)&, decltype(std::ignore)&>() {
    return std::tuple<decltype(std::ignore)&, decltype(std::ignore)&>(std::ignore, std::ignore);
  }
#endif
  T first;
  bool second;
};

} // CGAL

#endif
