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
  Pair_optional_adaptor(const std::optional<T>& obj) : std::optional<T>(obj), second(obj.has_value()), first(t_storage.t) {
    if (obj.has_value())
      t_storage.t = *obj;
  } //boost::tuples::detail::swallow_assign

  Pair_optional_adaptor(const std::nullopt_t& obj) : std::optional<T>(std::nullopt), second(false), first(t_storage.t) {}

  Pair_optional_adaptor(std::pair<T, bool>& p) : std::optional<T>(p.second ? p.first : std::optional<T>()), first(t_storage.t), second(p.second), t_storage(p.first) {}

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
    operator std::tuple<T&, std::_Ignore const&>() {
    return std::tuple<T&, std::_Ignore const&>(first, std::ignore);
  }

  CGAL_DEPRECATED_MSG("you are using the deprecated API, please update your code")
    operator std::tuple<std::_Ignore const&, bool&>() {
    return std::tuple<std::_Ignore const&, bool&>(std::ignore, second);
  }

  CGAL_DEPRECATED_MSG("you are using the deprecated API, please update your code")
    operator std::tuple<std::_Ignore const&, std::_Ignore const&>() {
    return std::tuple<std::_Ignore const&, std::_Ignore const&>(std::ignore, std::ignore);
  }
#endif

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
