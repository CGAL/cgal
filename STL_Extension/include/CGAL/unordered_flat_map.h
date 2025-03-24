// Copyright (c) 2025 GeometryFactory Sarl (France).
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_UNORDERED_FLAT_MAP_H
#define CGAL_UNORDERED_FLAT_MAP_H

#include <CGAL/config.h>

#include <boost/version.hpp>
#if BOOST_VERSION >= 108100 && !defined(CGAL_USE_BOOST_UNORDERED)
#  define CGAL_USE_BOOST_UNORDERED 1
#endif

#if CGAL_USE_BARE_STD_SET
#  define CGAL_USE_BARE_STD_MAP CGAL_USE_BARE_STD_SET
#endif

#if CGAL_USE_BARE_STD_MAP // to benchmark with the ordered std::map
#  include <map>
#  include <set>
#elif CGAL_USE_BOOST_UNORDERED
#  include <boost/unordered/unordered_flat_map.hpp>
#  include <boost/unordered/unordered_flat_set.hpp>
#else // Boost before 1.81.0, use the C++11 std::unordered_map
#  include <unordered_map>
#  include <unordered_set>
#endif

#include <functional> // for std::hash, std::equal_to
#include <memory> // for std::allocator

namespace CGAL {

template <
  typename Key,
  typename T,
  typename Hash = std::hash<Key>,
  typename KeyEqual = std::equal_to<Key>,
  typename Allocator = std::allocator<std::pair<const Key, T>>
  >
#if CGAL_USE_BARE_STD_MAP

  using unordered_flat_map = std::map<Key, T, std::less<Key>, Allocator>;

#elif CGAL_USE_BOOST_UNORDERED

  using unordered_flat_map = boost::unordered_flat_map<Key, T, Hash, KeyEqual, Allocator>;

#else // use the C++11 std::unordered_map

  using unordered_flat_map = std::unordered_map<Key, T, Hash, KeyEqual, Allocator>;

#endif

template <
  typename Key,
  typename Hash = std::hash<Key>,
  typename KeyEqual = std::equal_to<Key>,
  typename Allocator = std::allocator<Key>
  >
#if CGAL_USE_BARE_STD_MAP

  using unordered_flat_set = std::set<Key, std::less<Key>, Allocator>;

#elif CGAL_USE_BOOST_UNORDERED

  using unordered_flat_set = boost::unordered_flat_set<Key, Hash, KeyEqual, Allocator>;

#else // use the C++11 std::unordered_set

  using unordered_flat_set = std::unordered_set<Key, Hash, KeyEqual, Allocator>;

#endif

} // end namespace CGAL

#endif // CGAL_UNORDERED_FLAT_MAP_H
