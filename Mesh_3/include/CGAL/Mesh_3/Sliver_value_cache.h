// Copyright (c) 2026 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :

#ifndef CGAL_MESH_3_SLIVER_VALUE_CACHE_H
#define CGAL_MESH_3_SLIVER_VALUE_CACHE_H

#include <CGAL/Time_stamper.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/concurrent_hash_map.h>
#endif

#include <functional>
#include <unordered_map>

namespace CGAL {
namespace Mesh_3 {

// Cell must provide sliver_value(), set_sliver_value(), invalidate_sliver_value(), and sliver_cache_valid()
// FT is the floating point type, Cell_handle is the cell handle type

#ifndef CGAL_MESH_3_NO_SLIVER_VALUE
template <typename Cell_handle, typename FT>
struct Sliver_value_cache_in_cell
{
  bool has_value(const Cell_handle& c) const { return c->is_cache_valid(); }
  FT get(const Cell_handle& c) const { return c->sliver_value(); }
  void set(const Cell_handle& c, FT v) const { c->set_sliver_value(v); }
  void invalidate(const Cell_handle& c) const { c->reset_cache_validity(); }
  void clear() const {}
};
#endif

// Non-concurrent map-based cache
template <typename Cell_handle, typename FT>
struct Sliver_value_cache_map
{
  typedef CGAL::Hash_handles_with_or_without_timestamps Hash_fct;
  mutable std::unordered_map<Cell_handle, FT, Hash_fct> cache;
  bool has_value(const Cell_handle& c) const { return cache.find(c) != cache.end(); }
  FT get(const Cell_handle& c) const { return cache.at(c); }
  void set(const Cell_handle& c, FT v) const { cache[c] = v; }
  void invalidate(const Cell_handle& c) const { cache.erase(c); }
  void clear() const { cache.clear(); }
};

#ifdef CGAL_LINKED_WITH_TBB
// Concurrent map-based cache
template <typename Cell_handle, typename FT>
struct Sliver_value_cache_concurrent_map
{
  typedef CGAL::Hash_handles_with_or_without_timestamps Hash_fct;
  typedef tbb::concurrent_hash_map<Cell_handle, FT, Hash_fct> ConcurrentMap;
  mutable ConcurrentMap cache;
  bool has_value(const Cell_handle& c) const {
    typename ConcurrentMap::const_accessor acc;
    return cache.find(acc, c);
  }
  FT get(const Cell_handle& c) const {
    typename ConcurrentMap::const_accessor acc;
    if (cache.find(acc, c)) return acc->second;
    throw std::out_of_range("Cell_handle not found in concurrent cache");
  }
  void set(const Cell_handle& c, FT v) const {
    typename ConcurrentMap::accessor acc;
    cache.insert(acc, c);
    acc->second = v;
  }
  void invalidate(const Cell_handle& c) const {
    cache.erase(c);
  }
  void clear() const { cache.clear(); }
};
#endif

#ifndef CGAL_MESH_3_NO_SLIVER_VALUE
// Sanity checker cache: combines in-cell and map, checks consistency
template <typename Cell_handle, typename FT>
struct Sliver_value_cache_sanity_checker
{
  Sliver_value_cache_in_cell<Cell_handle, FT> in_cell;
  Sliver_value_cache_map<Cell_handle, FT> map;

  bool has_value(const Cell_handle& c) const {
    bool a = in_cell.has_value(c);
    bool b = map.has_value(c);
    if (a != b) {
      std::cerr << "[Sliver_value_cache_sanity_checker] has_value mismatch: in_cell=" << a << ", map=" << b << std::endl;
      std::abort();
    }
    return a;
  }

  FT get(const Cell_handle& c) const {
    FT a = in_cell.get(c);
    FT b = map.get(c);
    std::cout << "get sliver value of C#" << c->time_stamp() << " = " << a << std::endl;
    if (a != b) {
      std::cerr << "[Sliver_value_cache_sanity_checker] get mismatch: in_cell=" << a << ", map=" << b << std::endl;
      std::abort();
    }
    return a;
  }

  void set(const Cell_handle& c, FT v) const {
    std::cout << "set sliver value of C#" << c->time_stamp() << " to " << v << std::endl;
    in_cell.set(c, v);
    map.set(c, v);
  }

  void invalidate(const Cell_handle& c) const {
    in_cell.invalidate(c);
    map.invalidate(c);
  }

  void clear() const {
    in_cell.clear();
    map.clear();
  }
};
#endif

// Metafunction to select the default sliver cache type for a triangulation
template <typename Tr>
struct Default_sliver_cache
{
#ifdef CGAL_MESH_3_NO_SLIVER_VALUE
#  ifdef CGAL_LINKED_WITH_TBB
  using type = std::conditional_t<
    std::is_same<typename Tr::Concurrency_tag, CGAL::Parallel_tag>::value,
    Sliver_value_cache_concurrent_map<typename Tr::Cell_handle, typename Tr::Geom_traits::FT>,
    Sliver_value_cache_map<typename Tr::Cell_handle, typename Tr::Geom_traits::FT>>;
#  else
  using type = Sliver_value_cache_map<typename Tr::Cell_handle, typename Tr::Geom_traits::FT>;
#  endif
#else
  using type = Sliver_value_cache_sanity_checker<typename Tr::Cell_handle, typename Tr::Geom_traits::FT>;
#endif
};

template <typename Tr>
using Default_sliver_cache_t = typename Default_sliver_cache<Tr>::type;

} // namespace Mesh_3
} // namespace CGAL

#endif // CGAL_MESH_3_SLIVER_VALUE_CACHE_H
