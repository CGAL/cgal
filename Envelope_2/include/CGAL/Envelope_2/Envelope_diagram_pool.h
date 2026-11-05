// Copyright (c) 2026  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel  <efifogel@gmail.com>

#ifndef CGAL_ENVELOPE_2_ENVELOPE_DIAGRAM_POOL_H
#define CGAL_ENVELOPE_2_ENVELOPE_DIAGRAM_POOL_H

#include <CGAL/license/Envelope_2.h>

#include <vector>
#include <cstddef>

#include <CGAL/assertions.h>

namespace CGAL {
namespace Envelope_2 {
namespace internal {

// obtains the ceiling of the binary logarithm of `n`, which is the maximal
// depth of the recursion of the divide-and-conquer algorithm on `n` curves.
inline std::size_t ceil_log2(std::size_t n) {
  std::size_t depth = 0;
  for (--n; n > 0; n >>= 1) ++depth;
  return depth;
}

} // namespace internal

// ============================================================================
// Primary Template: Default fallback or primary definition
// ============================================================================
template <typename Diagram, bool Pooled = false, bool ValueBased = false>
class Envelope_diagram_pool;

// ============================================================================
// Specialization 1: Non-Pooled Strategy
// Zero data members, zero overhead.
// ============================================================================
template <typename Diagram_>
class Envelope_diagram_pool<Diagram_, false, false> {};

// ============================================================================
// Specialization 2: Pointer-Based Pool Strategy (std::vector<Diagram*>)
// High performance, zero move/copy overhead, stable heap addresses.
// ============================================================================
template <typename Diagram_>
class Envelope_diagram_pool<Diagram_, true, false> {
public:
  using Diagram = Diagram_;

  ~Envelope_diagram_pool() { clear_pool(); }

  void init_pool(std::size_t n) {
    clear_pool();
    if (n > 1) m_diagram_pool.reserve(2 * internal::ceil_log2(n));
  }

  void clear_pool() {
    for (auto* d : m_diagram_pool) delete d;
    m_diagram_pool.clear();
  }

  Diagram* acquire_diagram(Diagram& out_d) {
    if (! m_diagram_pool.empty()) {
      Diagram* d = m_diagram_pool.back();
      m_diagram_pool.pop_back();
      d->clear();
      return d;
    }
    return new Diagram(out_d.shared_geometry_traits_2());
  }

  void release_diagram(Diagram* d) { m_diagram_pool.push_back(d); }

protected:

  std::vector<Diagram*> m_diagram_pool;
};

// ============================================================================
// Specialization 3: Value-Based Pool Strategy (std::vector<Diagram>)
// Contiguous memory, zero raw pointers/new/delete operators.
// ============================================================================
template <typename Diagram_>
class Envelope_diagram_pool<Diagram_, true, true> {
public:
  using Diagram = Diagram_;

  void init_pool(std::size_t n) {
    m_active_count = 0;
    // At most two diagrams per recursion level are active at a time.  The
    // reserved capacity is never exceeded, so the vector never reallocates,
    // and the pointers to the active diagrams remain valid.
    if (n > 1) m_diagram_pool.reserve(2 * internal::ceil_log2(n));
  }

  void clear_pool() { m_diagram_pool.clear(); }

  Diagram* acquire_diagram(Diagram& out_d) {
    if (m_active_count < m_diagram_pool.size()) {
      Diagram& d = m_diagram_pool[m_active_count++];
      d.clear();
      return &d;
    }
    CGAL_assertion(m_diagram_pool.size() < m_diagram_pool.capacity());
    m_diagram_pool.emplace_back(out_d.shared_geometry_traits_2());
    ++m_active_count;
    return &m_diagram_pool.back();
  }

  void release_diagram(Diagram* /* d */) { --m_active_count; }

protected:
  std::vector<Diagram> m_diagram_pool;
  std::size_t m_active_count{0};
};

} // namespace Envelope_2
} // namespace CGAL

#endif
