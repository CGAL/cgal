// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel  <efifogel@gmail.com>

#ifndef CGAL_ENVELOPE_2_ENVELOPE_DIAGRAM_POOL_H
#define CGAL_ENVELOPE_2_ENVELOPE_DIAGRAM_POOL_H

#include <vector>
#include <cstddef>

namespace CGAL {
namespace Envelope_2 {

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
    if (n > 1) {
      const std::size_t max_depth = 64 - __builtin_clzll(n - 1);
      m_diagram_pool.reserve(2 * max_depth);
    }
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
    return new Diagram(out_d.shared_geometry_traits_1());
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
    if (n > 1) {
      const std::size_t max_depth = 64 - __builtin_clzll(n - 1);
      m_diagram_pool.reserve(2 * max_depth);
    }
  }

  void clear_pool() { m_diagram_pool.clear(); }

  Diagram* acquire_diagram(Diagram& out_d) {
    if (m_active_count < m_diagram_pool.size()) {
      Diagram& d = m_diagram_pool[m_active_count++];
      d.clear();
      return &d;
    }
    m_diagram_pool.emplace_back(out_d.shared_geometry_traits_1());
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
