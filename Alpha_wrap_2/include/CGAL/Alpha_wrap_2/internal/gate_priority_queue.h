// Copyright (c) 2019-2022 Google LLC (USA).
// Copyright (c) 2025 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez
//                 Cedric Portaneri,
//                 Mael Rouxel-Labbé
//                 Andreas Fabri
//                 Michael Hemmer
//
#ifndef CGAL_ALPHA_WRAP_2_INTERNAL_GATE_PRIORITY_QUEUE_H
#define CGAL_ALPHA_WRAP_2_INTERNAL_GATE_PRIORITY_QUEUE_H

#include <CGAL/license/Alpha_wrap_2.h>

#include <boost/property_map/property_map.hpp>

#include <cassert>
#include <iostream>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace internal {

#ifdef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE

// Represents an alpha-traversable edge in the mutable priority queue
template <typename Tr>
class Gate
{
  using Edge = typename Tr::Edge;
  using FT = typename Tr::Geom_traits::FT;

private:
  Edge m_edge;
  FT m_priority; // circumcircle sq_radius
  bool m_is_permissive_edge;

public:
  // Constructors
  Gate(const Edge& edge,
       const FT& priority,
       const bool is_permissive_edge)
    :
      m_edge(edge),
      m_priority(priority),
      m_is_permissive_edge(is_permissive_edge)
  {
    CGAL_assertion(priority >= 0);
  }

  // This overload is only used for contains() and erase(), priority and bbox flag are dummy value
  Gate(const Edge& edge)
    : Gate(edge, 0, false)
  { }

public:
  const Edge& edge() const { return m_edge; }
  const FT& priority() const { return m_priority; }
  bool is_permissive_edge() const { return m_is_permissive_edge; }
};

struct Less_gate
{
  template <typename Tr>
  bool operator()(const Gate<Tr>& a, const Gate<Tr>& b) const
  {
    // If one is permissive and the other is not, give priority to the permissive edge.
    //
    // The permissive edge are given highest priority because they need to be treated
    // regardless of their circumradius. Treating them first allow the part that depends
    // on alpha to be treated uniformly in a way: whatever the alpha, all permissive faces
    // will first be treated.
    if(a.is_permissive_edge() != b.is_permissive_edge())
      return a.is_permissive_edge();

    if(a.priority() == b.priority())
    {
      // arbitrary, the sole purpose is to make it a total order for determinism
      if(a.edge().first->time_stamp() == b.edge().first->time_stamp())
        return a.edge().second < b.edge().second;

      return a.edge().first->time_stamp() < b.edge().first->time_stamp();
    }

    return a.priority() > b.priority();
  }
};

#else // CGAL_AW2_USE_SORTED_PRIORITY_QUEUE

// Represents an alpha-traversable edge in the mutable priority queue
template <typename Tr>
class Gate
{
  using Edge = typename Tr::Edge;
  using FT = typename Tr::Geom_traits::FT;

private:
  Edge m_edge, m_mirror_edge;
  const unsigned int m_erase_counter_mem;
  const unsigned int m_mirror_erase_counter_mem;

public:
  // Constructors
  Gate(const Edge& edge,
       const Tr& tr)
    :
      m_edge(edge),
      m_mirror_edge(tr.mirror_edge(edge)),
      m_erase_counter_mem(m_edge.first->erase_counter()),
      m_mirror_erase_counter_mem(m_mirror_edge.first->erase_counter())
  {
  }

public:
  const Edge& edge() const { return m_edge; }

  bool is_zombie() const
  {
    return (m_edge.first->erase_counter() != m_erase_counter_mem) ||
           (m_mirror_edge.first->erase_counter() != m_mirror_erase_counter_mem);
  }
};

#endif // CGAL_AW3_USE_SORTED_PRIORITY_QUEUE

template <typename Tr>
struct Gate_ID_PM
{
  using key_type = Gate<Tr>;
  using value_type = std::size_t;
  using reference = std::size_t;
  using category = boost::readable_property_map_tag;

  inline friend value_type get(Gate_ID_PM, const key_type& k)
  {
    using Edge = typename Tr::Edge;

    const Edge& f = k.edge();
    return (3 * f.first->time_stamp() + f.second);
  }
};

} // namespace internal
} // namespace Alpha_wraps_2
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_2_INTERNAL_GATE_PRIORITY_QUEUE_H
