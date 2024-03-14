// Copyright (c) 2019 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSP_2_EVENT_H
#define CGAL_KSP_2_EVENT_H

#include <CGAL/license/Kinetic_space_partition.h>

#include <CGAL/KSP/utils.h>

namespace CGAL {
namespace KSP_2 {
namespace internal {

template <typename GeomTraits>
class Event_queue;

template <typename GeomTraits>
class Event
{
public:
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;

  typedef Event_queue<GeomTraits> Queue;
  friend Queue;

private:

  std::size_t m_vertex_idx;
  std::size_t m_meta_vertex_idx;
  FT m_time;

public:

  Event() { }

  Event(std::size_t vertex_idx, std::size_t meta_vertex_idx, FT time)
    : m_vertex_idx(vertex_idx), m_meta_vertex_idx(meta_vertex_idx), m_time(time)
  { }

  const std::size_t& vertex_idx() const { return m_vertex_idx; }
  std::size_t& vertex_idx() { return m_vertex_idx; }
  const std::size_t& meta_vertex_idx() const { return m_meta_vertex_idx; }
  std::size_t& meta_vertex_idx() { return m_meta_vertex_idx; }

  FT time() const { return m_time; }

  friend std::ostream& operator<< (std::ostream& os, const Event& ev)
  {
    os << "Event at t=" << ev.m_time << " between vertex " << ev.m_vertex_idx
      << " and meta vertex " << ev.m_meta_vertex_idx;
    return os;
  }

};

} // namespace internal
} // namespace KSP_2
} // namespace CGAL


#endif // CGAL_KSP_2_EVENT_H
