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

#ifndef CGAL_KSP_2_VERTEX_H
#define CGAL_KSP_2_VERTEX_H

#include <CGAL/license/Kinetic_space_partition.h>

#include <CGAL/KSP/utils.h>

namespace CGAL {
namespace KSP_2 {
namespace internal {

template <typename FT>
class Vertex
{
private:

  FT m_point;
  FT m_direction;
  std::size_t m_segment_idx;
  unsigned int m_remaining_intersections;
  std::size_t m_meta_vertex_idx;

public:

  Vertex() { }

  Vertex(FT point,
    std::size_t segment_idx = std::size_t(-1),
    unsigned int remaining_intersections = 0)
    : m_point(point)
    , m_direction(0)
    , m_segment_idx(segment_idx)
    , m_remaining_intersections(remaining_intersections)
    , m_meta_vertex_idx(std::size_t(-1))
  {
  }

  const std::size_t& segment_idx() const { return m_segment_idx; }
  std::size_t& segment_idx() { return m_segment_idx; }

  FT point(FT time) const { return m_point + time * m_direction; }
  const FT& direction() const { return m_direction; }
  FT& direction() { return m_direction; }
  FT speed() const { return CGAL::abs(m_direction); }

  const unsigned int& remaining_intersections() const { return m_remaining_intersections; }
  unsigned int& remaining_intersections() { return m_remaining_intersections; }

  const std::size_t& meta_vertex_idx() const { return m_meta_vertex_idx; }
  std::size_t& meta_vertex_idx() { return m_meta_vertex_idx; }

  bool is_frozen() const { return (m_direction == FT(0)); }
  void freeze(FT time)
  {
    m_point = m_point + time * m_direction;
    m_direction = FT(0);
    m_remaining_intersections = 0;
  }

  friend std::ostream& operator<< (std::ostream& os, const Vertex& vertex)
  {
    os << "vertex(" << vertex.m_point << "," << vertex.m_direction << ") on segment " << vertex.m_segment_idx
      << " and meta vertex " << vertex.meta_vertex_idx() << " with "
      << vertex.m_remaining_intersections << " remaining intersection(s)";
    return os;
  }

};

} // namespace internal
} // namespace KSP_2
} // namespace CGAL

#endif // CGAL_KSP_2_VERTEX_H
