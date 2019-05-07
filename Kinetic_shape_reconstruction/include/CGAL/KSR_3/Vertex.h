// Copyright (c) 2019 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_3_VERTEX_H
#define CGAL_KSR_3_VERTEX_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>

namespace CGAL
{

namespace KSR_3
{

template <typename Kernel>
class Vertex
{
private:

  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;

  Point_2 m_point;
  Vector_2 m_direction;
  KSR::size_t m_polygon_idx;
  unsigned int m_remaining_intersections;
  KSR::size_t m_meta_vertex_idx;
  KSR::size_t m_segment_idx;

public:

  Vertex () { }

  Vertex (const Point_2& point,
          KSR::size_t polygon_idx,
          unsigned int remaining_intersections = 0)
    : m_point (point)
    , m_direction (CGAL::NULL_VECTOR)
    , m_polygon_idx (polygon_idx)
    , m_remaining_intersections(remaining_intersections)
    , m_meta_vertex_idx (KSR::no_element())
    , m_segment_idx (KSR::no_element())
  {
  }

  const KSR::size_t& segment_idx() const { return m_segment_idx; }
  KSR::size_t& segment_idx() { return m_segment_idx; }
  
  const KSR::size_t& polygon_idx() const { return m_polygon_idx; }
  KSR::size_t& polygon_idx() { return m_polygon_idx; }

  Point_2 point (FT time) const { return m_point + time * m_direction; }
  const Vector_2& direction() const { return m_direction; }
  Vector_2& direction() { return m_direction; }
  FT speed() const { return CGAL::approximate_sqrt (m_direction * m_direction); }

  const unsigned int& remaining_intersections() const { return m_remaining_intersections; }
  unsigned int& remaining_intersections() { return m_remaining_intersections; }

  const KSR::size_t& meta_vertex_idx() const { return m_meta_vertex_idx; }
  KSR::size_t& meta_vertex_idx() { return m_meta_vertex_idx; }
  
  bool is_frozen() const { return (m_direction == CGAL::NULL_VECTOR); }
  void freeze(FT time)
  {
    m_point = m_point + time * m_direction;
    m_direction = CGAL::NULL_VECTOR;
    m_remaining_intersections = 0;
  }

  friend std::ostream& operator<< (std::ostream& os, const Vertex& vertex)
  {
    os << "vertex(" << vertex.m_point << "," << vertex.m_direction << ") on support plane " << vertex.m_support_plane_idx
       << " and meta vertex " << vertex.meta_vertex_idx() << " with "
       << vertex.m_remaining_intersections << " remaining intersection(s)";
    return os;
  }

};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_VERTEX_H
