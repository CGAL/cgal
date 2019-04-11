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

#ifndef CGAL_KSR_EVENT_H
#define CGAL_KSR_EVENT_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>

namespace CGAL
{

namespace KSR
{

template <typename GeomTraits>
class Event
{
public:
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;

  enum Type
  {
    REGULAR,
    BORDER,
    PARALLEL
  };

private:

  Type m_type;
  std::pair<KSR::size_t, KSR::size_t> m_indices;
  Point_2 m_intersection;

  FT m_time;

public:

  Event() { }

  Event (const Type& type, KSR::size_t index_0, KSR::size_t index_1, const Point_2& intersection, FT time)
    : m_type (type), m_indices(index_0, index_1), m_intersection (intersection), m_time (time)
  {
    // Avoid entering twice the same border event
    if (m_type == BORDER && m_indices.first > m_indices.second)
      std::swap (m_indices.first, m_indices.second);
  }


  const Type& type() const { return m_type; }

  const KSR::size_t& vertex_idx() const { return m_indices.first; }
  KSR::size_t& vertex_idx() { return m_indices.first; }
  const KSR::size_t& segment_idx() const { return m_indices.second; }
  KSR::size_t& segment_idx() { return m_indices.second; }
  const KSR::size_t& other_vertex_idx() const { return m_indices.second; }
  KSR::size_t& other_vertex_idx() { return m_indices.second; }

  bool is_vertex_to_vertex_event() const { return (m_type != REGULAR); }

  const Point_2& intersection() const { return m_intersection; }

  FT time() const { return m_time; }
  
  // Compare two events
  bool operator< (const Event& other) const
  {
    // Use lexicographic comparison of tuples
    return (std::make_tuple (this->m_time, this->m_type, this->m_indices)
            < std::make_tuple (other.m_time, other.m_type, other.m_indices));
  }

  friend std::ostream& operator<< (std::ostream& os, const Event& ev)
  {
    if (ev.m_type == REGULAR)
      os << "Regular ";
    else if (ev.m_type == BORDER)
      os << "Border ";
    else
      os << "Parallel ";
    
    os << "event at t=" << ev.m_time << " between vertex " << ev.m_indices.first << " and "
       << (ev.m_type == REGULAR ? "segment " : "vertex ")
       << ev.m_indices.second << " at point " << ev.m_intersection;
    return os;
  }

};


}} // namespace CGAL::KSR


#endif // CGAL_KSR_EVENT_H
