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

private:

  KSR::size_t m_vertex_idx;
  KSR::size_t m_intersection_line_idx;

  FT m_time;

public:

  Event(KSR::size_t vertex_idx, KSR::size_t intersection_line_idx, FT time)
    : m_vertex_idx (vertex_idx), m_intersection_line_idx (intersection_line_idx), m_time (time)
  { }

  KSR::size_t vertex_idx() const { return m_vertex_idx; }
  KSR::size_t intersection_line_idx() const { return m_intersection_line_idx; }
  FT time() const { return m_time; }

  // Compare two events
  bool operator< (const Event& other) const
  {
    if (this->m_time == other.m_time)
    {
      if (this->m_vertex_idx == other.m_vertex_idx)
        return this->m_intersection_line_idx < other.m_intersection_line_idx;
      return this->m_vertex_idx < other.m_vertex_idx;
    }
    return this->m_time < other.m_time;
  }

  friend std::ostream& operator<< (std::ostream& os, const Event& ev)
  {
    os << "Event at t=" << ev.m_time << " between vertex " << ev.m_vertex_idx << " and line " << ev.m_intersection_line_idx;
    return os;
  }

};


}} // namespace CGAL::KSR


#endif // CGAL_KSR_EVENT_H
