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

#ifndef CGAL_KSR_2_EVENT_H
#define CGAL_KSR_2_EVENT_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>

namespace CGAL
{

namespace KSR_2
{

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

  KSR::size_t m_vertex_idx;
  KSR::size_t m_meta_vertex_idx;
  FT m_time;

public:

  Event () { }

  Event (KSR::size_t vertex_idx, KSR::size_t meta_vertex_idx, FT time)
    : m_vertex_idx (vertex_idx), m_meta_vertex_idx (meta_vertex_idx), m_time (time)
  { }

  const KSR::size_t& vertex_idx() const { return m_vertex_idx; }
  KSR::size_t& vertex_idx() { return m_vertex_idx; }
  const KSR::size_t& meta_vertex_idx() const { return m_meta_vertex_idx; }
  KSR::size_t& meta_vertex_idx() { return m_meta_vertex_idx; }

  FT time() const { return m_time; }
  
  friend std::ostream& operator<< (std::ostream& os, const Event& ev)
  {
    os << "Event at t=" << ev.m_time << " between vertex " << ev.m_vertex_idx
       << " and meta vertex " << ev.m_meta_vertex_idx;
    return os;
  }

};


}} // namespace CGAL::KSR_2


#endif // CGAL_KSR_2_EVENT_H
