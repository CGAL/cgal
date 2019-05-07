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

#ifndef CGAL_KSR_3_SEGMENT_H
#define CGAL_KSR_3_SEGMENT_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

namespace CGAL
{

namespace KSR_3
{

class Segment
{
private:

  KSR::size_t m_intersection_line_idx;
  KSR::size_t m_source_idx;
  KSR::size_t m_target_idx;
  KSR::size_t m_other_source_idx;
  KSR::size_t m_other_target_idx;

public:

  Segment () { }

  Segment (KSR::size_t intersection_line_idx,
           KSR::size_t source_idx = KSR::no_element(),
           KSR::size_t target_idx = KSR::no_element(),
           KSR::size_t other_source_idx = KSR::no_element(),
           KSR::size_t other_target_idx = KSR::no_element())
    : m_intersection_line_idx (intersection_line_idx)
    , m_source_idx (source_idx)
    , m_target_idx (target_idx)
    , m_other_source_idx (other_source_idx)
    , m_other_target_idx (other_target_idx)
  { }

  const KSR::size_t& intersection_line_idx() const { return m_intersection_line_idx; }
  
  const KSR::size_t& source_idx() const { return m_source_idx; }
  KSR::size_t& source_idx() { return m_source_idx; }
  
  const KSR::size_t& target_idx() const { return m_target_idx; }
  KSR::size_t& target_idx() { return m_target_idx; }

  const KSR::size_t& other_source_idx() const { return m_other_source_idx; }
  KSR::size_t& other_source_idx() { return m_other_source_idx; }
  
  const KSR::size_t& other_target_idx() const { return m_other_target_idx; }
  KSR::size_t& other_target_idx() { return m_other_target_idx; }

  KSR::size_t mirror_vertex (KSR::size_t vertex_idx) const
  {
    if (vertex_idx == m_source_idx)
      return m_other_source_idx;
    if (vertex_idx == m_other_source_idx)
      return m_source_idx;
    if (vertex_idx == m_target_idx)
      return m_other_target_idx;
    CGAL_assertion (vertex_idx == m_other_target_idx);
    return m_target_idx;
  }
};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_POLYGON_H
