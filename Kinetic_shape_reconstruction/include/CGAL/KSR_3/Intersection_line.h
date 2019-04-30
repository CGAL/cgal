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

#ifndef CGAL_KSR_3_INTERSECTION_LINE_H
#define CGAL_KSR_3_INTERSECTION_LINE_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>

namespace CGAL
{

namespace KSR_3
{

template <typename Line_2>
class Intersection_line
{
private:

  Line_2 m_line;
  KSR::size_t m_support_plane_idx;
  KSR::size_t m_meta_line_idx;
  KSR::Idx_vector m_vertices_idx;

public:

  Intersection_line () { }

  Intersection_line (const Line_2& line, KSR::size_t support_plane_idx)
    : m_line (line), m_support_plane_idx (support_plane_idx), m_meta_line_idx (KSR::no_element())
  { }

  const Line_2& line() const { return m_line; }

  const KSR::size_t& support_plane_idx() const { return m_support_plane_idx; }

  const KSR::size_t& meta_line_idx() const { return m_meta_line_idx; }
  KSR::size_t& meta_line_idx() { return m_meta_line_idx; }

  const KSR::Idx_vector& vertices_idx() const { return m_vertices_idx; }
  KSR::Idx_vector& vertices_idx() { return m_vertices_idx; }

};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_INTERSECTION_LINE_H
