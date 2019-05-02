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

template <typename Line_3>
class Intersection_line
{
private:

  Line_3 m_line;
  
  KSR::Idx_vector m_support_planes_idx;
  KSR::Idx_vector m_meta_vertices_idx;
  KSR::Idx_vector m_segments_idx;

public:

  Intersection_line () { }

  Intersection_line (const Line_3& line)
    : m_line (line)
  { }

  const Line_3& line() const { return m_line; }

  const KSR::Idx_vector& support_planes_idx() const { return m_support_planes_idx; }
  KSR::Idx_vector& support_planes_idx() { return m_support_planes_idx; }

  const KSR::Idx_vector& meta_vertices_idx() const { return m_meta_vertices_idx; }
  KSR::Idx_vector& meta_vertices_idx() { return m_meta_vertices_idx; }

  const KSR::Idx_vector& segments_idx() const { return m_segments_idx; }
  KSR::Idx_vector& segments_idx() { return m_segments_idx; }


};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_INTERSECTION_LINE_H
