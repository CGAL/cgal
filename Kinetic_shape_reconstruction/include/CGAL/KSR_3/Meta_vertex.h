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

#ifndef CGAL_KSR_3_META_VERTEX_H
#define CGAL_KSR_3_META_VERTEX_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>
#include <set>

namespace CGAL
{

namespace KSR_3
{

template <typename Point_3>
class Meta_vertex
{
private:

  Point_3 m_point;
  
  std::set<KSR::size_t> m_support_planes_idx;

public:

  Meta_vertex () { }

  Meta_vertex (const Point_3& point) : m_point (point) { }

  const Point_3& point() const { return m_point; }

  const std::set<KSR::size_t>& support_planes_idx() const { return m_support_planes_idx; }
  std::set<KSR::size_t>& support_planes_idx() { return m_support_planes_idx; }

};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_META_VERTEX_H
