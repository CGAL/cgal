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

#ifndef CGAL_KSR_2_META_VERTEX_H
#define CGAL_KSR_2_META_VERTEX_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>

namespace CGAL
{

namespace KSR_2
{

template <typename GeomTraits>
class Meta_vertex
{
public:
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Ray_2 Ray_2;

private:

  Point_2 m_point;
  std::vector<KSR::size_t> m_vertices_idx;

public:

  Meta_vertex (const Point_2& point) : m_point (point) { }

  const Point_2& point() const { return m_point; }

  const std::vector<KSR::size_t>& vertices_idx() const { return m_vertices_idx; }
  std::vector<KSR::size_t>& vertices_idx() { return m_vertices_idx; }

};


}} // namespace CGAL::KSR_2


#endif // CGAL_KSR_2_META_VERTEX_H
