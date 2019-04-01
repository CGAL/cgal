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

#ifndef CGAL_KSR_3_POLYGON_H
#define CGAL_KSR_3_POLYGON_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

namespace CGAL
{

namespace KSR_3
{

template <typename GeomTraits>
class Polygon
{
public:
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Plane_3 Plane_3;

private:

  std::vector<KSR::size_t> m_vertices;
  KSR::size_t m_support_plane;

public:

  Polygon (KSR::size_t support_plane) : m_support_plane (support_plane) { }

  const std::vector<KSR::size_t>& vertices() const { return m_vertices; }
  std::vector<KSR::size_t>& vertices() { return m_vertices; }
  KSR::size_t support_plane() const { return m_support_plane; }

  void add_vertex (std::size_t idx) { m_vertices.push_back (KSR::size_t(idx)); }

  std::pair<KSR::size_t, KSR::size_t>
  previous_and_next_vertex (KSR::size_t idx)
  {
    std::size_t position = 0;
    while (m_vertices[position] != idx)
      ++ position;

    return std::make_pair (m_vertices[(position + m_vertices.size() - 1) % m_vertices.size()],
                           m_vertices[(position + 1) % m_vertices.size()]);
  }

};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_POLYGON_H
