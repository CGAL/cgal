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

template <typename GeomTraits>
class Vertex
{
public:
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Ray_2 Ray_2;

private:

  Point_2 m_point;
  Vector_2 m_direction;
  KSR::size_t m_polygon;
  KSR::size_t m_support_plane;
  KSR::size_t m_support_line;

public:

  Vertex (const Point_2& point,
          KSR::size_t polygon = KSR::no_element(),
          KSR::size_t support_plane = KSR::no_element())
    : m_point (point), m_direction (NULL_VECTOR)
    , m_polygon (polygon)
    , m_support_plane (support_plane)
    , m_support_line (KSR::no_element())
  { }

  KSR::size_t polygon() const { return m_polygon; }
  
  KSR::size_t support_plane() const { return m_support_plane; }
  
  const Point_2& point() const { return m_point; }
  Point_2& point() { return m_point; }
  const Vector_2& direction() const { return m_direction; }
  Vector_2& direction() { return m_direction; }

  FT speed() const { return CGAL::approximate_sqrt (m_direction.squared_length()); }

  Ray_2 ray() { return Ray_2 (m_point, m_direction); }

  bool is_frozen() const { return (m_direction == NULL_VECTOR); }
  bool is_constrained() const { return (m_support_line != KSR::no_element()); }

};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_POLYGON_H
