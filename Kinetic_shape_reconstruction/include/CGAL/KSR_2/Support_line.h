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

#ifndef CGAL_KSR_2_SUPPORT_LINE_H
#define CGAL_KSR_2_SUPPORT_LINE_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>
#include <CGAL/KSR_2/Vertex.h>

namespace CGAL
{

namespace KSR_2
{

template <typename GeomTraits>
class Support_line
{
public:
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Ray_2 Ray_2;
  typedef typename Kernel::Segment_2 Segment_2;

  typedef KSR_2::Vertex<GeomTraits> Vertex;

private:

  Point_2 m_origin;
  Vector_2 m_vector;
  std::vector<KSR::size_t> m_segments;

public:

  Support_line (const Segment_2& segment)
  {
    m_origin = CGAL::midpoint (segment.source(), segment.target());
    m_vector = KSR::normalize (Vector_2 (segment.source(), segment.target()));
  }

  Line_2 line() const { return Line_2 (m_origin, m_vector); }

  const std::vector<KSR::size_t>& segments() const { return m_segments; }
  std::vector<KSR::size_t>& segments() { return m_segments; }

  FT to_1d (const Point_2& point) const
  {
    return m_vector * Vector_2 (m_origin, point);
  }
  
  Point_2 to_2d (const FT& point) const { return m_origin + point * m_vector; }

  Ray_2 to_ray (const Vertex& vertex) const
  {
    return Ray_2 (to_2d(vertex.point()), m_vector * vertex.speed());
  }

};


}} // namespace CGAL::KSR_2


#endif // CGAL_KSR_2_SUPPORT_LINE_H
