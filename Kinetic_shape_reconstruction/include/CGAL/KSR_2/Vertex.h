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

#ifndef CGAL_KSR_2_VERTEX_H
#define CGAL_KSR_2_VERTEX_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>

namespace CGAL
{

namespace KSR_2
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

  FT m_point;
  FT m_speed;
  KSR::size_t m_segment;
  KSR::size_t m_support_line;
  unsigned int m_remaining_intersections;

public:

  Vertex (FT point,
          KSR::size_t segment = KSR::no_element(),
          KSR::size_t support_line = KSR::no_element())
    : m_point (point), m_speed (0)
    , m_segment (segment)
    , m_support_line (support_line)
    , m_remaining_intersections(0)
  { }

  KSR::size_t segment() const { return m_segment; }
  
  KSR::size_t support_line() const { return m_support_line; }
  
  const FT& point() const { return m_point; }
  FT& point() { return m_point; }
  const FT& speed() const { return m_speed; }
  FT& speed() { return m_speed; }

  const unsigned int& remaining_intersections() const { return m_remaining_intersections; }
  unsigned int& remaining_intersections() { return m_remaining_intersections; }

  bool is_frozen() const { return (m_speed == FT(0)); }

};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_2_POLYGON_H
