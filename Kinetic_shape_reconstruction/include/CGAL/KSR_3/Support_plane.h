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

#ifndef CGAL_KSR_3_SUPPORT_PLANE_H
#define CGAL_KSR_3_SUPPORT_PLANE_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>

namespace CGAL
{

namespace KSR_3
{

template <typename GeomTraits>
class Support_plane
{
public:
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Line_3 Line_3;
  typedef typename Kernel::Plane_3 Plane_3;

private:

  Plane_3 m_plane;
  std::vector<KSR::size_t> m_support_lines;

public:

  Support_plane (const Point_3& point, const Vector_3& normal)
    : m_plane (point, normal)
  {
    // Make sure transformations are isomorphic
    Vector_3 base1 = KSR::normalize(m_plane.base1());
    Vector_3 base2 = KSR::normalize(m_plane.base2());
    m_plane = Plane_3 (point, point + base1, point + base2);
  }

  const Plane_3& plane() const { return m_plane; }

  Line_2 to_2d (const Line_3& line) const
  {
    return Line_2 (m_plane.to_2d(line.point()),
                   m_plane.to_2d(line.point() + line.to_vector()));
  }
  
  Point_3 to_3d (const Point_2& point) const { return m_plane.to_3d (point); }
  

  const std::vector<KSR::size_t>& support_lines() const { return m_support_lines; }
  std::vector<KSR::size_t>& support_lines() { return m_support_lines; }
  

};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_SUPPORT_LINE_H
