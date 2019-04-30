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

  KSR::Idx_vector m_polygons_idx;
  KSR::Idx_vector m_meta_vertices_idx;

public:

  Support_plane () { }

  template <typename PointRange>
  Support_plane (const PointRange& points)
  {
    // Compute support plane
    Vector_3 normal = CGAL::NULL_VECTOR;

    //Newell's method
    for (std::size_t i = 0; i < points.size(); ++ i)
    {
      const Point_3& pa = points[i];
      const Point_3& pb = points[(i+1) % points.size()];
      FT x = normal.x() + (pa.y()-pb.y())*(pa.z()+pb.z());
      FT y = normal.y() + (pa.z()-pb.z())*(pa.x()+pb.x());
      FT z = normal.z() + (pa.x()-pb.x())*(pa.y()+pb.y());
      normal = Vector_3 (x,y,z);
    }
    CGAL_assertion_msg (normal != CGAL::NULL_VECTOR, "Polygon is flat");

    m_plane = Plane_3 (points[0], KSR::normalize(normal));
  }

  Support_plane (const Point_3& point, const Vector_3& normal)
    : m_plane (point, normal)
  {
    // Make sure transformations are isomorphic
    Vector_3 base1 = KSR::normalize(m_plane.base1());
    Vector_3 base2 = KSR::normalize(m_plane.base2());
    m_plane = Plane_3 (point, point + base1, point + base2);
  }

  const Plane_3& plane() const { return m_plane; }

  const KSR::Idx_vector& polygons_idx() const { return m_polygons_idx; }
  KSR::Idx_vector& polygons_idx() { return m_polygons_idx; }

  const KSR::Idx_vector& meta_vertices_idx() const { return m_meta_vertices_idx; }
  KSR::Idx_vector& meta_vertices_idx() { return m_meta_vertices_idx; }

  Point_2 to_2d (const Point_3& point) const
  {
    return m_plane.to_2d (point);
  }

  Line_2 to_2d (const Line_3& line) const
  {
    return Line_2 (m_plane.to_2d(line.point()),
                   m_plane.to_2d(line.point() + line.to_vector()));
  }
  
  Point_3 to_3d (const Point_2& point) const { return m_plane.to_3d (point); }
  

};

template <typename Kernel>
bool operator== (const Support_plane<Kernel>& a, const Support_plane<Kernel>& b)
{
  const typename Kernel::Plane_3& va = a.plane();
  const typename Kernel::Plane_3& vb = b.plane();

  if (CGAL::abs(va.orthogonal_vector() * vb.orthogonal_vector()) < CGAL_KSR_SAME_VECTOR_TOLERANCE)
    return false;

  return (CGAL::approximate_sqrt(CGAL::squared_distance (vb.point(), va)) < CGAL_KSR_SAME_POINT_TOLERANCE);
}


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_SUPPORT_LINE_H
