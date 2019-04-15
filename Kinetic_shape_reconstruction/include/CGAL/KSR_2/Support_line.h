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
  typedef typename Kernel::Segment_2 Segment_2;

private:

  Point_2 m_origin;
  Vector_2 m_vector;
  std::vector<KSR::size_t> m_segments_idx;
  std::vector<KSR::size_t> m_meta_vertices_idx;
  FT m_minimum;
  FT m_maximum;

public:

  Support_line (const Segment_2& segment)
    : m_minimum ((std::numeric_limits<FT>::max)())
    , m_maximum (-(std::numeric_limits<FT>::max)())
  {
    m_origin = CGAL::midpoint (segment.source(), segment.target());
    m_vector = KSR::normalize (Vector_2 (segment.source(), segment.target()));
  }

  Line_2 line() const { return Line_2 (m_origin, m_vector); }

  const Point_2& origin() const { return m_origin; }
  const Vector_2& vector() const { return m_vector; }

  const FT& minimum() const { return m_minimum; }
  FT& minimum() { return m_minimum; }
  const FT& maximum() const { return m_maximum; }
  FT& maximum() { return m_maximum; }

  CGAL::Bbox_2 bbox() const
  {
    Point_2 pmin = to_2d (m_minimum);
    Point_2 pmax = to_2d (m_maximum);
    return pmin.bbox() + pmax.bbox();
  }

  const std::vector<KSR::size_t>& segments_idx() const { return m_segments_idx; }
  std::vector<KSR::size_t>& segments_idx() { return m_segments_idx; }

  const std::vector<KSR::size_t>& meta_vertices_idx() const { return m_meta_vertices_idx; }
  std::vector<KSR::size_t>& meta_vertices_idx() { return m_meta_vertices_idx; }

  FT to_1d (const Point_2& point) const
  {
    return m_vector * Vector_2 (m_origin, point);
  }
  
  Point_2 to_2d (const FT& point) const { return m_origin + point * m_vector; }

};

template <typename Kernel>
bool operator== (const Support_line<Kernel>& a, const Support_line<Kernel>& b)
{
  const typename Kernel::Vector_2& va = a.vector();
  const typename Kernel::Vector_2& vb = b.vector();

  if (CGAL::abs(va * vb) < CGAL_KSR_SAME_VECTOR_TOLERANCE)
    return false;

  return (CGAL::approximate_sqrt(CGAL::squared_distance (b.origin(), a.line())) < CGAL_KSR_SAME_POINT_TOLERANCE);
}


#if 0
template <>
bool operator==<CGAL::Exact_predicates_exact_constructions_kernel>
(const Support_line<CGAL::Exact_predicates_exact_constructions_kernel>& a,
 const Support_line<CGAL::Exact_predicates_exact_constructions_kernel>& b)
{
  return (a.line() == b.line());
}
#endif

}} // namespace CGAL::KSR_2


#endif // CGAL_KSR_2_SUPPORT_LINE_H
