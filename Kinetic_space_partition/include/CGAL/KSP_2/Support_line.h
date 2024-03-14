// Copyright (c) 2019 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSP_2_SUPPORT_LINE_H
#define CGAL_KSP_2_SUPPORT_LINE_H

#include <CGAL/license/Kinetic_space_partition.h>

#include <CGAL/KSP/utils.h>
#include <CGAL/KSP_2/Vertex.h>

namespace CGAL {
namespace KSP_2 {
namespace internal {

template <typename GeomTraits>
class Support_line_2
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
  std::vector<std::size_t> m_segments_idx;
  std::vector<std::size_t> m_meta_vertices_idx;
  FT m_minimum;
  FT m_maximum;
  std::size_t m_connected_components;

public:

  Support_line_2() { }

  Support_line_2(const Segment_2& segment)
    : m_minimum((std::numeric_limits<FT>::max)())
    , m_maximum(-(std::numeric_limits<FT>::max)())
    , m_connected_components(1)
  {
    m_origin = CGAL::midpoint(segment.source(), segment.target());
    m_vector = KSP::internal::normalize(Vector_2(segment.source(), segment.target()));
  }

  Line_2 line() const { return Line_2(m_origin, m_vector); }

  const Point_2& origin() const { return m_origin; }
  const Vector_2& vector() const { return m_vector; }

  const FT& minimum() const { return m_minimum; }
  FT& minimum() { return m_minimum; }
  const FT& maximum() const { return m_maximum; }
  FT& maximum() { return m_maximum; }

  const std::size_t& connected_components() const { return m_connected_components; }
  std::size_t& connected_components() { return m_connected_components; }

  CGAL::Bbox_2 bbox() const
  {
    Point_2 pmin = to_2d(m_minimum);
    Point_2 pmax = to_2d(m_maximum);
    return pmin.bbox() + pmax.bbox();
  }

  Segment_2 segment_2() const
  {
    return Segment_2(to_2d(m_minimum), to_2d(m_maximum));
  }

  const std::vector<std::size_t>& segments_idx() const { return m_segments_idx; }
  std::vector<std::size_t>& segments_idx() { return m_segments_idx; }

  const std::vector<std::size_t>& meta_vertices_idx() const { return m_meta_vertices_idx; }
  std::vector<std::size_t>& meta_vertices_idx() { return m_meta_vertices_idx; }

  FT to_1d(const Point_2& point) const
  {
    return m_vector * Vector_2(m_origin, point);
  }

  Point_2 to_2d(const FT& point) const { return m_origin + point * m_vector; }

};

template <typename Kernel>
bool operator== (const Support_line_2<Kernel>& a, const Support_line_2<Kernel>& b)
{
  const typename Kernel::Vector_2& va = a.vector();
  const typename Kernel::Vector_2& vb = b.vector();

  if (CGAL::abs(va * vb) < 0.99999)
    return false;

  return (CGAL::approximate_sqrt(CGAL::squared_distance(b.origin(), a.line())) < 1e-10);
}


#if 0
template <>
bool operator==<CGAL::Exact_predicates_exact_constructions_kernel>
(const Support_line_2<CGAL::Exact_predicates_exact_constructions_kernel>& a,
  const Support_line_2<CGAL::Exact_predicates_exact_constructions_kernel>& b)
{
  return (a.line() == b.line());
}
#endif

} // namespace internal
} // namespace KSP_2
} // namespace CGAL


#endif // CGAL_KSP_2_SUPPORT_LINE_H
