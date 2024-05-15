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

#ifndef CGAL_KSP_2_META_VERTEX_H
#define CGAL_KSP_2_META_VERTEX_H

#include <CGAL/license/Kinetic_space_partition.h>

#include <CGAL/KSP/utils.h>
#include <set>

namespace CGAL {
namespace KSP_2 {
namespace internal {

template <typename Point_2>
class Meta_vertex
{
private:

  Point_2 m_point;
  std::set<std::size_t> m_support_lines_idx;

  std::set<std::size_t> m_deadends;

public:

  Meta_vertex() { }

  Meta_vertex(const Point_2& point) : m_point(point) { }

  const Point_2& point() const { return m_point; }

  const std::set<std::size_t>& support_lines_idx() const { return m_support_lines_idx; }
  std::set<std::size_t>& support_lines_idx() { return m_support_lines_idx; }

  void make_deadend_of(std::size_t support_line_idx)
  {
    m_deadends.insert(support_line_idx);
  }

  bool is_deadend_of(std::size_t support_line_idx) const
  {
    return m_deadends.find(support_line_idx) != m_deadends.end();
  }

  void make_no_longer_deadend_of(std::size_t support_line_idx)
  {
    m_deadends.erase(support_line_idx);
  }


};

} // namespace internal
} // namespace KSP_2
} // namespace CGAL


#endif // CGAL_KSP_2_META_VERTEX_H
