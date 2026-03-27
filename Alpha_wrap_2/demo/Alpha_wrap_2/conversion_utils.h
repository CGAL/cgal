// Copyright (c) 2019-2020 X, The Moonshot Factory (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later
//
//
// Author(s)     : Pierre Alliez pierre.alliez@inria.fr
//               : Michael Hemmer mhsaar@gmail.com
//               : Cedric Portaneri cportaneri@gmail.com
//
#ifndef CGAL_ALPHA_WRAP_2_DEMO_CONVERSION_UTILS_H
#define CGAL_ALPHA_WRAP_2_DEMO_CONVERSION_UTILS_H

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace internal {

template <typename Edge>
std::array<decltype(Edge().first->vertex(0)), 2>
delaunay_edge_to_delaunay_vertex_array(const Edge& edge)
{
  return { edge.first->vertex((edge.second + 1)%3),
           edge.first->vertex((edge.second + 2)%3) };
}

template <typename Kernel>
CGAL::Segment_2<Kernel> make_segment_2(const CGAL::Point_2<Kernel>& p0,
                                       const CGAL::Point_2<Kernel>& p1) {
  CGAL_precondition(p0 != p1);
  return CGAL::Segment_2<Kernel>(p0, p1);
}

template <typename Dt2_vertex_handle>
auto delaunay_vertex_array_to_segment(const std::array<Dt2_vertex_handle, 2>& va)
{
  return make_segment_2(va[0]->point(), va[1]->point());
}

template <typename Edge>
auto delaunay_edge_to_segment(const Edge& edge)
{
  return delaunay_vertex_array_to_segment(delaunay_edge_to_delaunay_vertex_array(edge));
}

} // namespace internal
} // namespace Alpha_wraps_2
} // namespace CGAL

#endif  // CGAL_ALPHA_WRAP_2_DEMO_CONVERSION_UTILS_H
