// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_FACE_COUNT_RATIO_STOP_PREDICATE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_FACE_COUNT_RATIO_STOP_PREDICATE_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/boost/graph/internal/helpers.h>
#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

namespace CGAL {
namespace Surface_mesh_simplification {

// Stops when the ratio of initial to current number of faces is below some value.
template<class TM_>
class Face_count_ratio_stop_predicate
{
public:
  typedef TM_                                                 TM;
  typedef typename boost::graph_traits<TM>::edges_size_type   size_type;

  Face_count_ratio_stop_predicate(const double ratio,
                                  const TM& tmesh)
    : m_ratio(ratio), m_initial_face_count(CGAL::internal::exact_num_faces(tmesh))
  {
    CGAL_warning(0. < ratio && ratio <= 1.);
  }

  template <typename F, typename Profile>
  bool operator()(const F& /*current_cost*/,
                  const Profile& profile,
                  size_type /*initial_edge_count*/,
                  size_type /*current_edge_count*/) const
  {
    const std::size_t current_face_count = CGAL::internal::exact_num_faces(profile.surface_mesh());
    return (static_cast<double>(current_face_count) / static_cast<double>(m_initial_face_count)) < m_ratio;
  }

private:
  const double m_ratio;
  const std::size_t m_initial_face_count;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_FACE_COUNT_RATIO_STOP_PREDICATE_H
