// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Baskin Burak Senbaslar
//

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_PLACEMENT_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_core.h>

#include <boost/optional/optional.hpp>

namespace CGAL {
namespace Surface_mesh_simplification {

template<class TM_, typename VCM_>
class GarlandHeckbert_placement
{
public:
  typedef TM_                                                     TM;
  typedef VCM_                                                    Vertex_cost_map;

  typedef typename internal::GarlandHeckbert_core<TM>             GH_core;
  typedef typename GH_core::Matrix4x4                             Matrix4x4;
  typedef typename GH_core::Col4                                  Col4;

  GarlandHeckbert_placement(const Vertex_cost_map& cost_matrices)
    : m_cost_matrices(cost_matrices)
  { }

  template <typename Profile>
  boost::optional<typename Profile::Point> operator()(const Profile& aProfile) const
  {
    CGAL_precondition(get(m_cost_matrices, aProfile.v0()) != Matrix4x4());
    CGAL_precondition(get(m_cost_matrices, aProfile.v1()) != Matrix4x4());

    // the combined matrix has already been computed in the evaluation of the cost...
    const Matrix4x4 combinedMatrix = GH_core::combine_matrices(
                                       get(m_cost_matrices, aProfile.v0()),
                                       get(m_cost_matrices, aProfile.v1()));

    const Col4 p0 = GH_core::point_to_homogenous_column(aProfile.p0());
    const Col4 p1 = GH_core::point_to_homogenous_column(aProfile.p1());
    const Col4 opt = GH_core::optimal_point(combinedMatrix, p0, p1);

    boost::optional<typename Profile::Point> pt = typename Profile::Point(opt(0) / opt(3),
                                                                          opt(1) / opt(3),
                                                                          opt(2) / opt(3));

    return pt;
  }

private:
  const Vertex_cost_map m_cost_matrices;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_GARLANDHECKBERT_PLACEMENT_H
