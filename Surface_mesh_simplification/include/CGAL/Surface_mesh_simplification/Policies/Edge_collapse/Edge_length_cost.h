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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_EDGE_LENGHT_COST_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_EDGE_LENGHT_COST_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

namespace CGAL {
namespace Surface_mesh_simplification {

// Edge-length cost: the squared length of the collapsing edge
template<class TM>
class Edge_length_cost
{
public:
  Edge_length_cost() {}

  template <typename Profile, typename T>
  boost::optional<typename Profile::FT> operator()(const Profile& profile, const T& /*placement*/) const
  {
    typedef boost::optional<typename Profile::FT>                     result_type;
    return result_type(profile.geom_traits().compute_squared_distance_3_object()(profile.p0(), profile.p1()));
  }
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_EDGE_LENGHT_COST_H
