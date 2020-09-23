// Copyright (c) 2014  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot <sebastien.loriot@cgal.org>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_CONSTRAINED_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_CONSTRAINED_PLACEMENT_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

namespace CGAL {
namespace Surface_mesh_simplification {

template<class BasePlacement, class EdgeIsConstrainedMap>
class Constrained_placement
  : public BasePlacement
{
public:
  Constrained_placement(const EdgeIsConstrainedMap map = EdgeIsConstrainedMap(),
                        const BasePlacement& base = BasePlacement())
    : BasePlacement(base),
      m_ecm(map)
  {}

  template <typename Profile>
  boost::optional<typename Profile::Point> operator()(const Profile& profile) const
  {
    typedef typename Profile::TM                                    TM;
    typedef typename boost::graph_traits<TM>::halfedge_descriptor   halfedge_descriptor;

    for(halfedge_descriptor h : halfedges_around_target(profile.v0(), profile.surface_mesh()))
    {
      if(get(m_ecm, edge(h, profile.surface_mesh())))
        return get(profile.vertex_point_map(), profile.v0());
    }

    for(halfedge_descriptor h : halfedges_around_target(profile.v1(), profile.surface_mesh()))
    {
      if(get(m_ecm, edge(h, profile.surface_mesh())))
        return get(profile.vertex_point_map(), profile.v1());
    }

    return static_cast<const BasePlacement*>(this)->operator()(profile);
  }

private:
  EdgeIsConstrainedMap m_ecm;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_CONSTRAINED_PLACEMENT_H
