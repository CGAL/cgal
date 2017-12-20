// Copyright (c) 2014  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Sebastien Loriot <sebastien.loriot@cgal.org>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_CONSTRAINED_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_CONSTRAINED_PLACEMENT_H

#include <CGAL/license/Surface_mesh_simplification.h>


#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

namespace CGAL {

namespace Surface_mesh_simplification
{

template<class BasePlacement, class EdgeIsConstrainedMap>
class Constrained_placement : public BasePlacement
{
public:

  EdgeIsConstrainedMap Edge_is_constrained_map;

public:
  Constrained_placement(
    EdgeIsConstrainedMap map=EdgeIsConstrainedMap(),
    BasePlacement base=BasePlacement() )
  : BasePlacement(base)
  , Edge_is_constrained_map(map)
  {}

  template <typename Profile> 
  optional<typename Profile::Point> operator()( Profile const& aProfile ) const
  {
    typedef typename Profile::ECM                                ECM;
    typedef typename CGAL::Halfedge_around_target_iterator<ECM>  in_edge_iterator;

    in_edge_iterator eb, ee ;
    for ( boost::tie(eb,ee) = halfedges_around_target(aProfile.v0(),aProfile.surface_mesh());
      eb != ee ; ++ eb )
    {
      if( get(Edge_is_constrained_map, edge(*eb,aProfile.surface_mesh())) )
        return get(aProfile.vertex_point_map(),
                   aProfile.v0());
    }
    for ( boost::tie(eb,ee) = halfedges_around_target(aProfile.v1(),aProfile.surface_mesh());
      eb != ee ; ++ eb )
    {
      if( get(Edge_is_constrained_map, edge(*eb,aProfile.surface_mesh())) )
        return get(aProfile.vertex_point_map(),
                   aProfile.v1());
    }

    return static_cast<const BasePlacement*>(this)->operator()(aProfile);
  }
};

} // namespace Surface_mesh_simplification

} //namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_CONSTRAINED_PLACEMENT_H
