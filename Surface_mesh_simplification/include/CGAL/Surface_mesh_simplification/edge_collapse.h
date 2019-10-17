// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/internal/Edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk.h>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

template<class TM,
         class GT,
         class ShouldStop,
         class VertexIndexMap,
         class VertexPointMap,
         class HalfedgeIndexMap,
         class EdgeIsConstrainedMap,
         class GetCost,
         class GetPlacement,
         class Visitor>
int edge_collapse(TM& tmesh,
                  const ShouldStop& should_stop,
                  // optional mesh information policies
                  const GT& traits,
                  const VertexIndexMap& vim, // defaults to get(vertex_index, tmesh)
                  const VertexPointMap& vpm, // defaults to get(vertex_point, tmesh)
                  const HalfedgeIndexMap& him, // defaults to get(edge_index, tmesh)
                  const EdgeIsConstrainedMap& ecm, // defaults to No_constrained_edge_map<TM>()
                  // optional strategy policies - defaults to LindstomTurk
                  const GetCost& get_cost,
                  const GetPlacement& get_placement,
                  Visitor visitor)
{
  typedef EdgeCollapse<TM, GT, ShouldStop,
                       VertexIndexMap, VertexPointMap, HalfedgeIndexMap, EdgeIsConstrainedMap,
                       GetCost, GetPlacement, Visitor> Algorithm;

  Algorithm algorithm(tmesh, traits, should_stop, vim, vpm, him, ecm, get_cost, get_placement, visitor);

  return algorithm.run();
}

struct Dummy_visitor
{
  template<class TM>
  void OnStarted(TM&) const {}
  template<class TM>
  void OnFinished(TM&) const {}
  template<class Profile>
  void OnStopConditionReached(const Profile&) const {}
  template<class Profile, class OFT>
  void OnCollected(const Profile&, const OFT&) const {}
  template<class Profile, class OFT, class Size_type>
  void OnSelected(const Profile&, const OFT&, Size_type, Size_type) const {}
  template<class Profile, class OPoint>
  void OnCollapsing(const Profile&, const OPoint&) const {}
  template<class Profile, class VH>
  void OnCollapsed(const Profile&, VH) const {}
  template<class Profile>
  void OnNonCollapsable(const Profile&) const {}
};

} // namespace internal

template<class TM, class ShouldStop, class NamedParameters>
int edge_collapse(TM& tmesh,
                  const ShouldStop& should_stop,
                  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetGeomTraits<TM, NamedParameters>::type                   Geom_traits;

  return internal::edge_collapse(tmesh, should_stop,
                                 choose_parameter(get_parameter(np, internal_np::geom_traits),
                                                  Geom_traits()),
                                 choose_parameter(get_parameter(np, internal_np::vertex_index),
                                                  get_const_property_map(boost::vertex_index, tmesh)),
                                 choose_parameter(get_parameter(np, internal_np::vertex_point),
                                                  get_property_map(vertex_point, tmesh)),
                                 choose_parameter(get_parameter(np, internal_np::halfedge_index),
                                                  get_const_property_map(boost::halfedge_index, tmesh)),
                                 choose_parameter(get_parameter(np, internal_np::edge_is_constrained),
                                                  No_constrained_edge_map<TM>()),
                                 choose_parameter(get_parameter(np, internal_np::get_cost_policy),
                                                  LindstromTurk_cost<TM>()),
                                 choose_parameter(get_parameter(np, internal_np::get_placement_policy),
                                                  LindstromTurk_placement<TM>()),
                                 choose_parameter(get_parameter(np, internal_np::graph_visitor),
                                                  internal::Dummy_visitor()));
}

template<class TM, class ShouldStop>
int edge_collapse(TM& tmesh, const ShouldStop& should_stop)
{
  return edge_collapse(tmesh, should_stop, CGAL::parameters::all_default());
}

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_H
