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
#define CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_H 1

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/named_function_params.h>

#include <CGAL/Surface_mesh_simplification/Detail/Edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk.h>

namespace CGAL {
namespace Surface_mesh_simplification {

template<class TM,
         class ShouldStop,
         class VertexIndexMap,
         class VertexPointMap,
         class EdgeIndexMap,
         class EdgeIsConstrainedMap,
         class GetCost,
         class GetPlacement,
         class Visitor>
int edge_collapse(TM& aSurface,
                  const ShouldStop& aShould_stop,
                  // optional mesh information policies
                  const VertexIndexMap& aVertex_index_map, // defaults to get(vertex_index, aSurface)
                  const VertexPointMap& aVertex_point_map, // defaults to get(vertex_point, aSurface)
                  const EdgeIndexMap& aEdge_index_map, // defaults to get(edge_index, aSurface)
                  const EdgeIsConstrainedMap& aEdge_is_constrained_map, // defaults to No_constrained_edge_map<TM>()
                  // optional strategy policies - defaults to LindstomTurk
                  const GetCost& aGet_cost,
                  const GetPlacement& aGet_placement,
                  Visitor aVisitor)
{
  typedef EdgeCollapse<TM, ShouldStop,
                       VertexIndexMap, VertexPointMap, EdgeIndexMap, EdgeIsConstrainedMap,
                       GetCost, GetPlacement, Visitor> Algorithm;

  Algorithm algorithm(aSurface, aShould_stop,
                      aVertex_index_map, aVertex_point_map, aEdge_index_map, aEdge_is_constrained_map,
                      aGet_cost, aGet_placement, aVisitor);

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
  void OnCollected(Profile const&, const OFT&) const {}
  template<class Profile, class OFT, class Size_type>
  void OnSelected(Profile const&, const OFT&, Size_type, Size_type) const {}
  template<class Profile, class OPoint>
  void OnCollapsing(const Profile&, const OPoint&) const {}
  template<class Profile, class VH>
  void OnCollapsed(const Profile&, VH) const {}
  template<class Profile>
  void OnNonCollapsable(Profile const&) const {}
};

template<class TM, class ShouldStop, class P, class T, class R>
int edge_collapse(TM& aSurface,
                  const ShouldStop& aShould_stop,
                  cgal_bgl_named_params<P,T,R> const& aParams)
{
  using boost::choose_param;
  using boost::choose_const_pmap;
  using boost::get_param;

  LindstromTurk_params lPolicyParams;
  internal_np::graph_visitor_t vis = internal_np::graph_visitor_t();

  return edge_collapse(aSurface, aShould_stop,
                       choose_const_pmap(get_param(aParams,internal_np::vertex_index),aSurface,boost::vertex_index),
                       choose_pmap(get_param(aParams,internal_np::vertex_point),aSurface,boost::vertex_point),
                       choose_const_pmap(get_param(aParams,internal_np::halfedge_index),aSurface,boost::halfedge_index),
                       choose_param(get_param(aParams,internal_np::edge_is_constrained),No_constrained_edge_map<TM>()),
                       choose_param(get_param(aParams,internal_np::get_cost_policy), LindstromTurk_cost<TM>()),
                       choose_param(get_param(aParams,internal_np::get_placement_policy), LindstromTurk_placement<TM>()),
                       choose_param(get_param(aParams,vis), Dummy_visitor()));
}

template<class TM, class ShouldStop, class GT, class P, class T, class R>
int edge_collapse(TM& aSurface,
                  ShouldStop const& aShould_stop,
                  cgal_bgl_named_params<P,T,R> const& aParams)
{
  using boost::choose_param;
  using boost::choose_const_pmap;
  using boost::get_param;

  LindstromTurk_params lPolicyParams;
  internal_np::graph_visitor_t vis = internal_np::graph_visitor_t();

  return edge_collapse(aSurface, aShould_stop,
                       choose_const_pmap(get_param(aParams,internal_np::vertex_index),aSurface,boost::vertex_index),
                       choose_const_pmap(get_param(aParams,internal_np::vertex_point),aSurface,boost::vertex_point),
                       choose_const_pmap(get_param(aParams,internal_np::halfedge_index),aSurface,boost::halfedge_index),
                       choose_param(get_param(aParams,internal_np::edge_is_constrained),No_constrained_edge_map<TM>()),
                       choose_param(get_param(aParams,internal_np::get_cost_policy), LindstromTurk_cost<TM>()),
                       choose_param(get_param(aParams,internal_np::get_placement_policy), LindstromTurk_placement<TM>()),
                       choose_param(get_param(aParams,vis), Dummy_visitor()));
}

template<class TM, class ShouldStop>
int edge_collapse(TM& aSurface, const ShouldStop& aShould_stop)
{
  return edge_collapse(aSurface,aShould_stop, CGAL::parameters::halfedge_index_map(get(boost::halfedge_index,aSurface)));
}

template<class TM, class ShouldStop, class GT>
int edge_collapse(TM& aSurface, const ShouldStop& aShould_stop)
{
  return edge_collapse(aSurface,aShould_stop, CGAL::parameters::halfedge_index_map(get(boost::halfedge_index,aSurface)));
}

} // namespace Surface_mesh_simplification

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_H
