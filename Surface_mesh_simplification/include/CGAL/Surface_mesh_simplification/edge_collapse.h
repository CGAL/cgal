// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_H 1

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/named_function_params.h>

#include <CGAL/Surface_mesh_simplification/Detail/Edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Vertex_is_fixed_property_map_always_false.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_params.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_set_cost_cache.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Cached_cost.h>

CGAL_BEGIN_NAMESPACE

namespace Surface_mesh_simplification
{


//
// Edge-collapse method:
//
//   Simplifies a triangulated surface mesh by iteratively selecting and collapsing edges.
//
//   Not all edges are selected for removal; those which cannot be removed are labled "fixed" 
//   and there are two kinds of fixed pairs: Intrisically fixed and explicitely fixed.
//   Intrinsically fixed edges are those which, if removed, would result in a topologically incosistent mesh.
//   (the algorithm automatically detects intrinsically fixed edges).
//   Explicitely fixed edges appear when the two vertices incident on the edge have been marked by the user as fixed,
//   via the "VertexIsFixedMap" property map.
//
//   For each non-fixed edge in the mesh, a "collapse data" record is constructed by calling the user-supplied function
//   "SetCollapseData".
//   The edge is then associated with its collapse data via the "EdgeExtraPtrMap" property map.
//
//   The user-supplied function "GetCost" is called, for each non-fixed edge.
//   This function returns a value which defines the collapsing cost of the edge. Edges with a lower cost are collapsed first.
//     
//   When a non-fixed edge is selected for removal, a user-supplied function "GetNewVertexPoint" is called.
//   This function returns a Point_3 which defines the coordinates of the single vertex that replaces the collapsed edge.
// 
//   The simplification continues until there are no more non-fixed edges to collapse or the user-defined function "ShouldStop"
//   returns true.
//
//   NOTE: The functions GetCost and GetNewVertexPoint return 'optional values'. This is to allow the function to reject an
//         edge becasue it's cost is too high or uncomputable or the new vertex point cannot be placed in any way that
//         satisfies the constriants required by method used in these functions.
//         Consequently, not all non-fixed edges are neccessarily removed.
//
//   This global function returns the number of edges collapsed.
//
template<class ECM
        ,class ShouldStop
        ,class VertexPointMap
        ,class VertexIsFixedMap
        ,class EdgeIndexMap
        ,class EdgeIsBorderMap
        ,class SetCache
        ,class GetCost
        ,class GetPlacement
        ,class CostParams
        ,class PlacementParams
        ,class Visitor
        >
int edge_collapse ( ECM&                    aSurface
                  , ShouldStop       const& aShould_stop
                  
                  // optional mesh information policies 
                  , VertexPointMap   const& aVertex_point_map     // defaults to get(Vertex_point,aSurface)
                  , VertexIsFixedMap const& aVertex_is_fixed_map  // defaults to Vertex_is_fixed_map_always_false<ECM>()
                  , EdgeIndexMap     const& aEdge_index_map       // defaults to get(Edge_index,aSurface) 
                  , EdgeIsBorderMap  const& aEdge_is_border_map   //  defaults to get(Edge_is_border,aSurface) 
                  
                  // optional strategy policies - defaults to LindstomTurk
                  , SetCache         const& aSet_cache
                  , GetCost          const& aGet_cost 
                  , GetPlacement     const& aGet_placement
                  , CostParams       const* aCostParams       // Can be NULL 
                  , PlacementParams  const* aPlacementParams  // Can be NULL 
                  
                  , Visitor*                aVisitor // Can be NULL
                  ) 
{
  typedef EdgeCollapse< ECM
                      , ShouldStop
                      , VertexPointMap
                      , VertexIsFixedMap
                      , EdgeIndexMap
                      , EdgeIsBorderMap
                      , SetCache
                      , GetCost
                      , GetPlacement
                      , CostParams
                      , PlacementParams
                      , Visitor
                      >
                      Algorithm 
                      ;
                      
  Algorithm algorithm( aSurface
                     , aShould_stop
                     , aVertex_point_map
                     , aVertex_is_fixed_map
                     , aEdge_index_map
                     , aEdge_is_border_map
                     , aSet_cache
                     , aGet_cost
                     , aGet_placement
                     , aCostParams
                     , aPlacementParams
                     , aVisitor
                     ) ;
                     
  return algorithm.run();
}                          


struct Dummy_visitor
{
  template<class ECM>
  void OnStarted( ECM& ) {} 
  
  template<class ECM>
  void OnFinished ( ECM& ) {} 
  
  template<class ECM>
  void OnStopConditionReached( ECM& ) {} 
  
  template<class Edge, class ECM>
  void OnCollected( Edge const&, bool, ECM& ) {}                
  
  template<class Edge, class ECM, class OFT, class Size_type>
  void OnSelected( Edge const&, ECM&, OFT const&, Size_type, Size_type ) {}                
  
  template<class Edge, class ECM, class OPoint>
  void OnCollapsing(Edge const&, ECM&, OPoint const& ) {}                
  
  template<class Edge, class ECM>
  void OnNonCollapsable(Edge const&, ECM& ) {}                
} ;

template<class ECM, class ShouldStop, class P, class T, class R>
int edge_collapse ( ECM& aSurface
                  , ShouldStop const& aShould_stop
                  , cgal_bgl_named_params<P,T,R> const& aParams 
                  ) 
{
  using boost::choose_param ;
  using boost::choose_const_pmap ;
  using boost::get_param ;
  
  LindstromTurk_params lPolicyParams ;
  
  boost::graph_visitor_t vis = boost::graph_visitor_t() ;
    
  return edge_collapse(aSurface
                      ,aShould_stop
                      ,choose_const_pmap     (get_param(aParams,vertex_point),aSurface,vertex_point)
                      ,choose_param          (get_param(aParams,vertex_is_fixed),Vertex_is_fixed_property_map_always_false<ECM>())
                      ,choose_const_pmap     (get_param(aParams,boost::edge_index),aSurface,boost::edge_index)
                      ,choose_const_pmap     (get_param(aParams,edge_is_border),aSurface,edge_is_border)
                      ,choose_param          (get_param(aParams,set_cache_policy), LindstromTurk_set_cost_cache<ECM>())
                      ,choose_param          (get_param(aParams,get_cost_policy), Cached_cost<ECM>())
                      ,choose_param          (get_param(aParams,get_placement_policy), LindstromTurk_placement<ECM>())
                      ,choose_param          (get_param(aParams,get_cost_policy_params), &lPolicyParams)
                      ,choose_param          (get_param(aParams,get_placement_policy_params), &lPolicyParams)
                      ,choose_param          (get_param(aParams,vis), ((Dummy_visitor*)0))
                      ) ;

}

template<class ECM, class ShouldStop>
int edge_collapse ( ECM& aSurface, ShouldStop const& aShould_stop ) 
{
  return edge_collapse(aSurface,aShould_stop, edge_index_map(get(boost::edge_index,aSurface)));
}

} // namespace Surface_mesh_simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_H //
// EOF //
 
