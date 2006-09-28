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

#include <CGAL/boost/graph/BGL_properties.h>

#include <CGAL/Surface_mesh_simplification/Detail/Edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Edge_extra_pointer_map_stored.h>
#include <CGAL/Surface_mesh_simplification/Vertex_is_fixed_map_always_false.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_params.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_set_partial_collapse_data.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Cached_cost.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification { namespace Edge_collapse
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
template<class TSM
        ,class ShouldStop
        ,class EdgeExtraPtrMap
        ,class VertexIsFixedMap
        ,class SetCollapseData
        ,class GetCost
        ,class GetPlacement
        ,class CostParams
        ,class PlacementParams
        ,class Visitor
        >
int edge_collapse ( TSM&                    aSurface
                  , ShouldStop       const& aShould_stop
                  
                  // optional mesh information policies 
                  , EdgeExtraPtrMap  const& aEdge_extra_ptr_map   // defaults to Edge_extra_pointer_map_stored<TSM>
                  , VertexIsFixedMap const& aVertex_is_fixed_map  // defaults to Vertex_is_fixed_map_always_false<TSM>
                  
                  // optional strategy policies - defaults to LindstomTurk
                  , SetCollapseData  const& aSet_collapse_data
                  , GetCost          const& aGet_cost 
                  , GetPlacement     const& aGet_placement
                  , CostParams       const* aCostParams       // Can be NULL 
                  , PlacementParams  const* aPlacementParams  // Can be NULL 
                  
                  , Visitor*                aVisitor // Can be NULL
                  ) 
{
  typedef EdgeCollapse<TSM
                      ,ShouldStop
                      ,EdgeExtraPtrMap
                      ,VertexIsFixedMap
                      ,SetCollapseData
                      ,GetCost
                      ,GetPlacement
                      ,CostParams
                      ,PlacementParams
                      ,Visitor
                      >
                      Algorithm 
                      ;
                      
  Algorithm algorithm(aSurface
                     ,aShould_stop
                     ,aEdge_extra_ptr_map
                     ,aVertex_is_fixed_map
                     ,aSet_collapse_data
                     ,aGet_cost
                     ,aGet_placement
                     ,aCostParams
                     ,aPlacementParams
                     ,aVisitor
                     ) ;
                     
  return algorithm.run();
}                          


struct Dummy_visitor
{
  template<class TSM>
  void OnStarted( TSM& ) {} 
  
  template<class TSM>
  void OnFinished ( TSM& ) {} 
  
  template<class TSM>
  void OnStopConditionReached( TSM& ) {} 
  
  template<class Edge, class TSM>
  void OnCollected( Edge const&, bool, TSM& ) {}                
  
  template<class Edge, class TSM, class OFT, class Size_type>
  void OnSelected( Edge const&, TSM&, OFT const&, Size_type, Size_type ) {}                
  
  template<class Edge, class TSM, class OPoint>
  void OnCollapsing(Edge const&, TSM&, OPoint const& ) {}                
  
  template<class Edge, class TSM>
  void OnNonCollapsable(Edge const&, TSM& ) {}                
} ;

template<class TSM
        ,class ShouldStop
        ,class EdgeExtraPtrMap
        ,class VertexIsFixedMap
        ,class SetCollapseData
        ,class GetCost
        ,class GetPlacement
        ,class CostParams
        ,class PlacementParams
        >
int edge_collapse ( TSM&                    aSurface
                  , ShouldStop       const& aShould_stop
                  , EdgeExtraPtrMap  const& aEdge_extra_ptr_map
                  , VertexIsFixedMap const& aVertex_is_fixed_map
                  , SetCollapseData  const& aSet_collapse_data
                  , GetCost          const& aGet_cost 
                  , GetPlacement     const& aGet_placement
                  , CostParams       const* aCostParams
                  , PlacementParams  const* aPlacementParams
                  ) 
{
  return edge_collapse(aSurface 
                      ,aShould_stop
                      ,aEdge_extra_ptr_map
                      ,aVertex_is_fixed_map
                      ,aSet_collapse_data
                      ,aGet_cost
                      ,aGet_placement
                      ,aCostParams
                      ,aPlacementParams
                      ,((Dummy_visitor*)0)
                      );
}                          

template<class TSM
        ,class ShouldStop
        ,class EdgeExtraPtrMap
        ,class VertexIsFixedMap
        ,class SetCollapseData
        ,class GetCost
        ,class GetPlacement
        >
int edge_collapse ( TSM&                    aSurface
                  , ShouldStop       const& aShould_stop
                  , EdgeExtraPtrMap  const& aEdge_extra_ptr_map
                  , VertexIsFixedMap const& aVertex_is_fixed_map
                  , SetCollapseData  const& aSet_collapse_data
                  , GetCost          const& aGet_cost 
                  , GetPlacement     const& aGet_placement
                  ) 
{
  typename ExtractCostParamsType     <GetCost     ,SetCollapseData>::type      cost_params ;
  typename ExtractPlacementParamsType<GetPlacement,SetCollapseData>::type placement_params ;
  
  return edge_collapse(aSurface 
                      ,aShould_stop
                      ,aEdge_extra_ptr_map
                      ,aVertex_is_fixed_map
                      ,aSet_collapse_data
                      ,aGet_cost
                      ,aGet_placement
                      ,&cost_params
                      ,&placement_params
                      ,((Dummy_visitor*)0)
                      );
}                          

template<class TSM, class ShouldStop, class EdgeExtraPtrMap, class VertexIsFixedMap>
int edge_collapse ( TSM&                    aSurface
                  , ShouldStop       const& aShould_stop
                  , EdgeExtraPtrMap  const& aEdge_extra_ptr_map
                  , VertexIsFixedMap const& aVertex_is_fixed_map
                  ) 
{
  LindstromTurk_params params ;
  
  return edge_collapse(aSurface 
                      ,aShould_stop
                      ,aEdge_extra_ptr_map
                      ,aVertex_is_fixed_map
                      ,Set_partial_collapse_data_LindstromTurk<TSM>()
                      ,Cached_cost<TSM>()
                      ,LindstromTurk_placement<TSM>()
                      ,&params
                      ,&params
                      ,((Dummy_visitor*)0)
                      );
}                          

template<class TSM, class ShouldStop, class EdgeExtraPtrMap>
int edge_collapse ( TSM& aSurface, ShouldStop const& aShould_stop, EdgeExtraPtrMap const& aEdge_extra_ptr_map ) 
{
  LindstromTurk_params params ;
  
  return edge_collapse(aSurface 
                      ,aShould_stop
                      ,aEdge_extra_ptr_map
                      ,Vertex_is_fixed_map_always_false<TSM>()
                      ,Set_partial_collapse_data_LindstromTurk<TSM>()
                      ,Cached_cost<TSM>()
                      ,LindstromTurk_placement<TSM>()
                      ,&params
                      ,&params
                      ,((Dummy_visitor*)0)
                      );
}                          

template<class TSM, class ShouldStop>
int edge_collapse ( TSM& aSurface, ShouldStop const& aShould_stop ) 
{
  LindstromTurk_params params ;
  
  return edge_collapse(aSurface 
                      ,aShould_stop
                      ,Edge_extra_pointer_map_stored<TSM>()
                      ,Vertex_is_fixed_map_always_false<TSM>()
                      ,Set_partial_collapse_data_LindstromTurk<TSM>()
                      ,Cached_cost<TSM>()
                      ,LindstromTurk_placement<TSM>()
                      ,&params
                      ,&params
                      ,((Dummy_visitor*)0)
                      );
}                          

} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_H //
// EOF //
 
