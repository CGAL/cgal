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

#include <CGAL/Surface_mesh_simplification/Detail/Edge_collapse.h>

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
//   via a property map.
//
//   For each non-fixed edge in the mesh, a "collapse data" record is constructed by calling the user-supplied function
//   "SetCollapseData", passing the edge and the specified parameters object (ParamsToSetCollapseData).
//   The edge is then associated with its collapse data via a property map.
//
//   The user-supplied function "GetCost" is called, for each non-fixed edge, passing its associated collapse data.
//   This function returns a value which defines the collapsing cost of the edge. Edges with a lower cost are collapsed first.
//     
//   When a non-fixed edge is selected for removal, a user-supplied function "GetNewVertexPoint" is called, 
//   passing its associated collapse data.
//   This function returns a Point_3 which defines the coordinates of the single vertex that "replaces" the collapsed edge.
// 
//   The simplification continues until there are no more non-fixed edges to collapse or the user-defined function "ShouldStop"
//   returns true.
//
//   NOTE: The functions GetCost and GetNewVertexPoint return 'optional values'. This is to allow the function to reject an
//         edge becasue it's cost is too high or uncomputable or the new vertex point cannot be placed in any way that
//         satisfies the constriants required by method used in these functions.
//         Consequently, not all non-fixed edges are neccessarily removed.
//
//   This global function returns the number of edges or -1 if there was an error 
//   (like the surface not being a valid triangulated surface mesh)
//       
template<class TSM,class Params,class SetCollapseData,class GetCost,class GetNewVertexPoint,class ShouldStop, class Visitor>
int edge_collapse ( TSM&                     aSurface
                  , Params const*            aParams // Can be NULL
                  , SetCollapseData const&   aSet_collapse_data
                  , GetCost         const&   aGet_cost 
                  , GetNewVertexPoint const& aGet_new_vertex_point
                  , ShouldStop      const&   aShould_stop
                  , Visitor*                 aVisitor = 0
                  ) 
{
  if ( is_valid_triangulated_surface_mesh(aSurface) )
  {
    typedef EdgeCollapse<TSM,Params,SetCollapseData,GetCost,GetNewVertexPoint,ShouldStop,Visitor> Algorithm ;
    Algorithm algorithm(aSurface
                       ,aParams
                       ,aSet_collapse_data
                       ,aGet_cost
                       ,aGet_new_vertex_point
                       ,aShould_stop
                       ,aVisitor
                       ) ;
    return algorithm.run();
  }
  else return -1 ;
}                          


} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_H //
// EOF //
 
