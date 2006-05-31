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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_H 1

#include <CGAL/Surface_mesh_simplification/Vertex_pair_collapse.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification 
{

//
// Vertex-pair-collapse method:
//
//   Simplifies a triangulated surface mesh by iteratively selecting and removing vertex-pairs.
//
//   Not all vertex pairs are selected for removal; those which cannot be removed are labled "fixed" 
//   and there are two kinds of fixed pairs: Intrisically fixed and explicitely fixed.
//   Intrinsically fixed pairs are edges that, if removed, would result in a topologically incosistent mesh.
//   (the algorithm automatically detects intrinsically fixed edges).
//   Explicitely fixed pairs appear when both vertices in the pair have been marked by the user as fixed, via an external property map.
//
//   For each non-fixed vertex pair in the mesh, a "collapse data" record is constructed by calling the user-supplied function
//   "GetCollapseData", passing the pair and the specified parameters object (ParamsToGetCollapseData).
//   The pair is then associated with its collapse data.
//
//   The user-supplied function "GetCost" is called, for each non-fixed pair, passing its associated collapse data.
//   This function returns a value which defines the priority of the pairs. Pairs with a lower cost are removed first.
//     
//   When a non-fixed pair is selected for removal, a user-supplied function "GetNewVertexPoint" is called, 
//   passing its associated collapse data.
//   This function returns a Point_3 which is defines the coordinates of the single vertex that "replaces" the collapsed pair.
// 
//   The simplification continues until there no more non-fixed pairs to collapse or the user-defined function "ShouldStop" returns true.
//
//   NOTE: The functions GetCost and GetNewVertexPoint return 'optional values'. This is to allow the function to reject a 
//         vertex-pair becasue it's cost is too high or uncomputable or the new vertex point cannot be placed in any way that
//         satisfies the constriants required by method used in these functions.
//         Consequently, not all non-fixed vertex-pairs are neccessarily removed.
//
//   This global function returns the number of vertex-pairs removed or -1 if there was an error 
//   (like the surface not being a valid triangulated surface mesh)
//       
template<class TSM,class GetCollapseData,class ParamsToGetCollapseData,class GetCost,class GetNewVertexPoint,class ShouldStop>
int vertex_pair_collapse ( TSM&                           aSurface
                         , GetCollapseData const&         aGet_collapse_data
                         , ParamsToGetCollapseData const* aParamsToGetCollapseData // Can be NULL
                         , GetCost         const&         aGet_cost 
                         , GetNewVertexPoint const&       aGet_new_vertex_point
                         , ShouldStop      const&         aShould_stop
                         , bool                           aIncludeNonEdgePairs = false
                         ) 
{
  if ( is_valid_triangulated_surface_mesh(aSurface) )
  {
    typedef VertexPairCollapse<TSM,GetCollapseData,GetCost,GetNewVertexPoint,ShouldStop> Algorithm ;
    Algorithm algorithm(aSurface
                       ,aGet_collapse_data
                       ,aParamsToGetCollapseData
                       ,aGet_cost
                       ,aGet_new_vertex_point
                       ,aShould_stop
                       ,aIncludeNonEdgePairs
                       ) ;
    return algorithm.run();
  }
  else return -1 ;
}                          


} } // namespace Triangulated_surface_mesh::Simplification


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_H //
// EOF //
 
