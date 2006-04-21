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
// $URL: $
// $Id: $
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
//   Candidate vertex-pairs are selected via a user-defined selection property map.
//   Selected vertex-pairs are sorted according to a user-defined cost property map.
//   Each vertex-pair removed is replaced by a new vertex whose location is given by a user-defined vertex placement function.
//   The simplification continues until a user-defined stoping conidition is verified.
// 
//   Template Parameters:
//
//   TSM : the triangulated surface mesh type (a Polyhedron_3 for example)
//   SelectionMap: the property map used to select candidate vertex-pairs.
//   CostMap: the property map providing the collapsing-cost of each selected vertex-pair.
//   VertexPlacement: the functor that constructs the 3D point of placement of the new vertex replacing each removed vertex-pair.
//   StopPred: the stopping condition predicate.
//
//   NOTE: The CostMap and the VertexPlacement functors return optional<> values. This allows any of those functions to reject a 
//         vertex-pair becasue it's cost is too high or uncomputable or the vertex cannot be placed in any way that satisfies the constriants
//         required. Consequently, an accepted vertex-pair (via the SelectMap) might not be removed in the end. 
//         Vertex-pairs are only removed if they are selected _and_ have a well defined collapsing cost and new vertex location.
//
//   Returns the number of vertex-pairs removed or -1 if there was an error (like the surface not being a valid triangulated surface mesh)
//       
template<class TSM,class SelectionMap,class CostMap,class VertexPlacement,class StopPred>
std::size_t vertex_pair_collapse ( TSM&                      aSurface
                                 , SelectionMap       const& aSelectMap
                                 , CostMap            const& aCostMap 
                                 , VertexPlacement    const& aVertexPlacement
                                 , StopPred           const& aStopPred
                                 , bool                      aIncludeNonEdgePairs = false
                                 ) 
{
  if ( is_valid_triangulated_surface_mesh(aSurface) )
  {
    typedef VertexPairCollapse<TSM,SelectionMap,CostMap,VertexPlacement,StopPred> Algorithm ;
    Algorithm algorithm(aSurface,aSelectMap,aCostMap,aVertexPlacement,aStopPred,aIncludeNonEdgePairs) ;
    return algorithm.run();
  }
  else return -1 ;
}                          


} } // namespace Triangulated_surface_mesh::Simplification


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_H //
// EOF //
 
