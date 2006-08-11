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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_OPERATOR_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_OPERATOR_H 1

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification 
{


//
// collapse_triangulation_edge euler operator (topological operation only).
//
// This operator removes from a triangulation 'aDS' 2 facets, 1 vertex and 3 edges, as follows:
//
// Consider 2 vertices 'p', 'q', 't' and 'b', such that the edge 'p-q' is shared by two facets
// 'p-q-t' and 'q-p-b', called top and bottom facets resp. (their orientation is unimportant)
//
// Consider additional vertices 'l' and 'r' such that 
// 'p-t-l' and 'q-r-t' are the other 2 facets adjacent to the top facet, called the top-left and top-right facets, resp.,
// and 
// 'l-b-p' and 'b-r-q' are the other 2 facets adjacent to the bottom facet, called the bottom-left and bottom-right facets, resp.
//
// In such a triangulated path there is a cycle of vertices around 'p' and 'q' called the "link" of the edge 'p-q'.
// In the minimal case, the link is the unordered cycle of vertices: l-b-r-t. It could also be that 'l' and 'r' are really a
// sequence of vertices instead of just 1 vertex.
// As the link has at least 4 vertices there are at least 6 facets inside the link: 
//   bottom-left  (l-b-p)
//   bottom       (q-p-b)
//   bottom-right (q-b-r)
//   top-right    (q-r-t) 
//   top          (p-q-t)
//   top-left     (l-p-t)
// (the vertices are enumerated in an arbiraty order as that is unimportant here)
// 
// The operator is passed 3 DIRECTED edges: 'p-q', 'p-t' and 'q-b' as proceeds with the following 3 steps:
//
// (1) Merges the top facet with the top-left facet keeping the top-left facet; that is,
//     removes the top-facet 'q-p-t', the edge 'p-t' and 
//     redefines the top-left facet to be 'l-p-q-t' instead of 'l-p-t'.
//
// (2) Merges the bottom facet with the bottom-right facet keeping the bottom-right facet; that is,
//     removes the botton-facet 'q-p-t', the edge 'q-b' and 
//     redefines the bottom-right facet to be 'r-q-p-b' instead of 'r-q-b'.
//
// (3) Joins the vertices 'p' and 'q' keeping the vertex 'q'; that is
//     removes the edge 'p-q' and the vertex 'p' redefining the top-left facet to be 'l-q-t' instead of 'l-p-q-t'
//     and the bottom-right facet to be 'r-q-b' instead of 'r-q-p-b'
//
// The net result is a valid triangulation but the intermediate results from steps 1 and 2 are not as the top-left
// and bottom-right facets are temporarily 4-sided until 'p' and 'q' is joint.
//
// The operator merges the bottom facets with the bottom-right facet instead of the bottom-left facet becasue in the
// minimal link case (that is, just 'l-b-r-t'), the top-left and bottom-left facets are adjacent, so merging one of 
// them prevents the other to be merged becuase doing so would introduce degree-2 vertex (which is illegal in most DS)
//
// The operator is required to erase from the aDS the following: vertex 'p', edges 'pq', 'pt' and 'qb', and
// facets 'pqt' and 'qpb'. It is also required NOT to erase anything else.
// 
// The code in this primary template simply forwards the operations to the DS. 
// Thus, the DS must support the join_facet(edge) and join_vertex(edge) operations. If not, an specialization is required.
//
template<class DS_>
struct Collapse_triangulation_edge
{
  typedef DS_ DS ;
  
  typedef typename boost::graph_traits<DS>::edge_descriptor edge_descriptor ;
  
  void operator() ( edge_descriptor const& pq
                  , edge_descriptor const& pt
                  , edge_descriptor const& qb
                  , DS&                    aDS 
                  ) const
  {
    CGAL_precondition( pt->vertex()->vertex_degree() >= 3 || qb->vertex()->vertex_degree() >= 3 ) ;
  
    if ( pt->vertex()->vertex_degree() >= 3 )
      aDS.join_facet (pt);
    
    if ( qb->vertex()->vertex_degree() >= 3 )
      aDS.join_facet (qb);
    
    aDS.join_vertex(pq);
    
    CGAL_expensive_postcondition(aDS.is_valid());
  }
  
} ;

} } // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_OPERATOR_H //
// EOF //
 
