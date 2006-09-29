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

namespace Surface_mesh_simplification
{


//
// collapse_triangulation_edge euler operator (topological operation only).
//
// This operator collapses the edge p-q replacing it with one single vertex connected to the link of p-q.
//
// The net effect of the operator is equivalent to removing one of the vertices and re-triangulating the resulting hole.
//
// The actual collapse operation only removes up to 2 faces, 3 edges and 1 vertex, as follows.
//
// The top face, p-q-t, to the left of p->q, if any, is removed.
// The bottom face, q-p-b, to the left of q->p (right of p->q), if any, is removed.
//
// If there is a top-left face, that is, if the edge p->t is NOT a border edge, the top face is removed by joining it
// with the top-left face (this creates a 4 sided face: p-q-t-l).
// If p->t IS a border edge then the top face is simply erased.
//
// If there is a bottom-right face, that is, if the edge q->b is NOT a border edge, the bottom face is removed by joining it
// with the bottom-right face (this creates a 4 sided face: q-p-b-r)
// If q->b IS a border edge then the bottom face is simply erased.
//
// If there is a top face edge p->t is removed. 
// If there is a bottom face edge q->b is removed. 
// 
// One of the vertices (p or q) is removed. 
// If there is no top-left face so the top face is directly erased AND there is no bottom face,
// that face erasure automatically removes vertex p. Likewise, if there is no bottom-right face 
// so the bottom face is directly erased AND there is no top face, that erasure automatically 
// removes vertex q.
// Directly erasing the top/bottom faces when the opposite face exists does not removes any vertex
// automatically, in which case vertex p is removed by joining p->q
//
// NOTES:
//
// This operator can only be called by collapsable edges (which satisfy the link condition), hence,
// there must exist at least the top face or a bottom face, and, if there is no top-left face there is a top-right face
// and likewise, if there is no bottom-right face there is a bottom-left face. (IOW vertices t/b must have degree >=3
// unless the top/bottom faces do not exists)
//
// The operator doesn't join the top face with the top-left face and the bottom-face with the bottom-left face
// (or both to the right) because if the top-left and bottom-left (or both right) faces are themselve adjacent,
// the first joint would introduce a degree 2 vertex.
// That is why the operator alternates left and right to join the top and bottom faces.
//
// PARAMETERS:
//  pq : the edge to collapse
//  pt : the edge shared between the top face and the top-left face (if any). If there is no top face this parameter is a null handle.
//  qb : the edge shared between the bottom face and the bottom-right face. If there is no bottom face this parameter is a null handle.
//
// RETURN VALUE: A handle to the vertex that IS NOT removed.
//
template<class ECM> struct Collapse_triangulation_edge ;

} // namespace Surface_mesh_simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_OPERATOR_H //
// EOF //
 
