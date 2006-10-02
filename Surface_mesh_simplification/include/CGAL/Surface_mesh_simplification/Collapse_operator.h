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
// Returns a handle to the vertex that IS NOT removed.
//
template<class ECM> struct Collapse_triangulation_edge ;

} // namespace Surface_mesh_simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_COLLAPSE_OPERATOR_H //
// EOF //
 
