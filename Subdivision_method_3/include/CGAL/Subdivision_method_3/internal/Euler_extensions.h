// ======================================================================
//
// Copyright (c) 2005-2017 GeometryFactory (France).  All Rights Reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s): Le-Jeng Shiue <Andy.Shiue@gmail.com>
//
// ======================================================================

#ifndef CGAL_SUBDIVISION_METHOD_3_INTERNAL_EULER_EXTENSIONS_H
#define CGAL_SUBDIVISION_METHOD_3_INTERNAL_EULER_EXTENSIONS_H

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/iterator.h>

namespace CGAL {

namespace Subdivision_method_3 {

namespace internal {

template<typename PolygonMesh>
typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
insert_vertex_return_edge(PolygonMesh& p,
                          typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h) {
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor hopp = opposite(h, p);
  halfedge_descriptor r = Euler::split_vertex(prev(hopp, p), h, p);
  if (! is_border(h, p))
    set_halfedge(face(h, p), r, p);
  if (! is_border(hopp, p))
    set_halfedge(face(hopp, p), hopp, p);
  return opposite(r, p);
}

/** Insert a new vertex into a helfedge h (a--b)

    Precondition:
         h
    a <-----> b
        -h
    h is the halfedge connecting vertex a to b
    -h is the opposite halfedge connecting b to a

    Postcondition:
         h         r
    a <-----> V <-----> b
        -h         -r
    V is the return vertex whose geometry is UNDEFINED.
    -r is the returned halfedge that is pointing to V
*/
template<typename PolygonMesh>
typename boost::graph_traits<PolygonMesh>::vertex_descriptor
insert_vertex(PolygonMesh& p,
              typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h) {
  return target(insert_vertex_return_edge(p, h), p);
}

/** Insert a new edge (two halfedges) between the two vertices

    Precondition:
    vertex a and b are in the SAME facet and do NOT connect to each other

    Postcondition:
           H
    a <----------> b
    H is the return halfedge connecting vertex a to b.
*/
template<typename PolygonMesh>
typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
insert_edge(PolygonMesh& p,
            typename boost::graph_traits<PolygonMesh>::halfedge_descriptor a,
            typename boost::graph_traits<PolygonMesh>::halfedge_descriptor b) {
  return Euler::split_face(a, b, p);
}

} // namespace internal

} // namespace Subdivision_method_3

} //namespace CGAL

#endif // CGAL_SUBDIVISION_METHOD_3_INTERNAL_EULER_EXTENSIONS_H
