// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
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
//
// Author(s)     : Xiang Gao <gaox@ethz.ch>
//

#ifndef CGAL_MCFSKEL_COLLAPSE_H
#define CGAL_MCFSKEL_COLLAPSE_H

/// @cond CGAL_DOCUMENT_INTERNAL

/** 
 * @file Collapse.h
 * @brief This file contains the helper functions for collapsing short edges.
 * 
 * The approach is based on the functions from surface_mesh.
 */

#include <boost/graph/graph_traits.hpp>

namespace CGAL {
namespace internal {

/**
* Test if it is topologicaly feasible to collapse the given edge.
*
* @pre the polyhedron is a triangular mesh
* @param polyhedron the mesh containing the given edge
* @param v0v1 the edge to be collapsed
*/
template<class HalfedgeGraph>
bool is_collapse_ok(HalfedgeGraph& polyhedron, 
                    typename boost::graph_traits<HalfedgeGraph>::edge_descriptor v0v1)
{
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor	         vertex_descriptor;
  typedef typename boost::graph_traits<HalfedgeGraph>::edge_descriptor	           edge_descriptor;
  typedef typename boost::graph_traits<HalfedgeGraph>::in_edge_iterator            in_edge_iterator;

  edge_descriptor v1v0 = v0v1->opposite();
  vertex_descriptor v0 = target(v1v0, polyhedron);
  vertex_descriptor v1 = source(v1v0, polyhedron);

  vertex_descriptor vv, vl, vr;
  edge_descriptor  h1, h2;

  // the edges v1-vl and vl-v0 must not be both boundary edges
  if (!(is_border(v0v1, polyhedron)))
  {
    vl = target(v0v1->next(), polyhedron);
    h1 = v0v1->next();
    h2 = h1->next();
    if ( is_border(opposite(h1, polyhedron), polyhedron) &&
         is_border(opposite(h2, polyhedron), polyhedron) )
    {
      return false;
    }
  }

  // the edges v0-vr and vr-v1 must not be both boundary edges
  if (!is_border(v1v0, polyhedron))
  {
    vr = target(v1v0->next(), polyhedron);
    h1 = v1v0->next();
    h2 = h1->next();
    if ( is_border(opposite(h1, polyhedron), polyhedron) &&
         is_border(opposite(h2, polyhedron), polyhedron) )
    {
      return false;
    }
  }

  // if vl and vr are equal or both invalid -> fail
  if (vl == vr)
  {
    return false;
  }

  // edge between two boundary vertices should be a boundary edge
  if (is_border(v0, polyhedron) 
   && is_border(v1, polyhedron)
   && !(is_border(v0v1, polyhedron))
   && !(is_border(v1v0,polyhedron)))
  {
    return false;
  }

  // test intersection of the one-rings of v0 and v1
  in_edge_iterator eb, ee;
  for (boost::tie(eb, ee) = in_edges(v0, polyhedron); eb != ee; ++eb)
  {
    vv = source(*eb, polyhedron);
    if (vv != v1 && vv != vl && vv != vr)
    {
      if (find_halfedge(polyhedron, vv, v1))
      {
        return false;
      }
    }
  }

  // passed all tests
  return true;
}

/**
* Find the edge which containg given vertices.
*
* @param polyhedron the mesh containing the edge.
* @param vi one vertex incident to the edge.
* @param vj the other vertex incident to the edge.
*/
template<class HalfedgeGraph>
bool find_halfedge(HalfedgeGraph& polyhedron, 
                   typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor vi,
                   typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor vj)
{
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor          vertex_descriptor;
  typedef typename boost::graph_traits<HalfedgeGraph>::edge_descriptor            edge_descriptor;
  typedef typename boost::graph_traits<HalfedgeGraph>::in_edge_iterator           in_edge_iterator;

  in_edge_iterator eb, ee;
  for (boost::tie(eb, ee) = in_edges(vj, polyhedron); eb != ee; ++eb)
  {
    vertex_descriptor vv = source(*eb, polyhedron);
    if (vv == vi)
    {
      return true;
    }
  }
  return false;
}

} //namespace internal
} //namespace CGAL

/// @endcond

#endif //CGAL_MCFSKEL_COLLAPSE_H
