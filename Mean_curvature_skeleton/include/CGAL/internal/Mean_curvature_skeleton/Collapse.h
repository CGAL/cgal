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
#include <boost/foreach.hpp>

namespace CGAL {
namespace internal {

/**
* Test if it is topologicaly feasible to collapse the given edge.
*
* @pre the hg is a triangular mesh
* @param hg the mesh containing the given edge
* @param v0v1 the edge to be collapsed
*/
/// \todo see the difference with is_collapse_topologically_valid
template<class TriangleMesh>
bool is_collapse_ok(TriangleMesh& hg,
                    typename boost::graph_traits<TriangleMesh>::halfedge_descriptor v0v1)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor          vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor        halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor            edge_descriptor;

  halfedge_descriptor v1v0 = opposite(v0v1, hg);
  vertex_descriptor v0 = target(v1v0, hg);
  vertex_descriptor v1 = source(v1v0, hg);

  vertex_descriptor vv, vl, vr;
  halfedge_descriptor  h1, h2;

  // the edges v1-vl and vl-v0 must not be both boundary edges
  if (!(is_border(v0v1, hg)))
  {
    vl = target(v0v1->next(), hg);
    h1 = next(v0v1, hg);
    h2 = next(h1, hg);
    if ( is_border(opposite(h1, hg), hg) &&
         is_border(opposite(h2, hg), hg) )
    {
      return false;
    }
  }

  // the edges v0-vr and vr-v1 must not be both boundary edges
  if (!is_border(v1v0, hg))
  {
    vr = target( next(v1v0, hg), hg);
    h1 = next(v1v0, hg);
    h2 = next(h1, hg);
    if ( is_border(opposite(h1, hg), hg) &&
         is_border(opposite(h2, hg), hg) )
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
  if (is_border(v0, hg)
   && is_border(v1, hg)
   && !(is_border(v0v1, hg))
   && !(is_border(v1v0,hg)))
  {
    return false;
  }

  // test intersection of the one-rings of v0 and v1
  BOOST_FOREACH(edge_descriptor ed, in_edges(v0, hg))
  {
    vv = source(ed, hg);
    if (vv != v1 && vv != vl && vv != vr)
    {
      if (find_halfedge(hg, vv, v1))
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
* @param hg the mesh containing the edge.
* @param vi one vertex incident to the edge.
* @param vj the other vertex incident to the edge.
*/
template<class TriangleMesh>
bool find_halfedge(TriangleMesh& hg,
                   typename boost::graph_traits<TriangleMesh>::vertex_descriptor vi,
                   typename boost::graph_traits<TriangleMesh>::vertex_descriptor vj)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor          vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor            edge_descriptor;

  BOOST_FOREACH(edge_descriptor ed, in_edges(vj, hg))
  {
    vertex_descriptor vv = source(ed, hg);
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
