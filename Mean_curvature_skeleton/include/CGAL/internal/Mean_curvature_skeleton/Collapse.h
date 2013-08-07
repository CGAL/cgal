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

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <boost/graph/graph_traits.hpp>

namespace CGAL {
namespace internal {

/**
* Test if it is topologicaly feasible to collapse the given edge
* @pre the polyhedron is a triangular mesh
* @param polyhedron the mesh containing the given edge
* @param v0v1 the edge to be collapsed
*/
template<class Polyhedron>
bool is_collapse_ok(Polyhedron& polyhedron, 
                    typename boost::graph_traits<Polyhedron>::edge_descriptor v0v1)
{
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor	         vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor	           edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::in_edge_iterator           in_edge_iterator;

  edge_descriptor v1v0 = v0v1->opposite();
  vertex_descriptor v0 = boost::target(v1v0, polyhedron);
  vertex_descriptor v1 = boost::source(v1v0, polyhedron);

  vertex_descriptor vv, vl, vr;
  edge_descriptor  h1, h2;

  // the edges v1-vl and vl-v0 must not be both boundary edges
  if (!(v0v1->is_border()))
  {
    vl = boost::target(v0v1->next(), polyhedron);
    h1 = v0v1->next();
    h2 = h1->next();
    if (h1->opposite()->is_border() && h2->opposite()->is_border())
    {
      return false;
    }
  }

  // the edges v0-vr and vr-v1 must not be both boundary edges
  if (!(v1v0->is_border()))
  {
    vr = boost::target(v1v0->next(), polyhedron);
    h1 = v1v0->next();
    h2 = h1->next();
    if (h1->opposite()->is_border() && h2->opposite()->is_border())
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
  if (is_border(polyhedron, v0) 
   && is_border(polyhedron, v1)
   && !(v0v1->is_border()) 
   && !(v1v0->is_border()))
  {
    return false;
  }

  // test intersection of the one-rings of v0 and v1
  in_edge_iterator eb, ee;
  for (boost::tie(eb, ee) = boost::in_edges(v0, polyhedron); eb != ee; ++eb)
  {
    vv = boost::source(*eb, polyhedron);
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

template<class Polyhedron>
bool find_halfedge(Polyhedron& polyhedron, 
                   typename boost::graph_traits<Polyhedron>::vertex_descriptor vi,
                   typename boost::graph_traits<Polyhedron>::vertex_descriptor vj)
{
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor          vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor            edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::in_edge_iterator           in_edge_iterator;

  in_edge_iterator eb, ee;
  for (boost::tie(eb, ee) = boost::in_edges(vj, polyhedron); eb != ee; ++eb)
  {
    vertex_descriptor vv = boost::source(*eb, polyhedron);
    if (vv == vi)
    {
      return true;
    }
  }
  return false;
}

template<class Polyhedron>
bool is_border(Polyhedron& polyhedron,
               typename boost::graph_traits<Polyhedron>::vertex_descriptor aV)
{
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor          vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor            edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::in_edge_iterator           in_edge_iterator;

  bool rR = false;

  in_edge_iterator eb, ee;
  for (boost::tie(eb, ee) = boost::in_edges(aV, polyhedron); eb != ee; ++eb)
  {
    edge_descriptor lEdge = *eb;
    if (is_undirected_edge_a_border(lEdge))
    {
      rR = true;
      break;
    }
  }

  return rR;
}

template<class Edge_descriptor>
bool is_undirected_edge_a_border(Edge_descriptor aEdge)
{
  return aEdge->is_border() || aEdge->opposite()->is_border();
}

} //namespace internal
} //namespace CGAL

/// @endcond

#endif //CGAL_MCFSKEL_COLLAPSE_H
