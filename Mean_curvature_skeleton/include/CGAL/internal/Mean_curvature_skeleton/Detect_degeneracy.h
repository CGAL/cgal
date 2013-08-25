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

#ifndef CGAL_MCFSKEL_DETECT_DEGENERACY_H
#define CGAL_MCFSKEL_DETECT_DEGENERACY_H

/// @cond CGAL_DOCUMENT_INTERNAL

/** 
 * @file Detect_degeneracy.h
 * @brief This file contains functions to detect degeneracy at a given vertex.
 *
 * The approach is to count the Euler characteristics within a small geodesic 
 * distance at the given vertex. If it is not equal to one, which is the case
 * for disk topology, the vertex is considered to be degenerate.
 */

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <boost/graph/graph_traits.hpp>

#include <queue>

namespace CGAL {
namespace internal {

template<class HalfedgeGraph>
bool is_vertex_degenerate(HalfedgeGraph& polyhedron,
                          typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor root,
                          double edgelength_TH)
{
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor          vertex_descriptor;
  typedef typename boost::graph_traits<HalfedgeGraph>::edge_descriptor            edge_descriptor;
  typedef typename boost::graph_traits<HalfedgeGraph>::out_edge_iterator		      out_edge_iterator;
  typedef typename HalfedgeGraph::Face_handle                                     Face_handle;
  typedef typename HalfedgeGraph::Halfedge_around_facet_circulator                Halfedge_facet_circulator;

  std::set<vertex_descriptor> vertices_in_disk;
  std::set<edge_descriptor> edges_in_disk;
  std::set<Face_handle> faces_in_disk;

  vertices_in_disk.clear();
  search_vertices_in_disk(polyhedron, root, vertices_in_disk, edgelength_TH);

  typename std::set<vertex_descriptor>::iterator v_iter;
  for (v_iter = vertices_in_disk.begin(); v_iter != vertices_in_disk.end(); ++v_iter)
  {
    vertex_descriptor vd = *v_iter;
    out_edge_iterator e, e_end;
    for (boost::tie(e, e_end) = boost::out_edges(vd, polyhedron); e != e_end; ++e)
    {
      edge_descriptor ed = *e;
      edge_descriptor ed_op = ed->opposite();
      vertex_descriptor target = boost::target(ed, polyhedron);
      if (vertices_in_disk.find(target) != vertices_in_disk.end())
      {
        edges_in_disk.insert(ed);
        edges_in_disk.insert(ed_op);
      }
      Face_handle f = ed->face();
      Halfedge_facet_circulator j = f->facet_begin();
      bool in = true;
      do
      {
        vertex_descriptor v = j->vertex();
        if (vertices_in_disk.find(v) == vertices_in_disk.end())
        {
          in = false;
          break;
        }
      } while (++j != f->facet_begin());

      if (in)
      {
        faces_in_disk.insert(f);
      }
    }
  }

  int V = vertices_in_disk.size();
  int E = edges_in_disk.size() / 2;
  int F = faces_in_disk.size();
  int euler = V + F - E;
  if (euler != 1)
  {
    return true;
  }
  return false;
}

template<class HalfedgeGraph>
void search_vertices_in_disk(HalfedgeGraph& polyhedron,
                             typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor root,
                             std::set<typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor>& vertices_in_disk,
                             double edgelength_TH)
{
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor          vertex_descriptor;
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_iterator            vertex_iterator;
  typedef typename boost::graph_traits<HalfedgeGraph>::edge_descriptor            edge_descriptor;
  typedef typename boost::graph_traits<HalfedgeGraph>::out_edge_iterator		      out_edge_iterator;

  std::map<vertex_descriptor, bool> vertex_visited;

  std::queue<vertex_descriptor> Q;
  Q.push(root);
  vertices_in_disk.insert(root);
  vertex_visited[root] = true;

  double dist_TH = edgelength_TH;
  while (!Q.empty())
  {
    vertex_descriptor v = Q.front();
    Q.pop();

    out_edge_iterator e, e_end;
    for(boost::tie(e, e_end) = boost::out_edges(v, polyhedron); e != e_end; ++e)
    {
      edge_descriptor ed = *e;

      vertex_descriptor new_v = boost::target(ed, polyhedron);
      if (vertex_visited.find(new_v) == vertex_visited.end())
      {
        double distance = sqrtf(squared_distance(new_v->point(), root->point()));
        if (distance < dist_TH)
        {
          vertex_visited[new_v] = true;
          Q.push(new_v);
          vertices_in_disk.insert(new_v);
        }
      }
    }
  }
}

} //namespace internal
} //namespace CGAL

/// @endcond

#endif //CGAL_MCFSKEL_DETECT_DEGENERACY_H
