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
 */

#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <CGAL/boost/graph/iterator.h>
#include <cmath>
#include <queue>

namespace CGAL {
namespace internal {

/**
* Test if a given vertex is degenerate.
*
* The approach is to count the Euler characteristics within a small geodesic
* distance at the given vertex. If it is not equal to one, which is the case
* for disk topology, the vertex is considered to be degenerate.

* @param hg the mesh containing the given vertex
* @param root the given vertex
* @param min_edge_length the diameter of the geodesic disk
*/
template<class TriangleMesh, class TriangleMeshPointPMap, class Traits>
bool is_vertex_degenerate(TriangleMesh& hg,
                          TriangleMeshPointPMap& hg_point_pmap,
                          typename boost::graph_traits<TriangleMesh>::vertex_descriptor root,
                          double min_edge_length,
                          const Traits& traits)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor          vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor        halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor            edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor            face_descriptor;

  std::set<vertex_descriptor> vertices_in_disk;
  std::set<halfedge_descriptor> edges_in_disk;
  std::set<face_descriptor> faces_in_disk;

  vertices_in_disk.clear();
  search_vertices_in_disk(hg, hg_point_pmap, root, vertices_in_disk, min_edge_length, traits);

  BOOST_FOREACH(vertex_descriptor vd, vertices_in_disk)
  {
    BOOST_FOREACH(edge_descriptor ed, out_edges(vd, hg))
    {
      halfedge_descriptor hd = halfedge(ed, hg);
      halfedge_descriptor hd_op = opposite(hd, hg);
      vertex_descriptor tgt = target(hd, hg);
      if (vertices_in_disk.find(tgt) != vertices_in_disk.end())
      {
        edges_in_disk.insert(hd);
        edges_in_disk.insert(hd_op);
      }

      bool in = true;
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(hd, hg))
      {
        vertex_descriptor v = target(hd,hg);
        if (vertices_in_disk.find(v) == vertices_in_disk.end())
        {
          in = false;
          break;
        }
      }

      if (in)
      {
        faces_in_disk.insert(face(hd,hg));
      }
    }
  }

  std::size_t V = vertices_in_disk.size();
  std::size_t E = edges_in_disk.size() / 2;
  std::size_t F = faces_in_disk.size();
  std::size_t euler = V + F - E;
  if (euler != 1)
  {
    return true;
  }
  return false;
}

/**
* Find all the vertices within a geodesic disk.
*
* @param hg the mesh containing the vertices
* @param root the center of the geodesic disk
* @param vertices_in_disk containing the found vertices within the disk
* @param min_edge_length the diameter of the geodesic disk
*/
template<class TriangleMesh, class TriangleMeshPointPMap, class Traits>
void search_vertices_in_disk(TriangleMesh& hg,
                             TriangleMeshPointPMap& hg_point_pmap,
                             typename boost::graph_traits<TriangleMesh>::vertex_descriptor root,
                             std::set<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>& vertices_in_disk,
                             double min_edge_length,
                             const Traits& traits)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor          vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor            edge_descriptor;

  std::queue<vertex_descriptor> Q;
  Q.push(root);
  vertices_in_disk.insert(root);

  double dist_TH = min_edge_length;
  while (!Q.empty())
  {
    vertex_descriptor v = Q.front();
    Q.pop();

    BOOST_FOREACH(edge_descriptor ed, out_edges(v, hg))
    {
      vertex_descriptor new_v = target(ed, hg);
      if (!vertices_in_disk.count(new_v))
      {
        double distance = std::sqrt(traits.compute_squared_distance_3_object()(
                                      get(hg_point_pmap, new_v),
                                      get(hg_point_pmap, root)) );
        if (distance < dist_TH)
        {
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
