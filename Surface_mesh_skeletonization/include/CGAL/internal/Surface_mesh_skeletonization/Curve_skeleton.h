// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Xiang Gao <gaox@ethz.ch>
//

#ifndef CGAL_MCFSKEL_CURVE_SKELETON_H
#define CGAL_MCFSKEL_CURVE_SKELETON_H

#include <CGAL/license/Surface_mesh_skeletonization.h>


/// @cond CGAL_DOCUMENT_INTERNAL

/**
 * @file Curve_skeleton.h
 * @brief This file contains the class used to turn a contracted mesh to a
 * curve skeleton.
 *
 */
#include <cmath>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/graph/copy.hpp>

// For debugging macro
#include <CGAL/internal/Surface_mesh_skeletonization/Debug.h>

namespace CGAL {
namespace internal {

template <class TriangleMesh, class VertexIndexMap,
          class HalfedgeIndexMap, class TriangleMeshPointPMap>
class Curve_skeleton
{
// Public types
public:

  // Geometric types
  typedef typename TriangleMesh::Traits      Kernel;
  typedef typename Kernel::Vector_3           Vector;
  typedef typename Kernel::Point_3            Point;

  // Repeat TriangleMesh types
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor           face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor           edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor       halfedge_descriptor;

// Data members
private:
  std::vector<std::vector<int> > edge_to_face;
  std::vector<std::vector<int> > edge_to_vertex;
  std::vector<std::vector<int> > vertex_to_edge;
  std::vector<std::vector<int> > face_to_edge;
  std::vector<bool> is_vertex_deleted;
  std::vector<bool> is_edge_deleted;
  std::vector<bool> is_face_deleted;

  std::vector<int> surface_vertex_id;
  // records the vertices collapsed to a given vertex
  std::vector<std::vector<int> > record;
  // vertex id mapped to vertex descriptor
  std::vector<vertex_descriptor> id_to_descriptor;

  TriangleMesh& hg;

  VertexIndexMap vertex_id_pmap;
  HalfedgeIndexMap hedge_id_pmap;
  TriangleMeshPointPMap hg_point_pmap;

  std::vector<double> edge_squared_lengths;

  class EdgeCompareFunctor
  {
  public:
      std::vector<double> edge_squared_lengths;

      EdgeCompareFunctor(std::vector<double>& squared_lengths)
      {
          edge_squared_lengths = squared_lengths;
      }

      bool operator() (const int& e0, const int& e1) const
      {
          double p0 = edge_squared_lengths[e0];
          double p1 = edge_squared_lengths[e1];
          return (p0 == p1) ? (e0 < e1) : (p0 < p1);
      }
  };

  // type of priority queue for edges
  typedef std::set<int, EdgeCompareFunctor>                                       Edge_queue;

// Public methods
public:
  Curve_skeleton(TriangleMesh& hg,
                 VertexIndexMap vertex_id_pmap,
                 HalfedgeIndexMap hedge_id_pmap,
                 TriangleMeshPointPMap hg_point_pmap)
    : hg(hg)
    , vertex_id_pmap(vertex_id_pmap)
    , hedge_id_pmap(hedge_id_pmap)
    , hg_point_pmap(hg_point_pmap)
  {}

  // Extracting the skeleton to a boost::graph data structure.
  template <class Graph>
  void extract_skeleton(Graph& curve)
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor                  vertex_desc;
    typedef typename boost::graph_traits<Graph>::edge_descriptor                    edge_desc;

    init();
    collapse();

    // orig_vertex_id maps the new id for a vertex to its original id.
    // new_vertex_id only contains id for vertices that have not been deleted.
    std::vector<int> new_vertex_id;
    std::vector<int> orig_vertex_id;
    new_vertex_id.clear();
    new_vertex_id.resize(vertex_to_edge.size(), -1);
    orig_vertex_id.clear();
    orig_vertex_id.resize(vertex_to_edge.size(), -1);

    int id = 0;
    for (size_t i = 0; i < is_vertex_deleted.size(); ++i)
    {
      if (!is_vertex_deleted[i])
      {
        orig_vertex_id[id] = static_cast<int>(i);
        new_vertex_id[i] = id++;
      }
    }

    // Mapping a skeleton vertex id to its descriptor.
    std::vector<vertex_desc> id_to_vd;
    id_to_vd.clear();
    id_to_vd.resize(id);

    for (int i = 0; i < id; ++i)
      id_to_vd[i] = boost::add_vertex(curve);

    for (int i = 0; i < id; ++i)
    {
      int orig_id = orig_vertex_id[i];
      vertex_desc vd = id_to_vd[i];

      /// code that is not working
      for(int vid : record[orig_id])
      {
        vertex_descriptor ovd = id_to_descriptor[vid];
        curve[vd].vertices.insert(
          curve[vd].vertices.end(),
          ovd->vertices.begin(),
          ovd->vertices.end()
        );
      }
    }

    for (size_t i = 0; i < is_edge_deleted.size(); ++i)
    {
      if (!is_edge_deleted[i])
      {
        int p1 = edge_to_vertex[i][0];
        int p2 = edge_to_vertex[i][1];
        int p1_id = new_vertex_id[p1];
        int p2_id = new_vertex_id[p2];
        vertex_desc p1_vd = id_to_vd[p1_id];
        vertex_desc p2_vd = id_to_vd[p2_id];

        bool exist;
        edge_desc edge;
        boost::tie(edge, exist) = boost::edge(p1_vd, p2_vd, curve);
        if (!exist)
        {
          boost::add_edge(p1_vd, p2_vd, curve);
        }
      }
    }

    for(vertex_descriptor vd : vertices(hg))
    {
      int id = static_cast<int>(get(vertex_id_pmap, vd));
      int new_id = new_vertex_id[id];
      if (new_id == -1)
      {
        continue;
      }

      // move to the centroid
      Point pos = Point(0, 0, 0);
      for (size_t i = 0; i < record[id].size(); ++i)
      {
        vertex_descriptor vd = id_to_descriptor[record[id][i]];
        Point pv = get(hg_point_pmap, vd);
        pos = Point(pos.x() + pv.x(), pos.y() + pv.y(), pos.z() + pv.z());
      }
      double num = static_cast<double>(record[id].size());
      curve[id_to_vd[new_id]].point = Point(pos.x() / num, pos.y() / num, pos.z() / num);
    }
  }

// Private methods
private:
  void init()
  {
    MCFSKEL_DEBUG( std::cerr <<"init" << std::endl; )

    int nb_edges = static_cast<int>(num_edges(hg));
    int num_faces = static_cast<int>(hg.size_of_facets());
    int nb_vertices = static_cast<int>(num_vertices(hg));
    edge_to_face.resize(nb_edges);
    edge_to_vertex.resize(nb_edges);
    vertex_to_edge.resize(nb_vertices);
    face_to_edge.resize(num_faces);

    is_vertex_deleted.resize(nb_vertices, false);
    is_edge_deleted.resize(nb_edges, false);
    is_face_deleted.resize(num_faces, false);

    record.resize(nb_vertices);
    for (size_t i = 0; i < record.size(); ++i)
    {
      record[i].push_back(static_cast<int>(i));
    }

    id_to_descriptor.resize(nb_vertices);
    edge_squared_lengths.resize(nb_edges);

    // assign vertex id
    surface_vertex_id.resize(nb_vertices);
    int idx = 0;
    for(vertex_descriptor vd : vertices(hg))
    {
      surface_vertex_id[idx] = static_cast<int>(get(vertex_id_pmap, vd));
      put(vertex_id_pmap, vd, idx++);
    }

    // assign edge id
    // the two halfedges representing the same edge get the same id
    idx = 0;
    for(edge_descriptor ed : edges(hg))
    {
      halfedge_descriptor hd = halfedge(ed, hg);
      put(hedge_id_pmap, hd, idx);
      halfedge_descriptor hd_opposite = opposite(hd,hg);
      put(hedge_id_pmap, hd_opposite, idx);

      // also cache the length of the edge
      vertex_descriptor v1 = target(hd,hg);
      vertex_descriptor v2 = source(hd,hg);
      Point source = get(hg_point_pmap, v1);
      Point target = get(hg_point_pmap, v2);
      edge_squared_lengths[idx] = squared_distance(source, target);

      ++idx;
    }

    // assign face id and compute edge-face connectivity
    int face_id = 0;
    for(face_descriptor fd : faces(hg))
    {
      for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd,hg), hg))
      {
        int id = static_cast<int>(get(hedge_id_pmap, hd));
        face_to_edge[face_id].push_back(id);
        edge_to_face[id].push_back(face_id);
      }
      ++face_id;
    }

    // compute vertex-edge connectivity
    for(vertex_descriptor vd : vertices(hg))
    {
      int vid = static_cast<int>(get(vertex_id_pmap, vd));
      for(edge_descriptor ed : in_edges(vd, hg))
      {
        halfedge_descriptor hd = halfedge(ed, hg);
        int eid = static_cast<int>(get(hedge_id_pmap, hd));
        vertex_to_edge[vid].push_back(eid);
        edge_to_vertex[eid].push_back(vid);
      }

      // save the vertex descriptor
      id_to_descriptor[vid] = vd;
    }
  }

  void init_queue(Edge_queue& queue)
  {
    // put all the edges into a priority queue
    // shorter edge has higher priority
    for(edge_descriptor ed : edges(hg))
    {
      int id = static_cast<int>(get(hedge_id_pmap, halfedge(ed, hg)));
      queue.insert(id);
    }
  }

  // iteratively collapse short edges until no edges have incident faces
  void collapse()
  {
    EdgeCompareFunctor edge_comparator(edge_squared_lengths);
    Edge_queue queue(edge_comparator);

    init_queue(queue);

    // start collapsing edges until all the edges have no incident faces
    while (!queue.empty())
    {
      int eid = *(queue.begin());
      queue.erase(queue.begin());

      // skip already deleted edges and edges with no face
      if (is_edge_deleted[eid]) continue;
      if (edge_to_face[eid].size() == 0) continue;

      // mark the incident faces as deleted
      remove_incident_faces(eid);

      // p1 to be deleted
      int p1 = edge_to_vertex[eid][0];
      int p2 = edge_to_vertex[eid][1];
      is_vertex_deleted[p1] = true;

      // merge vertices collapsed on p1 to p2
      update_record(p1, p2);

      // delete the edge from p1 and p2's incident edges
      delete_edge(p1, p2, eid);

      // add the incident edges of p1 to p2
      add_edge(queue, p1, p2);

      // remove duplicate edges
      std::vector<int> vertex_to_edge_p2(vertex_to_edge[p2]);
      for (size_t i = 0; i < vertex_to_edge_p2.size(); ++i)
      {
        // ei to be removed
        int ei = vertex_to_edge_p2[i];
        for (size_t j = i + 1; j < vertex_to_edge_p2.size(); ++j)
        {
          int ej = vertex_to_edge_p2[j];
          if (is_same_edge(ei, ej) || is_edge_deleted[ei])
          {
            // look for ei from p2's incident edges
            bool found;
            int ind;
            boost::tie(found, ind) = find_edge(vertex_to_edge[p2], ei);
            if (!found)
            {
              continue;
            }

            // migrate faces from ei to ej
            move_face(ei, ej);

            // finally remove ei from p2
            remove_edge(p2, ei, ind);
            break;
          }
        }
      }
    }

    // for debugging purpose
    MCFSKEL_INFO(print_stat();)
    MCFSKEL_INFO(check_edge();)
  }

  void add_edge(Edge_queue& queue, int p1, int p2)
  {
    for (size_t i = 0; i < vertex_to_edge[p1].size(); ++i)
    {
      int edge = vertex_to_edge[p1][i];
      if (is_edge_deleted[edge])
      {
        continue;
      }
      vertex_to_edge[p2].push_back(edge);

      // after change the incident vertex of an edge,
      // we need to update the length of the edge
      update_edge_length(queue, edge, p1, p2);

      // change the incident vertex to p2
      for (size_t j = 0; j < edge_to_vertex[edge].size(); ++j)
      {
        if (edge_to_vertex[edge][j] == p1)
        {
          edge_to_vertex[edge][j] = p2;
        }
      }
    }
  }

  void remove_incident_faces(int eid)
  {
    for(int fid : edge_to_face[eid])
    {
      is_face_deleted[fid] = true;
      // remove the face from the container of the other incident edges
      for (size_t j = 0; j < face_to_edge[fid].size(); ++j)
      {
        int e = face_to_edge[fid][j];
        if (e==eid) continue;
        for (size_t k = 0; k < edge_to_face[e].size(); ++k)
        {
          if (edge_to_face[e][k] == fid)
          {
            edge_to_face[e].erase(edge_to_face[e].begin() + k);
            break;
          }
        }
      }
    }
    edge_to_face[eid].clear();
  }

  void delete_edge(int p1, int p2, int eid)
  {
    for (size_t i = 0; i < vertex_to_edge[p1].size(); ++i)
    {
      if (vertex_to_edge[p1][i] == eid)
      {
        vertex_to_edge[p1].erase(vertex_to_edge[p1].begin() + i);
        break;
      }
    }
    for (size_t i = 0; i < vertex_to_edge[p2].size(); ++i)
    {
      if (vertex_to_edge[p2][i] == eid)
      {
        vertex_to_edge[p2].erase(vertex_to_edge[p2].begin() + i);
        break;
      }
    }
    is_edge_deleted[eid] = true;
  }

  bool is_same_edge(int ei, int ej)
  {
    if (edge_to_vertex[ei][0] == edge_to_vertex[ej][0]
     && edge_to_vertex[ei][1] == edge_to_vertex[ej][1])
    {
      return true;
    }
    if (edge_to_vertex[ei][1] == edge_to_vertex[ej][0]
     && edge_to_vertex[ei][0] == edge_to_vertex[ej][1])
    {
      return true;
    }
    return false;
  }

  std::pair<bool, int> find_edge(std::vector<int>& edges, int eid)
  {
    for (size_t i = 0; i < edges.size(); ++i)
    {
      if (eid == edges[i])
      {
        return std::make_pair(true, static_cast<int>(i));
      }
    }
    return std::make_pair(false, -1);
  }

  void move_face(int ei, int ej)
  {
    for (size_t i = 0; i < edge_to_face[ei].size(); ++i)
    {
      int fid = edge_to_face[ei][i];
      if (!is_face_deleted[fid])
      {
        if (std::find(edge_to_face[ej].begin(),
                      edge_to_face[ej].end(),
                      fid)
            == edge_to_face[ej].end())
        {
          edge_to_face[ej].push_back(fid);
          for (size_t j = 0; j < face_to_edge[fid].size(); ++j)
          {
            if (face_to_edge[fid][j] == ei)
            {
              face_to_edge[fid][j] = ej;
              break;
            }
          }
        }
      }
    }
  }

  void remove_edge(int v, int e, int ind)
  {
    vertex_to_edge[v].erase(vertex_to_edge[v].begin() + ind);
    // and also remove ei from the other end point
    for (size_t i = 0; i < edge_to_vertex[e].size(); ++i)
    {
      int vid = edge_to_vertex[e][i];
      if (vid != v)
      {
        for (size_t j = 0; j < vertex_to_edge[vid].size(); ++j)
        {
          if (vertex_to_edge[vid][j] == e)
          {
            vertex_to_edge[vid].erase(vertex_to_edge[vid].begin() + j);
            break;
          }
        }
      }
    }
    is_edge_deleted[e] = true;
  }

  void update_record(int p1, int p2)
  {
    for (size_t i = 0; i < record[p1].size(); ++i)
    {
      record[p2].push_back(record[p1][i]);
    }
    record[p1].clear();
  }

  void update_edge_length(Edge_queue& queue, int eid, int p1, int p2)
  {
    int vid1 = edge_to_vertex[eid][0];
    int vid2 = edge_to_vertex[eid][1];
    if (vid1 == p1)
    {
      vid1 = p2;
    }
    else
    {
      vid2 = p2;
    }
    if (queue.find(eid) != queue.end())
    {
      vertex_descriptor v1 = id_to_descriptor[vid1];
      vertex_descriptor v2 = id_to_descriptor[vid2];

      Point source = get(hg_point_pmap, v1);
      Point target = get(hg_point_pmap, v2);
      double new_len = squared_distance(source, target);

      edge_squared_lengths[eid] = new_len;
      queue.insert(eid);
    }
  }

  void check_edge()
  {
    for(halfedge_descriptor hd : halfedges(hg))
    {
      int id = get(hedge_id_pmap, hd);
      if (!is_edge_deleted[id])
      {
        if (edge_to_face[id].size() > 0)
        {
          std::cerr << "edge should not have faces " << edge_to_face[id].size() << "\n";
        }
      }
    }
  }

  void print_stat()
  {
    int cnt = 0;
    for (size_t i = 0; i < is_vertex_deleted.size(); ++i)
    {
      if (!is_vertex_deleted[i])
      {
        ++cnt;
      }
    }
    std::cerr << "num of vertices " << cnt << "\n";

    cnt = 0;
    for (size_t i = 0; i < is_edge_deleted.size(); ++i)
    {
      if (!is_edge_deleted[i])
      {
        ++cnt;
      }
    }
    std::cerr << "num of edges " << cnt << "\n";

    cnt = 0;
    for (size_t i = 0; i < is_face_deleted.size(); ++i)
    {
      if (!is_face_deleted[i])
      {
        ++cnt;
      }
    }
    std::cerr << "num of faces " << cnt << "\n";
  }
};

} // namespace internal
} // namespace CGAL

/// @endcond

#endif // CGAL_MCFSKEL_CURVE_SKELETON_H
