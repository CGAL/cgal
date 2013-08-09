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

#ifndef CGAL_MCFSKEL_CURVE_SKELETON_H
#define CGAL_MCFSKEL_CURVE_SKELETON_H

/// @cond CGAL_DOCUMENT_INTERNAL

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/graph/copy.hpp>

// For debugging macro
#include <CGAL/internal/Mean_curvature_skeleton/Debug.h>

namespace CGAL {
namespace internal {

template <class HalfedgeGraph, class Graph,
          class VertexIndexMap, class EdgeIndexMap>
class Curve_skeleton
{
// Public types
public:

  // Geometric types
  typedef typename HalfedgeGraph::Traits      Kernel;
  typedef typename Kernel::Vector_3           Vector;
  typedef typename Kernel::Point_3            Point;

  // Repeat HalfedgeGraph types
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor	        vertex_descriptor;
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_iterator            vertex_iterator;
  typedef typename boost::graph_traits<HalfedgeGraph>::edge_descriptor            edge_descriptor;
  typedef typename boost::graph_traits<HalfedgeGraph>::edge_iterator              edge_iterator;
  typedef typename boost::graph_traits<HalfedgeGraph>::in_edge_iterator           in_edge_iterator;
  typedef typename HalfedgeGraph::Facet_iterator                                  Facet_iterator;
  typedef typename HalfedgeGraph::Halfedge_around_facet_circulator                Halfedge_facet_circulator;

  // Repeat Graph types
  typedef typename boost::graph_traits<Graph>::edge_descriptor                    edge_desc;

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

  VertexIndexMap vertex_id_pmap;
  EdgeIndexMap edge_id_pmap;

  HalfedgeGraph& polyhedron;

  std::vector<double> edge_lengths;

  class EdgeCompareFunctor
  {
  public:
      std::vector<double> edge_lengths;

      EdgeCompareFunctor(std::vector<double>& lengths)
      {
          edge_lengths = lengths;
      }

      bool operator() (const int& e0, const int& e1) const
      {
          double p0 = edge_lengths[e0];
          double p1 = edge_lengths[e1];
          return (p0 == p1) ? (e0 < e1) : (p0 < p1);
      }
  };

  // type of priority queue for edges
  typedef std::set<int, EdgeCompareFunctor>                                    Edge_queue;

// Public methods
public:
  Curve_skeleton(HalfedgeGraph& polyhedron) : polyhedron(polyhedron)
  {
  }

  // extract the skeleton to a boost::graph data structure
  void extract_skeleton(Graph& graph, std::vector<Point>& points,
                        std::vector<std::vector<int> >& corr)
  {
    init();
    collapse();

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
        orig_vertex_id[id] = i;
        new_vertex_id[i] = id++;
      }
    }

    Graph curve(id);
    corr.clear();
    corr.resize(id);

    for (int i = 0; i < id; ++i)
    {
      int orig_id = orig_vertex_id[i];
      corr[i] = record[orig_id];
      for (size_t j = 0; j < corr[i].size(); ++j)
      {
        corr[i][j] = surface_vertex_id[corr[i][j]];
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

        bool exist;
        edge_desc edge;
        boost::tie(edge, exist) = boost::edge(p1_id, p2_id, curve);
        if (!exist)
        {
          boost::add_edge(p1_id, p2_id, curve);
        }
      }
    }

    vertex_iterator vb, ve;
    points.resize(id);
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
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
        Point pv = vd->point();
        pos = Point(pos.x() + pv.x(), pos.y() + pv.y(), pos.z() + pv.z());
      }
      double num = record[id].size();
      pos = Point(pos.x() / num, pos.y() / num, pos.z() / num);
      points[new_id] = pos;
    }
    boost::copy_graph(curve, graph);
  }

// Private methods
private:
  void init()
  {
    int num_edges = boost::num_edges(polyhedron) / 2;
    int num_faces = polyhedron.size_of_facets();
    int num_vertices = boost::num_vertices(polyhedron);
    edge_to_face.clear();
    edge_to_face.resize(num_edges);
    edge_to_vertex.clear();
    edge_to_vertex.resize(num_edges);
    vertex_to_edge.clear();
    vertex_to_edge.resize(num_vertices);
    face_to_edge.clear();
    face_to_edge.resize(num_faces);

    is_vertex_deleted.clear();
    is_vertex_deleted.resize(num_vertices, false);
    is_edge_deleted.clear();
    is_edge_deleted.resize(num_edges, false);
    is_face_deleted.clear();
    is_face_deleted.resize(num_faces, false);

    record.clear();
    record.resize(num_vertices);
    for (size_t i = 0; i < record.size(); ++i)
    {
      record[i].push_back(i);
    }

    id_to_descriptor.clear();
    id_to_descriptor.resize(num_vertices);

    edge_lengths.clear();
    edge_lengths.resize(num_edges);

    // assign vertex id
    surface_vertex_id.resize(num_vertices);
    vertex_iterator vb, ve;
    int idx = 0;
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      surface_vertex_id[idx] = boost::get(vertex_id_pmap, *vb);
      boost::put(vertex_id_pmap, *vb, idx++);
    }

    // assign edge id
    // the two halfedges representing the same edge get the same id
    edge_iterator eb, ee;
    idx = 0;
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, -1);
    }
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      edge_descriptor ed = *eb;
      int id = boost::get(edge_id_pmap, ed);
      if (id == -1)
      {
        boost::put(edge_id_pmap, ed, idx);
        edge_descriptor ed_opposite = ed->opposite();
        boost::put(edge_id_pmap, ed_opposite, idx);

        // also cache the length of the edge
        vertex_descriptor v1 = ed->vertex();
        vertex_descriptor v2 = ed->opposite()->vertex();
        Point source = v1->point();
        Point target = v2->point();
        edge_lengths[idx] = sqrtf(squared_distance(source, target));

        idx++;
      }
    }

    // assign face id and compute edge-face connectivity
    int face_id = 0;
    for (Facet_iterator i = polyhedron.facets_begin(); i != polyhedron.facets_end(); ++i)
    {
      Halfedge_facet_circulator j = i->facet_begin();
      // Facets in polyhedral surfaces are at least triangles.
      CGAL_assertion(CGAL::circulator_size(j) >= 3);
      do
      {
        int id = j->id();
        face_to_edge[face_id].push_back(id);
        edge_to_face[id].push_back(face_id);
      } while (++j != i->facet_begin());
      face_id++;
    }

    // compute vertex-edge connectivity
    for (boost::tie(vb, ve) = boost::vertices(polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor vd = *vb;
      int vid = boost::get(vertex_id_pmap, vd);
      in_edge_iterator e, e_end;
      for (boost::tie(e, e_end) = boost::in_edges(*vb, polyhedron); e != e_end; ++e)
      {
        edge_descriptor ed = *e;
        int eid = boost::get(edge_id_pmap, ed);
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
    edge_iterator eb, ee;
    std::vector<bool> is_edge_inserted;
    is_edge_inserted.clear();
    is_edge_inserted.resize(edge_to_face.size(), false);
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      edge_descriptor ed = *eb;
      int id = boost::get(edge_id_pmap, ed);

      if (is_edge_inserted[id])
      {
        continue;
      }

      queue.insert(id);
      is_edge_inserted[id] = true;
    }
  }

  // iteratively collapse short edges until no edges have incident faces
  void collapse()
  {
    EdgeCompareFunctor edge_comparator(edge_lengths);
    Edge_queue queue(edge_comparator);

    queue.clear();

    init_queue(queue);

    // start collapsing edges until all the edges have no incident faces
    while (!queue.empty())
    {
      int eid = *(queue.begin());
      queue.erase(queue.begin());

      if (is_edge_deleted[eid])
      {
        continue;
      }
      if (edge_to_face[eid].size() == 0)
      {
        continue;
      }

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
    std::vector<int> faces(edge_to_face[eid]);
    for (size_t i = 0; i < faces.size(); ++i)
    {
      int fid = faces[i];
      is_face_deleted[fid] = true;
      // remove face from the incident edges
      for (size_t j = 0; j < face_to_edge[fid].size(); ++j)
      {
        int e = face_to_edge[fid][j];
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
        return std::make_pair(true, i);
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
    vertex_descriptor v1 = id_to_descriptor[vid1];
    vertex_descriptor v2 = id_to_descriptor[vid2];

    Point source = v1->point();
    Point target = v2->point();
    double new_len = sqrtf(squared_distance(source, target));

    if (queue.find(eid) != queue.end())
    {
      edge_lengths[eid] = new_len;
      queue.insert(eid);
    }
  }

  void check_edge()
  {
    edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::edges(polyhedron); eb != ee; ++eb)
    {
      edge_descriptor ed = *eb;
      int id = boost::get(edge_id_pmap, ed);
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
        cnt++;
      }
    }
    std::cerr << "num of vertices " << cnt << "\n";

    cnt = 0;
    for (size_t i = 0; i < is_edge_deleted.size(); ++i)
    {
      if (!is_edge_deleted[i])
      {
        cnt++;
      }
    }
    std::cerr << "num of edges " << cnt << "\n";

    cnt = 0;
    for (size_t i = 0; i < is_face_deleted.size(); ++i)
    {
      if (!is_face_deleted[i])
      {
        cnt++;
      }
    }
    std::cerr << "num of faces " << cnt << "\n";
  }
};

} // namespace internal
} // namespace CGAL

/// @endcond

#endif // CGAL_MCFSKEL_CURVE_SKELETON_H
