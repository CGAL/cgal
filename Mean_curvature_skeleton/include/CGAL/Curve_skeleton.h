#ifndef CURVE_SKELETON_H
#define CURVE_SKELETON_H

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/graph/copy.hpp>

namespace CGAL {

template <class Polyhedron, class Graph,
          class PolyhedronVertexIndexMap, class PolyhedronEdgeIndexMap>
class Curve_skeleton
{
// Public types
public:

  // Geometric types
  typedef typename Polyhedron::Traits         Kernel;
  typedef typename Kernel::Vector_3           Vector;
  typedef typename Kernel::Point_3            Point;

  // Repeat Polyhedron types
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor	         vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::vertex_iterator            vertex_iterator;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor            edge_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::edge_iterator              edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::in_edge_iterator           in_edge_iterator;
  typedef typename Polyhedron::Facet_iterator                                  Facet_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator                Halfedge_facet_circulator;

// Data members
private:
  std::vector<std::vector<int> > edge_to_face;
  std::vector<std::vector<int> > edge_to_vertex;
  std::vector<std::vector<int> > vertex_to_edge;
  std::vector<std::vector<int> > face_to_edge;
  std::vector<bool> is_vertex_deleted;
  std::vector<bool> is_edge_deleted;
  std::vector<bool> is_face_deleted;

  // records the vertices collapsed to a given vertex
  std::vector<std::vector<int> > record;
  // vertex id mapped to vertex descriptor
  std::vector<vertex_descriptor> id_to_descriptor;

  PolyhedronVertexIndexMap vertex_id_pmap;
  PolyhedronEdgeIndexMap edge_id_pmap;

  Polyhedron *polyhedron;

  struct Edge
  {
    int id;
    double length;

    bool operator <(const Edge& rhs) const
    {
      return (length == rhs.length) ? (id < rhs.id) : length < rhs.length;
    }
  };

// Public methods
public:
  Curve_skeleton(Polyhedron* polyhedron) : polyhedron(polyhedron)
  {
  }

  // extract the skeleton to a boost::graph data structure
  void extract_skeleton(Graph& graph, std::vector<Point>& points)
  {
    init();
    collapse();

    std::vector<int> new_vertex_id;
    new_vertex_id.clear();
    new_vertex_id.resize(vertex_to_edge.size(), -1);

    int id = 0;
    for (size_t i = 0; i < is_vertex_deleted.size(); i++)
    {
      if (!is_vertex_deleted[i])
      {
        new_vertex_id[i] = id++;
      }
    }

    Graph curve(id);

    for (size_t i = 0; i < is_edge_deleted.size(); i++)
    {
      if (!is_edge_deleted[i])
      {
        int p1 = edge_to_vertex[i][0];
        int p2 = edge_to_vertex[i][1];
        int p1_id = new_vertex_id[p1];
        int p2_id = new_vertex_id[p2];
        boost::add_edge(p1_id, p2_id, curve);
      }
    }

    vertex_iterator vb, ve;
    points.resize(id);
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      int id = boost::get(vertex_id_pmap, *vb);
      int new_id = new_vertex_id[id];
      if (new_id == -1)
      {
        continue;
      }
      Point pos = Point(0, 0, 0);
      for (size_t i = 0; i < record[id].size(); i++)
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
    int num_edges = boost::num_edges(*polyhedron) / 2;
    int num_faces = polyhedron->size_of_facets();
    int num_vertices = boost::num_vertices(*polyhedron);
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
    for (size_t i = 0; i < record.size(); i++)
    {
      record[i].push_back(i);
    }
    id_to_descriptor.clear();
    id_to_descriptor.resize(num_vertices);

    // assign vertex id
    vertex_iterator vb, ve;
    int idx = 0;
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      boost::put(vertex_id_pmap, *vb, idx++);
    }
    std::cerr << "vertex num " << idx << "\n";

    // assign edge id
    // the two halfedges representing the same edge get the same id
    edge_iterator eb, ee;
    idx = 0;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, -1);
    }
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      edge_descriptor ed = *eb;
      int id = boost::get(edge_id_pmap, ed);
      if (id == -1)
      {
        boost::put(edge_id_pmap, ed, idx);
        edge_descriptor ed_opposite = ed->opposite();
        boost::put(edge_id_pmap, ed_opposite, idx);
        idx++;
      }
    }
    std::cerr << "edge num " << idx << "\n";

    // assign face id and compute edge-face connectivity
    int face_id = 0;
    for (Facet_iterator i = polyhedron->facets_begin(); i != polyhedron->facets_end(); ++i)
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
    std::cerr << "face num " << face_id << "\n";

    // compute vertex-edge connectivity
    for (boost::tie(vb, ve) = boost::vertices(*polyhedron); vb != ve; ++vb)
    {
      vertex_descriptor vd = *vb;
      int vid = boost::get(vertex_id_pmap, vd);
      in_edge_iterator e, e_end;
      for (boost::tie(e, e_end) = boost::in_edges(*vb, *polyhedron); e != e_end; e++)
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

  void init_queue(std::set<Edge>& queue)
  {
    // put all the edges into a priority queue
    // shorter edge has higher priority
    edge_iterator eb, ee;
    std::vector<bool> is_edge_inserted;
    is_edge_inserted.clear();
    is_edge_inserted.resize(edge_to_face.size(), false);
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
    {
      Edge edge;

      edge_descriptor ed = *eb;
      edge.id = boost::get(edge_id_pmap, ed);

      if (is_edge_inserted[edge.id])
      {
        continue;
      }
      vertex_descriptor v1 = ed->vertex();
      vertex_descriptor v2 = ed->opposite()->vertex();

      Point source = v1->point();
      Point target = v2->point();
      edge.length = sqrtf(squared_distance(source, target));

      queue.insert(edge);
      is_edge_inserted[edge.id] = true;
    }
  }

  void collapse()
  {
    std::set<Edge> queue;
    queue.clear();

    init_queue(queue);

    // start collapsing edges until all the edges have no incident faces
    while (!queue.empty())
    {
      Edge edge = *(queue.begin());
      queue.erase(queue.begin());

      int eid = edge.id;

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
      add_edge(p1, p2);

      // remove duplicate edges
      std::vector<int> vertex_to_edge_p2(vertex_to_edge[p2]);
      for (size_t i = 0; i < vertex_to_edge_p2.size(); i++)
      {
        // ei to be removed
        int ei = vertex_to_edge_p2[i];
        for (size_t j = i + 1; j < vertex_to_edge_p2.size(); j++)
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
//            std::cerr << "delete edge " << ei << "\n";
            break;
          }
        }
      }
    }

    // for debugging purpose
    std::cerr << "finish collapse\n";
    print_stat();
    check_edge();
  }

  void add_edge(int p1, int p2)
  {
    for (size_t i = 0; i < vertex_to_edge[p1].size(); i++)
    {
      int edge = vertex_to_edge[p1][i];
      if (is_edge_deleted[edge])
      {
        continue;
      }
      vertex_to_edge[p2].push_back(edge);
      // change the incident vertex to p2
      for (size_t j = 0; j < edge_to_vertex[edge].size(); j++)
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
    for (size_t i = 0; i < faces.size(); i++)
    {
      int fid = faces[i];
      is_face_deleted[fid] = true;
      // remove face from the incident edges
      for (size_t j = 0; j < face_to_edge[fid].size(); j++)
      {
        int e = face_to_edge[fid][j];
        for (size_t k = 0; k < edge_to_face[e].size(); k++)
        {
          if (edge_to_face[e][k] == fid)
          {
            edge_to_face[e].erase(edge_to_face[e].begin() + k);
            break;
          }
        }
      }
//      std::cerr << "delete face " << fid << "\n";
    }
  }

  void delete_edge(int p1, int p2, int eid)
  {
    for (size_t i = 0; i < vertex_to_edge[p1].size(); i++)
    {
      if (vertex_to_edge[p1][i] == eid)
      {
        vertex_to_edge[p1].erase(vertex_to_edge[p1].begin() + i);
        break;
      }
    }
    for (size_t i = 0; i < vertex_to_edge[p2].size(); i++)
    {
      if (vertex_to_edge[p2][i] == eid)
      {
        vertex_to_edge[p2].erase(vertex_to_edge[p2].begin() + i);
        break;
      }
    }
    is_edge_deleted[eid] = true;
//    std::cerr << "delete edge " << eid << "\n";
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
    for (size_t i = 0; i < edges.size(); i++)
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
    for (size_t i = 0; i < edge_to_face[ei].size(); i++)
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
          for (size_t j = 0; j < face_to_edge[fid].size(); j++)
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
    for (size_t i = 0; i < edge_to_vertex[e].size(); i++)
    {
      int vid = edge_to_vertex[e][i];
      if (vid != v)
      {
        for (size_t j = 0; j < vertex_to_edge[vid].size(); j++)
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
    for (size_t i = 0; i < record[p1].size(); i++)
    {
      record[p2].push_back(record[p1][i]);
    }
    record[p1].clear();
  }

  void check_edge()
  {
    edge_iterator eb, ee;
    for (boost::tie(eb, ee) = boost::edges(*polyhedron); eb != ee; ++eb)
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
    for (size_t i = 0; i < is_vertex_deleted.size(); i++)
    {
      if (!is_vertex_deleted[i])
      {
        cnt++;
      }
    }
    std::cerr << "num of vertices " << cnt << "\n";

    cnt = 0;
    for (size_t i = 0; i < is_edge_deleted.size(); i++)
    {
      if (!is_edge_deleted[i])
      {
        cnt++;
      }
    }
    std::cerr << "num of edges " << cnt << "\n";

    cnt = 0;
    for (size_t i = 0; i < is_face_deleted.size(); i++)
    {
      if (!is_face_deleted[i])
      {
        cnt++;
      }
    }
    std::cerr << "num of faces " << cnt << "\n";
  }
};

} // namespace CGAL

#endif // CURVE_SKELETON_H

