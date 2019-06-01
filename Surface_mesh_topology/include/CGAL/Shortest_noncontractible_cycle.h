#ifndef CGAL_SHORTEST_NONCONTRACTIBLE_CYCLE_H
#define CGAL_SHORTEST_NONCONTRACTIBLE_CYCLE_H

#include <queue>
#include <CGAL/Path_on_surface.h>

namespace CGAL {

template <class GeneralizedMap>
class Shortest_noncontractible_cycle {
public:
  using Gmap = GeneralizedMap;
  using Dart_handle = typename Gmap::Dart_handle;
  using size_type = typename Gmap::size_type;
  using Dart_const_handle = typename Gmap::Dart_const_handle;
  using Path = CGAL::Path_on_surface<Gmap>;

  Shortest_noncontractible_cycle(Gmap& gmap) :
    m_gmap(gmap) {  }
  
  Path find_cycle(Dart_const_handle root) {
    m_root = root;
    BFS();
    find_noncon_edges();
  }

private:

  /// Mark every dart belonging to a specific cell

  template <unsigned int i>
  void mark_cell(Dart_const_handle dh, size_type gmap_mask) {
    for (auto it = m_gmap.darts_of_cell<i>(dh).begin(), itend = m_gmap.darts_of_cell<i>(dh).end(); it != itend; ++it) {
      m_gmap.mark(dh, gmap_mask);
    }
  }


  /// Create a spanning tree using BFS

  void BFS() {
    size_type vertex_visited;
    try {
      vertex_visited = m_gmap.get_new_mark();
    } catch (typename Gmap::Exception_no_more_available_mark) {
      std::cerr << "No more free mark, exit." << std::endl;
      exit(-1);
    }
    std::queue<Dart_const_handle> q;
    q.push(m_root);
    mark_cell<0>(m_root, vertex_visited);
    m_spanning_tree.push_back(m_root);
    while (q.size()) {
      Dart_const_handle u = q.front();
      q.pop();
      for (auto it = m_gmap.template one_dart_per_incident_cell<1,0>(u).begin(), 
                itend = m_gmap.template one_dart_per_incident_cell<1,0>(u).end();
                it != itend; ++it) {
        Dart_const_handle v = m_gmap.alpha<0>(it);
        if (!m_gmap.is_marked(v, vertex_visited)) {
          mark_cell<0>(v, vertex_visited);
          q.push(v);
          m_spanning_tree.push_back(it);
        }
      }
    }
    m_gmap.free_mark(vertex_visited);
  }


  /// Check if a face (containing dh_face) has degree 1.
  /// If it does, let dh_only_edge = the edge separating it and its only adjacent face.

  bool is_degree_one_face(Dart_const_handle dh_face, Dart_const_handle& dh_only_edge, size_type edge_deleted) {
    int degree = 0;
    Dart_const_handle dh_edge = dh_face;
    for (auto dh = m_gmap.one_dart_per_incident_cell<1,2>(dh_face).begin(), dhend = m_gmap.one_dart_per_incident_cell<1,2>(dh_face).end(); dh != dhend; ++dh) {
      if (m_gmap.is_marked(dh, edge_deleted)) continue;
      dh_edge = dh;
      ++degree;
    }
    if (degree == 1) {
      dh_only_edge = dh_edge;
      return true;
    }
    return false;
  }


  /// Find E_nc

  void find_noncon_edges() {
    size_type face_deleted, edge_deleted;
    try {
      face_deleted = m_gmap.get_new_mark();
      edge_deleted = m_gmap.get_new_mark();
    } catch (typename Gmap::Exception_no_more_available_mark) {
      std::cerr << "No more free mark, exit." << std::endl;
      exit(-1);
    }
    std::queue<Dart_const_handle> degree_one_faces;
    for (auto dh : m_spanning_tree) {
      mark_cell<1>(dh, edge_deleted);
    }
    for (auto it = m_gmap.template one_dart_per_cell<2>().begin(), itend = m_gmap.template one_dart_per_cell<2>().end(); it != itend; ++it) {
      Dart_const_handle dh_only_edge = it;
      if (is_degree_one_face(it, dh_only_edge, edge_deleted)) 
        degree_one_faces.push_back(dh_only_edge);
    }
    while (degree_one_faces.size()) {
      Dart_const_handle dh_face = degree_one_faces.front();
      degree_one_faces.pop();
      mark_cell<2>(dh_face, face_deleted);
      mark_cell<1>(dh_face, edge_deleted);
      Dart_const_handle dh_adj_face = m_gmap.alpha<2>(dh_face);
      Dart_const_handle dh_only_edge = dh_adj_face;
      if (is_degree_one_face(dh_adj_face, dh_only_edge, edge_deleted)) 
        degree_one_faces.push_back(dh_only_edge);
    }
    for (auto it = m_gmap.template one_dart_per_cell<1>().begin(), itend = m_gmap.template one_dart_per_cell<1>().end(); it != itend; ++it) {
      if (!m_gmap.is_marked(it, edge_deleted)) {
        m_noncon_edges.push_back(it);
      }
    }
    m_gmap.free_mark(edge_deleted);
    m_gmap.free_mark(face_deleted);
  }

  Gmap m_gmap;
  Dart_const_handle m_root;
  std::vector<Dart_const_handle> m_spanning_tree, m_noncon_edges;
};

}

#endif