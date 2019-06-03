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
  using Dart_container = std::vector<Dart_const_handle>;
  using Path = CGAL::Path_on_surface<Gmap>;
  using Distance_type = int;

  Shortest_noncontractible_cycle(Gmap& gmap) :
    m_gmap(gmap) {  }
  
  Path find_cycle(Dart_const_handle root) {
    Dart_container spanning_tree, noncon_edges;
    std::vector<Distance_type> distance_from_root;
    std::vector<size_type> index_marks;
    std::vector<int> trace_index;

    BFS(root, spanning_tree, distance_from_root, trace_index);
    find_noncon_edges(spanning_tree, noncon_edges);
    mark_vertices_with_indicies(spanning_tree, index_marks);

    Distance_type min_distance = -1;
    Dart_const_handle min_noncon_edge;
    int min_a = -1, min_b = -1;
    for (auto dh : noncon_edges) {
      Dart_const_handle a = dh, b = m_gmap.alpha<0>(a);
      int index_a = 0, index_b = 0;
      for (int i = 0; i < index_marks.size(); ++i) {
        if (m_gmap.is_marked(a, index_marks[i])) index_a ^= (1 << i);
        if (m_gmap.is_marked(b, index_marks[i])) index_b ^= (1 << i);
      }
      Distance_type sum_distance = distance_from_root[index_a] + distance_from_root[index_b];
      if (min_distance < 0 || min_distance > sum_distance) {
        min_distance = sum_distance;
        min_noncon_edge = dh;
        min_a = index_a;
        min_b = index_b;
      }
    }

    Path cycle(m_gmap);
    if (min_distance < 0) return cycle; // empty cycle;
    // Trace back the path from `a` to root
    int ind;
    do {
      ind = min_a;
      cycle.push_back(m_gmap.alpha<0>(spanning_tree[ind]));
      ind = trace_index[ind];
    } while (ind != -1);
    // Reverse: now it is the path from root to `a`
    cycle.reverse();
    cycle.push_back(min_noncon_edge);
    // Trace back the path from `b` to root
    do {
      ind = min_b;
      cycle.push_back(m_gmap.alpha<0>(spanning_tree[ind]));
      ind = trace_index[ind];
    } while (ind != -1);

    CGAL_assertion(cycle.is_closed());

    return cycle;
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

  void BFS(Dart_const_handle root, Dart_container& spanning_tree,
           std::vector<int>& distance_from_root, std::vector<int>& trace_index) {
    size_type vertex_visited;
    try {
      vertex_visited = m_gmap.get_new_mark();
    } catch (typename Gmap::Exception_no_more_available_mark) {
      std::cerr << "No more free mark, exit." << std::endl;
      exit(-1);
    }
    // The first of pair is a dart of a 0-cell,
    // the second of pair is the index of the dart leading to this 0-cell in spanning_tree
    std::queue<std::pair<Dart_const_handle, int> > q;

    // spanning_tree is empty (so far), so the index is -1
    q.push(std::make_pair(root, -1));
    mark_cell<0>(root, vertex_visited);
    distance_from_root.push_back(0);
    // Note: distance_from_root will have n (= #0-cells) elements while spanning_tree only has n-1
    while (q.size()) {
      Dart_const_handle u = q.front().first;
      int ind = q.front().second;
      q.pop();
      for (auto it = m_gmap.template one_dart_per_incident_cell<1,0>(u).begin(), 
                itend = m_gmap.template one_dart_per_incident_cell<1,0>(u).end();
                it != itend; ++it) {
        Dart_const_handle v = m_gmap.alpha<0>(it);
        if (!m_gmap.is_marked(v, vertex_visited)) {
          mark_cell<0>(v, vertex_visited);
          if (ind == -1) 
            distance_from_root.push_back(1);
          else
            distance_from_root.push_back(1 + distance_from_root[ind+1]);
          spanning_tree.push_back(it);
          // `it` will lead to v
          q.push(make_pair(v, spanning_tree.size() - 1));
          trace_index.push_back(ind);
          // v will later be marked with the binary code of spanning_tree.size()
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
    for (auto dh = m_gmap.template one_dart_per_incident_cell<1,2>(dh_face).begin(), dhend = m_gmap.template one_dart_per_incident_cell<1,2>(dh_face).end(); dh != dhend; ++dh) {
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

  void find_noncon_edges(const Dart_container& spanning_tree, Dart_container noncon_edges) {
    size_type face_deleted, edge_deleted;
    try {
      face_deleted = m_gmap.get_new_mark();
      edge_deleted = m_gmap.get_new_mark();
    } catch (typename Gmap::Exception_no_more_available_mark) {
      std::cerr << "No more free mark, exit." << std::endl;
      exit(-1);
    }
    std::queue<Dart_const_handle> degree_one_faces;
    for (auto dh : spanning_tree) {
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
        noncon_edges.push_back(it);
      }
    }
    m_gmap.free_mark(edge_deleted);
    m_gmap.free_mark(face_deleted);
  }

  void mark_vertices_with_indicies(const Dart_container& spanning_tree, std::vector<size_type>& index_marks) {
    int n = spanning_tree.size()+1, bit_count = 0;
    while (n) ++bit_count, n >> 1;
    index_marks.resize(bit_count);
    try {
      for (int i = 0; i < index_marks.size(); ++i)
        index_marks[i] = m_gmap.get_new_mark();
    } catch (typename Gmap::Exception_no_more_available_mark) {
      std::cerr << "No more free mark, exit." << std::endl;
      exit(-1);
    }
    n = 0;
    for (auto dh : spanning_tree) {
      ++n;
      Dart_const_handle v = m_gmap.alpha<0>(dh);
      for (int i = 0; i < index_marks.size(); ++i)
        if (n & (1 << i)) mark_cell<0>(v, index_marks[i]);
    }
  }

  Gmap m_gmap;
};

}

#endif