#ifndef CGAL_SHORTEST_NONCONTRACTIBLE_CYCLE_H
#define CGAL_SHORTEST_NONCONTRACTIBLE_CYCLE_H

#include <queue>
#include <CGAL/Generalized_map.h>

namespace CGAL {

template <class GeneralizedMap>
class Shortest_noncontractible_cycle {

public:

  using Gmap_origin = GeneralizedMap;

  struct Attributes {
    template <class GMap>
    struct Dart_wrapper {
      using Vertex_attribute = CGAL::Cell_attribute<GMap, int>;
      using Edge_attribute = void;
      using Face_attribute = CGAL::Cell_attribute<GMap, void>;
      using Attributes = CGAL::cpp11::tuple<Vertex_attribute, Edge_attribute, Face_attribute>;
    };
  };

  using Gmap = CGAL::Generalized_map<2, Attributes>;
  using Dart_handle = typename Gmap::Dart_handle;
  using size_type = typename Gmap::size_type;
  using Dart_const_handle = typename Gmap::Dart_const_handle;
  using Dart_container = std::vector<Dart_handle>;
  using Path = std::vector<typename Gmap_origin::Dart_const_handle>; // Consider: CGAL::Path_on_surface<Gmap>;
  using Distance_type = int;
  
  Shortest_noncontractible_cycle(const Gmap_origin& gmap) {
    m_gmap.copy(gmap, &m_origin_to_copy, &m_copy_to_origin);
    for (auto it = m_gmap.template one_dart_per_cell<2>().begin(), itend = m_gmap.template one_dart_per_cell<2>().end(); it != itend; ++it) {
      m_gmap.template set_attribute<2>(it, m_gmap.template create_attribute<2>());
    }
    m_gmap.template close<2>();
    m_gmap.display_characteristics(std::cout);
    std::cout << '\n';
  }
  
  Path find_cycle(typename Gmap_origin::Dart_const_handle root_vertex) {
    Dart_handle root = m_origin_to_copy[root_vertex];
    Dart_container spanning_tree, noncon_edges;
    std::vector<Distance_type> distance_from_root;
    std::vector<int> trace_index;

    BFS(root, spanning_tree, distance_from_root, trace_index);
    std::cerr << "Done BFS. spanning_tree.size() = " << spanning_tree.size() << '\n';
    find_noncon_edges(spanning_tree, noncon_edges);
    std::cerr << "Done find_noncon_edges. noncon_edges.size() = " << noncon_edges.size() << '\n';

    Distance_type min_distance = -1;
    Dart_handle min_noncon_edge;
    int min_a = -1, min_b = -1;
    for (auto dh : noncon_edges) {
      Dart_handle a = dh, b = m_gmap.template alpha<0>(a);
      int index_a = m_gmap.template info<0>(a), index_b = m_gmap.template info<0>(b);
      Distance_type sum_distance = distance_from_root[index_a] + distance_from_root[index_b];
      if (min_distance < 0 || min_distance > sum_distance) {
        min_distance = sum_distance;
        min_noncon_edge = dh;
        min_a = index_a;
        min_b = index_b;
      }
    }

    Path cycle;
    if (min_distance < 0) return cycle; // empty cycle;
    // Trace back the path from `a` to root
    for (int ind = min_a - 1; ind != -1; ind = trace_index[ind])
      cycle.push_back(m_copy_to_origin[spanning_tree[ind]]);
      // If use Path_on_surface: cycle.push_back(m_gmap.template alpha<0>(spanning_tree[ind]));
    // Reverse: now it is the path from root to `a`
    std::reverse(cycle.begin(), cycle.end());
    cycle.push_back(m_copy_to_origin[min_noncon_edge]);
    // Trace back the path from `b` to root
    for (int ind = min_b - 1; ind != -1; ind = trace_index[ind])
      cycle.push_back(m_copy_to_origin[m_gmap.template alpha<0>(spanning_tree[ind])]);

    // CGAL_assertion(cycle.is_closed());

    return cycle;
  }

private:

  /// Create a spanning tree using BFS

  void BFS(Dart_handle root, Dart_container& spanning_tree,
           std::vector<int>& distance_from_root, std::vector<int>& trace_index) {
    // The first of pair is a dart of a 0-cell,
    // the second of pair is the index of the dart leading to this 0-cell in spanning_tree
    std::queue<std::pair<Dart_handle, int> > q;
    int vertex_index = 0;

    // spanning_tree is empty (so far), so the index is -1
    q.push(std::make_pair(root, -1));
    m_gmap.template set_attribute<0>(root, m_gmap.template create_attribute<0>());
    m_gmap.template info<0>(root) = vertex_index++;
    distance_from_root.push_back(0);
    // Note: distance_from_root will have n (= #0-cells) elements while spanning_tree only has n-1
    while (q.size()) {
      Dart_handle u = q.front().first;
      int ind = q.front().second;
      q.pop();
      // TODO: This iterator (one_dart_per_incident_cell) can have some overhead for time complexity
      //       comparing to a direct traversal using next/opposite operator. 
      for (auto it = m_gmap.template one_dart_per_incident_cell<1,0>(u).begin(), 
                itend = m_gmap.template one_dart_per_incident_cell<1,0>(u).end();
                it != itend; ++it) {
        Dart_handle v = m_gmap.template alpha<0>(it);
        if (m_gmap.template attribute<0>(v) == NULL) {
          m_gmap.template set_attribute<0>(v, m_gmap.template create_attribute<0>());
          m_gmap.template info<0>(v) = vertex_index++;
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
  }


  /// Check if a face (containing dh_face) has degree 1.
  /// If it does, let dh_only_edge = the edge separating it and its only adjacent face.

  bool is_degree_one_face(Dart_const_handle dh_face, Dart_const_handle& dh_only_edge, size_type edge_deleted) {
    int degree = 0;
    Dart_const_handle dh_edge = NULL;
    for (auto dh = m_gmap.template one_dart_per_incident_cell<1,2>(dh_face).begin(), dhend = m_gmap.template one_dart_per_incident_cell<1,2>(dh_face).end(); dh != dhend; ++dh) {
      if (!m_gmap.is_marked(dh, edge_deleted)) {
        if (dh_edge!=NULL) return false;
        dh_edge=dh;
      }
    }
    dh_only_edge = dh_edge;
    return true;
  }


  /// Find E_nc

  void find_noncon_edges(const Dart_container& spanning_tree, Dart_container& noncon_edges) {
    size_type face_deleted, edge_deleted;
    // Turn out that face_deleted is a redundant mark. Consider removing it.
    try {
      face_deleted = m_gmap.get_new_mark();
      edge_deleted = m_gmap.get_new_mark();
    } catch (typename Gmap::Exception_no_more_available_mark) {
      std::cerr << "No more free mark, exit." << std::endl;
      exit(-1);
    }
    for (auto it = m_gmap.template one_dart_per_cell<2>().begin(), itend = m_gmap.template one_dart_per_cell<2>().end(); it != itend; ++it) {
      if (m_gmap.template attribute<2>(it) == NULL) {
        for (auto dh = m_gmap.template one_dart_per_incident_cell<1,2>(it).begin(), dhend = m_gmap.template one_dart_per_incident_cell<1,2>(it).end(); dh != dhend; ++dh) {
          m_gmap.template mark_cell<1>(dh, edge_deleted);
        }
      }
      m_gmap.template mark_cell<2>(it, face_deleted);
    }
    for (auto dh : spanning_tree) {
      if (m_gmap.is_marked(dh, edge_deleted)) continue;
      m_gmap.template mark_cell<1>(dh, edge_deleted);
    }
    std::queue<Dart_const_handle> degree_one_faces;
    for (auto it = m_gmap.template one_dart_per_cell<2>().begin(), itend = m_gmap.template one_dart_per_cell<2>().end(); it != itend; ++it) {
      if (m_gmap.is_marked(it, face_deleted)) continue;
      Dart_const_handle dh_only_edge = it;
      if (is_degree_one_face(it, dh_only_edge, edge_deleted)) 
        degree_one_faces.push(dh_only_edge);
    }
    while (degree_one_faces.size()) {
      Dart_const_handle dh_face = degree_one_faces.front();
      degree_one_faces.pop();
      m_gmap.template mark_cell<2>(dh_face, face_deleted);
      m_gmap.template mark_cell<1>(dh_face, edge_deleted);
      Dart_const_handle dh_adj_face = m_gmap.template alpha<2>(dh_face);
      Dart_const_handle dh_only_edge = dh_adj_face;
      if (is_degree_one_face(dh_adj_face, dh_only_edge, edge_deleted))
        degree_one_faces.push(dh_only_edge);
    }
    for (auto it = m_gmap.template one_dart_per_cell<1>().begin(), itend = m_gmap.template one_dart_per_cell<1>().end(); it != itend; ++it) {
      if (!m_gmap.is_marked(it, edge_deleted)) {
        noncon_edges.push_back(it);
      }
    }
    m_gmap.free_mark(edge_deleted);
    m_gmap.free_mark(face_deleted);
  }

  Gmap m_gmap;
  boost::unordered_map<typename Gmap_origin::Dart_const_handle, Dart_handle> m_origin_to_copy;
  boost::unordered_map<Dart_handle, typename Gmap_origin::Dart_const_handle> m_copy_to_origin;
};

}

#endif