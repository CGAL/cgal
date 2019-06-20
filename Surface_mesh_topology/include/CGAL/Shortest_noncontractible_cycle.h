#ifndef CGAL_SHORTEST_NONCONTRACTIBLE_CYCLE_H
#define CGAL_SHORTEST_NONCONTRACTIBLE_CYCLE_H

#include <queue>
#include <CGAL/Generalized_map.h>

namespace CGAL {

template <class GeneralizedMap, class WeightFunctor = void>
class Shortest_noncontractible_cycle {

public:

  using Gmap_origin = GeneralizedMap;
  using Dart_const_handle_orig = typename Gmap_origin::Dart_const_handle;

  template <class T>
  struct Weight_functor { using Weight = T; };

  template <bool, class T, class>
  struct Weight_functor_selector : Weight_functor<T> { };

  template <class T, class F>
  struct Weight_functor_selector<false, T, F> : Weight_functor<F> { }; 
  
  struct Default_weight_functor {
    using Weight_t = unsigned int;
    Weight_t operator() (Dart_const_handle_orig) const { return 1; }
  };

  using Weight = typename Weight_functor_selector<std::is_same<WeightFunctor, void>::value,
                                                  Default_weight_functor,
                                                  WeightFunctor>::Weight;
  using Distance_type = typename Weight::Weight_t;

  struct Attributes {
    template <class GMap>
    struct Dart_wrapper {
      using Vertex_attribute = CGAL::Cell_attribute<GMap, int>;
      using Edge_attribute = CGAL::Cell_attribute<GMap, Distance_type>;
      using Face_attribute = CGAL::Cell_attribute<GMap, void>;
      using Attributes = CGAL::cpp11::tuple<Vertex_attribute, Edge_attribute, Face_attribute>;
    };
  };

  using Gmap = CGAL::Generalized_map<2, Attributes>;
  using Dart_handle = typename Gmap::Dart_handle;
  using size_type = typename Gmap::size_type;
  using Dart_container = std::vector<Dart_handle>;
  using Path = std::vector<Dart_const_handle_orig>; // Consider: CGAL::Path_on_surface<Gmap>;
  
  Shortest_noncontractible_cycle(const Gmap_origin& gmap, const Weight& wf = Weight())
  {
    m_gmap.copy(gmap, &m_origin_to_copy, &m_copy_to_origin);
    // m_gmap.display_characteristics(std::cerr);
    // std::cerr << '\n';
    // Initialize 2-attributes
    for (auto it = m_gmap.darts().begin(), itend = m_gmap.darts().end(); it != itend; ++it) {
      m_gmap.template set_attribute<2>(it, m_gmap.template create_attribute<2>());
    }
    // Remove all boundary by adding faces
    m_gmap.template close<2>();
    for (auto it = m_gmap.darts().begin(), itend = m_gmap.darts().end(); it != itend; ++it) {
      m_gmap.template set_attribute<1>(it, m_gmap.template create_attribute<1>());
      m_gmap.template set_attribute<0>(it, m_gmap.template create_attribute<0>());
    }
    // Initialize 1-attributes
    for (auto it = gmap.template one_dart_per_cell<1>().begin(), itend = gmap.template one_dart_per_cell<1>().end(); it != itend; ++it) {
      Dart_handle img_dart = m_origin_to_copy[it];
      m_gmap.template info<1>(img_dart) = wf(it);
    }
    // Initialize 0-attributes
    for (auto it = m_gmap.template one_dart_per_cell<0>().begin(), itend = m_gmap.template one_dart_per_cell<0>().end(); it != itend; ++it) {
      ++m_nb_of_vertices;
      // m_gmap.template info<0>(it) = -1;
    }
    m_spanning_tree.reserve(m_nb_of_vertices-1);
    m_distance_from_root.reserve(m_nb_of_vertices);
    m_trace_index.reserve(m_nb_of_vertices-1);
    // m_gmap.display_characteristics(std::cerr);
    // std::cerr << '\n';
  }
  
  bool find_cycle(Dart_const_handle_orig root_vertex, Path& cycle, Distance_type* length = NULL) {
    Dart_handle root = m_origin_to_copy[root_vertex];
    return this->find_cycle(root, cycle, length);
  }

  void edge_width(Path& cycle, Distance_type* length = NULL) {
    cycle.clear();
    bool first_check = true;
    Distance_type min_length = 0;
    for (auto it = m_gmap.template one_dart_per_cell<0>().begin(), itend = m_gmap.template one_dart_per_cell<0>().end(); it != itend; ++it) {
      Distance_type temp_length;
      if (first_check) {
        if (!find_cycle(it, cycle, &temp_length)) continue;
        min_length = temp_length;
        first_check = false;
      } else {
        if (find_cycle(it, cycle, &temp_length, &min_length))
          min_length = temp_length;
      }
    }
    if (length != NULL) *length = min_length;
  }

private:

  void find_spanning_tree(Dart_handle root, Dart_container& spanning_tree,
                          std::vector<Distance_type>& distance_from_root, std::vector<int>& trace_index) {
    if (std::is_same<WeightFunctor, void>::value)
      find_BFS_tree(root, spanning_tree, distance_from_root, trace_index);
    else
      find_Dijkstra_tree(root, spanning_tree, distance_from_root, trace_index);
  }

  struct Dijkstra_comparator {
    Dijkstra_comparator(const std::vector<Distance_type>& distance_from_root) : m_distance(distance_from_root) {}
    bool operator()(const int x, const int y) const { return m_distance[x] > m_distance[y]; }
  private:
    const std::vector<Distance_type>& m_distance;
  };

  /// Create a spanning tree using Dijkstra
  void find_Dijkstra_tree(Dart_handle root, Dart_container& spanning_tree,
                          std::vector<Distance_type>& distance_from_root, std::vector<int>& trace_index) {
    // Preparation
    Dijkstra_comparator dc (distance_from_root);
    std::priority_queue<int, std::vector<int>, Dijkstra_comparator> pq(dc);
    int vertex_index = 0;
    size_type vertex_visited;
    try {
      vertex_visited = m_gmap.get_new_mark();
    } catch (typename Gmap::Exception_no_more_available_mark) {
      std::cerr << "No more free mark, exit." << std::endl;
      exit(-1);
    }
    // Begin Dijkstra
    pq.push(0);
    m_gmap.template info<0>(root) = vertex_index;
    m_gmap.template mark_cell<0>(root, vertex_visited);
    distance_from_root.push_back(0);

    while (pq.size()) {
      int u_index = pq.top();
      pq.pop();
      Dart_handle u = (u_index == 0) ? root : m_gmap.template alpha<0>(spanning_tree[u_index - 1]);
      CGAL_assertion(u_index == m_gmap.template info<0>(u));
      bool first_run = true;
      for (auto it = u; first_run || it != u; it = m_gmap.template alpha<2,1>(it)) {
        first_run = false;
        Dart_handle v = m_gmap.template alpha<0>(it);
        Distance_type w = m_gmap.template info<1>(it);
        if (!m_gmap.is_marked(v, vertex_visited)) {
          int v_index = ++vertex_index;
          CGAL_assertion(v_index == distance_from_root.size());
          distance_from_root.push_back(distance_from_root[u_index] + w);
          spanning_tree.push_back(it);
          trace_index.push_back(u_index - 1);
          m_gmap.template info<0>(v) = v_index;
          m_gmap.template mark_cell<0>(v, vertex_visited);
          pq.push(v_index);
        } else {
          int v_index = m_gmap.template info<0>(v);
          if (distance_from_root[v_index] > distance_from_root[u_index] + w) {
            CGAL_assertion(v_index > 0);
            distance_from_root[v_index] = distance_from_root[u_index] + w;
            spanning_tree[v_index - 1] = it;
            trace_index[v_index - 1] = u_index - 1;
            pq.push(v_index);
          }
        }
      }
    }
    m_gmap.free_mark(vertex_visited);
    // End Dijkstra
  }

  /// Create a spanning tree using BFS
  void find_BFS_tree(Dart_handle root, Dart_container& spanning_tree,
                     std::vector<Distance_type>& distance_from_root, std::vector<int>& trace_index) {
    // Preparation
    std::queue<int> q;
    int vertex_index = 0;
    size_type vertex_visited;
    try {
      vertex_visited = m_gmap.get_new_mark();
    } catch (typename Gmap::Exception_no_more_available_mark) {
      std::cerr << "No more free mark, exit." << std::endl;
      exit(-1);
    }
    // Begin BFS
    q.push(0);
    m_gmap.template info<0>(root) = vertex_index;
    m_gmap.template mark_cell<0>(root, vertex_visited);
    distance_from_root.push_back(0);
    while (q.size()) {
      int u_index = q.front();
      q.pop();
      Dart_handle u = (u_index == 0) ? root : m_gmap.template alpha<0>(spanning_tree[u_index-1]);
      bool first_run = true;
      for (auto it = u; first_run || it != u; it = m_gmap.template alpha<2,1>(it)) {
        first_run = false;
        Dart_handle v = m_gmap.template alpha<0>(it);
        if (!m_gmap.is_marked(v, vertex_visited)) {
          int v_index = ++vertex_index;
          distance_from_root.push_back(1 + distance_from_root[u_index]);
          spanning_tree.push_back(it);
          // `it` will lead to v
          q.push(v_index);
          trace_index.push_back(u_index-1);
          m_gmap.template info<0>(v) = v_index;
          m_gmap.template mark_cell<0>(v, vertex_visited);
        }
      }
    }
    m_gmap.free_mark(vertex_visited);
    // End BFS
  }


  /// Check if there is only one unmarked edge around a face.
  /// If there is, let dh_adjacent_edge = the edge separating it and its only adjacent face.

  bool is_degree_one_face(Dart_handle dh_face, Dart_handle& dh_only_edge, size_type edge_deleted) {
    Dart_handle dh_edge = NULL;
    for (auto dh = m_gmap.template one_dart_per_incident_cell<1,2>(dh_face).begin(), dhend = m_gmap.template one_dart_per_incident_cell<1,2>(dh_face).end(); dh != dhend; ++dh) {
      if (!m_gmap.is_marked(dh, edge_deleted)) {
        if (dh_edge!=NULL) return false;
        dh_edge=dh;
      }
    }
    if (dh_edge == NULL) return false;
    dh_only_edge = dh_edge;
    return true;
  }


  /// Find E_nc

  void find_noncon_edges(const Dart_container& spanning_tree, Dart_container& noncon_edges) {
    noncon_edges.clear();
    size_type face_deleted, edge_deleted;
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
          if (m_gmap.is_marked(dh, edge_deleted)) continue;
          m_gmap.template mark_cell<1>(dh, edge_deleted);
        }
        m_gmap.template mark_cell<2>(it, face_deleted);
      }
    }
    for (auto dh : spanning_tree) {
      if (m_gmap.is_marked(dh, edge_deleted)) continue;
      m_gmap.template mark_cell<1>(dh, edge_deleted);
    }
    std::queue<Dart_handle> degree_one_faces;
    // Add to queue the degree-1 faces
    for (auto it = m_gmap.template one_dart_per_cell<2>().begin(), itend = m_gmap.template one_dart_per_cell<2>().end(); it != itend; ++it) {
      if (m_gmap.is_marked(it, face_deleted)) continue;
      Dart_handle dh_only_edge = NULL;
      if (is_degree_one_face(it, dh_only_edge, edge_deleted)) 
        degree_one_faces.push(dh_only_edge);
    }
    // Remove the degree-1 faces
    while (degree_one_faces.size()) {
      Dart_handle dh_face = degree_one_faces.front();
      degree_one_faces.pop();
      if (!m_gmap.is_marked(dh_face, face_deleted))
        m_gmap.template mark_cell<2>(dh_face, face_deleted);
      if (!m_gmap.is_marked(dh_face, edge_deleted))
        m_gmap.template mark_cell<1>(dh_face, edge_deleted);
      Dart_handle dh_adj_face = m_gmap.template alpha<2>(dh_face);
      if (m_gmap.is_marked(dh_adj_face, face_deleted)) continue;
      Dart_handle dh_only_edge = NULL;
      if (is_degree_one_face(dh_adj_face, dh_only_edge, edge_deleted))
        degree_one_faces.push(dh_only_edge);
    }
    for (auto it = m_gmap.template one_dart_per_cell<1>().begin(), itend = m_gmap.template one_dart_per_cell<1>().end(); it != itend; ++it) {
      if (m_gmap.template info<0>(it) >= 0 && !m_gmap.is_marked(it, edge_deleted)) {
        noncon_edges.push_back(it);
      }
    }
    m_gmap.free_mark(edge_deleted);
    m_gmap.free_mark(face_deleted);
  }

  void add_to_cycle(Dart_handle dh, Path& cycle) {
    CGAL_assertion(dh != NULL);
    if (m_gmap.template attribute<2>(dh) == NULL)
      dh = m_gmap.template alpha<2>(dh);
    CGAL_assertion(m_gmap.template attribute<2>(dh) != NULL);
    cycle.push_back(m_copy_to_origin[dh]);
  }

  bool find_cycle(Dart_handle root, Path& cycle, Distance_type* length = NULL, const Distance_type* max_length = NULL) {
    m_spanning_tree.clear();
    m_distance_from_root.clear();
    m_trace_index.clear();
    for (auto it = m_gmap.template one_dart_per_cell<0>().begin(), itend = m_gmap.template one_dart_per_cell<0>().end(); it != itend; ++it)
      m_gmap.template info<0>(it) = -1;
    find_spanning_tree(root, m_spanning_tree, m_distance_from_root, m_trace_index);
    find_noncon_edges(m_spanning_tree, m_noncon_edges);
    // std::cerr << "Done find_noncon_edges. noncon_edges.size() = " << m_noncon_edges.size() << '\n';

    bool first_check = true;
    Distance_type min_distance = 0;
    Dart_handle min_noncon_edge;
    int min_a = -1, min_b = -1;
    for (auto dh : m_noncon_edges) {
      Dart_handle a = dh, b = m_gmap.template alpha<0>(dh);
      int index_a = m_gmap.template info<0>(a), index_b = m_gmap.template info<0>(b);
      Distance_type sum_distance = m_distance_from_root[index_a] + m_distance_from_root[index_b] 
                                   + m_gmap.template info<1>(dh);
      if (first_check || min_distance > sum_distance) {
        min_distance = sum_distance;
        min_noncon_edge = dh;
        min_a = index_a;
        min_b = index_b;
        first_check = false;
      }
    }
    if (first_check) return false; // no cycle found
    if (length != NULL) *length = min_distance < 0 ? 0 : min_distance;
    if (max_length != NULL)
      if (min_distance >= *max_length) return false; // abort

    cycle.clear();
    // Trace back the path from `a` to root
    for (int ind = min_a - 1; ind != -1; ind = m_trace_index[ind])
      add_to_cycle(m_spanning_tree[ind], cycle);
    // Reverse: now it is the path from root to `a`
    std::reverse(cycle.begin(), cycle.end());
    add_to_cycle(min_noncon_edge, cycle);
    // Trace back the path from `b` to root
    for (int ind = min_b - 1; ind != -1; ind = m_trace_index[ind])
      add_to_cycle(m_gmap.template alpha<0>(m_spanning_tree[ind]), cycle);
    // CGAL_assertion(cycle.is_closed());

    return true;
  }

  Gmap m_gmap;
  boost::unordered_map<Dart_const_handle_orig, Dart_handle> m_origin_to_copy;
  boost::unordered_map<Dart_handle, Dart_const_handle_orig> m_copy_to_origin;
  unsigned int m_nb_of_vertices = 0;
  Dart_container m_spanning_tree, m_noncon_edges;
  std::vector<Distance_type> m_distance_from_root;
  std::vector<int> m_trace_index;
};

}

#endif