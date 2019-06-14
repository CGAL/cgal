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
  struct Weight_functor_selector {
    using Weight = T;
  };

  template <>
  struct Weight_functor_selector<void> {
    struct Weight {
      using Weight_t = unsigned int;
      Weight() {}
      Weight_t operator() (Dart_const_handle_orig) { return 1; }
    };
  };

  using Weight = typename Weight_functor_selector<WeightFunctor>::Weight;
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
    m_gmap.display_characteristics(std::cerr);
    std::cerr << '\n';
    // Initialize 2-attributes
    for (auto it = m_gmap.template one_dart_per_cell<2>().begin(), itend = m_gmap.template one_dart_per_cell<2>().end(); it != itend; ++it) {
      m_gmap.template set_attribute<2>(it, m_gmap.template create_attribute<2>());
    }
    // Remove all boundary by adding faces
    m_gmap.template close<2>();
    // Initialize 1-attributes
    for (auto it = gmap.template one_dart_per_cell<1>().begin(), itend = gmap.template one_dart_per_cell<1>().end(); it != itend; ++it) {
      Dart_handle img_dart = m_origin_to_copy[it];
      m_gmap.template set_attribute<1>(img_dart, m_gmap.template create_attribute<1>());
      m_gmap.template info<1>(img_dart) = wf(it);
    }
    // Initialize 0-attributes
    for (auto it = m_gmap.template one_dart_per_cell<0>().begin(), itend = m_gmap.template one_dart_per_cell<0>().end(); it != itend; ++it) {
      m_gmap.template set_attribute<0>(it, m_gmap.template create_attribute<0>());
      ++m_nb_of_vertices;
    }
    m_spanning_tree.resize(m_nb_of_vertices-1);
    m_distance_from_root.resize(m_nb_of_vertices);
    m_trace_index.resize(m_nb_of_vertices-1);
    m_gmap.display_characteristics(std::cerr);
    std::cerr << '\n';
  }
  
  bool find_cycle(Dart_const_handle_orig root_vertex, Path& cycle, Distance_type* length = NULL, const Distance_type* max_length = NULL) {
    Dart_handle root = m_origin_to_copy[root_vertex];
    return this->find_cycle(root, cycle, length, max_length);
  }

  void edge_width(Path& cycle, Distance_type* length = NULL) {
    cycle.clear();
    bool first_check = true;
    Distance_type min_length = 0;
    int cnt = 0;
    for (auto it = m_gmap.template one_dart_per_cell<0>().begin(), itend = m_gmap.template one_dart_per_cell<0>().end(); it != itend; ++it) {
      // if (cnt % 100 == 0) std::cerr << "Processing vertex #" << cnt + 1 << ": \n";
      Distance_type temp_length;
      if (first_check) {
        find_cycle(it, cycle, &temp_length);
        min_length = temp_length;
        first_check = false;
      } else {
        find_cycle(it, cycle, &temp_length, &min_length);
        if (temp_length < min_length) {
          min_length = temp_length;
        }
      }
      // if (cnt % 100 == 0) std::cerr << temp_length << '\n';
      ++cnt;
    }
    if (length != NULL) *length = min_length;
  }

private:

  struct Dijkstra_comparator {
    Dijkstra_comparator(const std::vector<Distance_type>& distance_from_root) : m_distance(distance_from_root) {}
    bool operator()(const int x, const int y) const { return m_distance[x] > m_distance[y]; }
  private:
    const std::vector<Distance_type>& m_distance;
  };

  /// Create a spanning tree using Dijkstra
  template <class>
  void find_spanning_tree(Dart_handle root, Dart_container& spanning_tree,
                          std::vector<Distance_type>& distance_from_root, std::vector<int>& trace_index) {
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
    CGAL_assertion(spanning_tree.size() == m_nb_of_vertices - 1);
    CGAL_assertion(distance_from_root.size() == m_nb_of_vertices);
    CGAL_assertion(trace_index.size() == m_nb_of_vertices - 1);

    pq.push(0);
    m_gmap.template info<0>(root) = vertex_index;
    m_gmap.template mark_cell<0>(root, vertex_visited);
    distance_from_root[vertex_index] = 0;
    // std::cerr << "Begin finding spanning tree\n";

    while (pq.size()) {
      int u_index = pq.top();
      // std::cerr << '(' << u_index << ' ';
      pq.pop();
      Dart_handle u = (u_index == 0) ? root : m_gmap.template alpha<0>(spanning_tree[u_index - 1]);
      CGAL_assertion(u_index == m_gmap.template info<0>(u));

      // for (auto it = m_gmap.template one_dart_per_incident_cell<1,0>(u).begin(), itend = m_gmap.template one_dart_per_incident_cell<1,0>(u).end(); it != itend; ++it) {
      for (auto it = m_gmap.template alpha<2,1>(u); it != u; it = m_gmap.template alpha<2,1>(it)) {
        Dart_handle v = m_gmap.template alpha<0>(it);
        // std::cerr << '.';
        Distance_type w = m_gmap.template info<1>(it);
        // std::cerr << '-';
        if (!m_gmap.is_marked(v, vertex_visited)) {
          int v_index = ++vertex_index;
          // std::cerr << v_index << ' ';
          distance_from_root[v_index] = distance_from_root[u_index] + w;
          spanning_tree[v_index - 1] = it;
          trace_index[v_index - 1] = u_index - 1;
          m_gmap.template info<0>(v) = v_index;
          m_gmap.template mark_cell<0>(v, vertex_visited);
          pq.push(v_index);
          // if (pq.size() > m_nb_of_vertices) std::cerr << "Too long queue. Please halt.\n";
        } else {
          int v_index = m_gmap.template info<0>(v);
          // std::cerr << v_index << "* ";
          if (distance_from_root[v_index] > distance_from_root[u_index] + w) {
            // std::cerr << "+ ";
            CGAL_assertion(v_index > 0);
            distance_from_root[v_index] = distance_from_root[u_index] + w;
            spanning_tree[v_index - 1] = it;
            trace_index[v_index - 1] = u_index - 1;
            pq.push(v_index);
            // if (pq.size() > m_nb_of_vertices) std::cerr << "Too long queue. Please halt.\n";
          }
        }
      }
      // std::cerr << ") ";
    }
    // std::cerr << "\n";
    m_gmap.free_mark(vertex_visited);
  }

  /// Create a spanning tree using BFS
  template <>
  void find_spanning_tree<void>(Dart_handle root, Dart_container& spanning_tree,
                                std::vector<Distance_type>& distance_from_root, std::vector<int>& trace_index) {
    std::queue<int> q;
    int vertex_index = 0;
    size_type vertex_visited;
    try {
      vertex_visited = m_gmap.get_new_mark();
    } catch (typename Gmap::Exception_no_more_available_mark) {
      std::cerr << "No more free mark, exit." << std::endl;
      exit(-1);
    }
    CGAL_assertion(spanning_tree.size() == m_nb_of_vertices - 1);
    CGAL_assertion(distance_from_root.size() == m_nb_of_vertices);
    CGAL_assertion(trace_index.size() == m_nb_of_vertices - 1);
    q.push(0);
    m_gmap.template info<0>(root) = vertex_index;
    m_gmap.template mark_cell<0>(root, vertex_visited);
    distance_from_root[vertex_index] = 0;
    // Note: distance_from_root will have n (= #0-cells) elements while spanning_tree only has n-1
    while (q.size()) {
      int u_index = q.front();
      q.pop();
      Dart_handle u = (u_index == 0) ? root : m_gmap.template alpha<0>(spanning_tree[u_index-1]);
      for (auto it = m_gmap.template alpha<2,1>(u); it != u; it = m_gmap.template alpha<2,1>(it)) {
        Dart_handle v = m_gmap.template alpha<0>(it);
        if (!m_gmap.is_marked(v, vertex_visited)) {
          int v_index = ++vertex_index;
          distance_from_root[v_index] = 1 + distance_from_root[u_index];
          spanning_tree[v_index - 1] = it;
          // `it` will lead to v
          q.push(v_index);
          trace_index[v_index-1] = u_index-1;
          m_gmap.template info<0>(v) = v_index;
          m_gmap.template mark_cell<0>(v, vertex_visited);
        }
      }
    }
    m_gmap.free_mark(vertex_visited);
  }


  /// Check if a face (containing dh_face) has degree 1.
  /// If it does, let dh_only_edge = the edge separating it and its only adjacent face.

  bool is_degree_one_face(Dart_handle dh_face, Dart_handle& dh_only_edge, size_type edge_deleted) {
    int degree = 0;
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
        m_gmap.template mark_cell<2>(it, face_deleted);
      }
    }
    for (auto dh : spanning_tree) {
      if (m_gmap.is_marked(dh, edge_deleted)) continue;
      m_gmap.template mark_cell<1>(dh, edge_deleted);
    }
    std::queue<Dart_handle> degree_one_faces;
    for (auto it = m_gmap.template one_dart_per_cell<2>().begin(), itend = m_gmap.template one_dart_per_cell<2>().end(); it != itend; ++it) {
      if (m_gmap.is_marked(it, face_deleted)) continue;
      Dart_handle dh_only_edge = it;
      if (is_degree_one_face(it, dh_only_edge, edge_deleted)) 
        degree_one_faces.push(dh_only_edge);
    }
    // std::cerr << "degree_one_faces.size() = " << degree_one_faces.size() << std::endl;
    while (degree_one_faces.size()) {
      Dart_handle dh_face = degree_one_faces.front();
      degree_one_faces.pop();
      m_gmap.template mark_cell<2>(dh_face, face_deleted);
      m_gmap.template mark_cell<1>(dh_face, edge_deleted);
      Dart_handle dh_adj_face = m_gmap.template alpha<2>(dh_face);
      Dart_handle dh_only_edge = dh_adj_face;
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

  bool find_cycle(Dart_handle root, Path& cycle, Distance_type* length = NULL, Distance_type* max_length = NULL) {
    find_spanning_tree<WeightFunctor>(root, m_spanning_tree, m_distance_from_root, m_trace_index);
    // std::cerr << "Done find_spanning_tree. spanning_tree.size() = " << m_spanning_tree.size() << '\n';
    find_noncon_edges(m_spanning_tree, m_noncon_edges);
    // std::cerr << "Done find_noncon_edges. noncon_edges.size() = " << m_noncon_edges.size() << '\n';

    bool first_check = true;
    Distance_type min_distance = 0;
    Dart_handle min_noncon_edge;
    int min_a = -1, min_b = -1;
    for (auto dh : m_noncon_edges) {
      Dart_handle a = dh, b = m_gmap.template alpha<0>(a);
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
      cycle.push_back(m_copy_to_origin[m_spanning_tree[ind]]);
      // If use Path_on_surface: cycle.push_back(m_gmap.template alpha<0>(m_spanning_tree[ind]));
    // Reverse: now it is the path from root to `a`
    std::reverse(cycle.begin(), cycle.end());
    cycle.push_back(m_copy_to_origin[min_noncon_edge]);
    // Trace back the path from `b` to root
    for (int ind = min_b - 1; ind != -1; ind = m_trace_index[ind])
      cycle.push_back(m_copy_to_origin[m_gmap.template alpha<0>(m_spanning_tree[ind])]);

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