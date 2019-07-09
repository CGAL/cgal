#ifndef CGAL_SHORTEST_NONCONTRACTIBLE_CYCLE_H
#define CGAL_SHORTEST_NONCONTRACTIBLE_CYCLE_H

#include <queue>
#include <CGAL/Generalized_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>

namespace CGAL {
namespace Surface_mesh_topology {

template <class Mesh_, class Weight_ = void>
class Shortest_noncontractible_cycle {

public:

  using Mesh_original = Mesh_;

  template <class T>
  struct Weight_functor { using Weight = T; };

  template <bool, class T, class>
  struct Weight_functor_selector : Weight_functor<T> { };

  template <class T, class F>
  struct Weight_functor_selector<false, T, F> : Weight_functor<F> { }; 
  
  struct Default_weight_functor {
    using Weight_t = unsigned int;
    template <class T>
    Weight_t operator() (T) const { return 1; }
  };

  using Weight = typename Weight_functor_selector<std::is_same<Weight_, void>::value,
                                                  Default_weight_functor,
                                                  Weight_>::Weight;
  using Distance_type = typename Weight::Weight_t;

  struct Attributes {
    template <class GenericMap>
    struct Dart_wrapper {
      using Vertex_attribute = CGAL::Cell_attribute<GenericMap, int>;
      using Edge_attribute   = CGAL::Cell_attribute<GenericMap, Distance_type>;
      using Face_attribute   = CGAL::Cell_attribute<GenericMap, void>;
      using Attributes       = CGAL::cpp11::tuple<Vertex_attribute, Edge_attribute, Face_attribute>;
    };
  };

  struct SNC_for_generalized_map {
    using Generic_map                = CGAL::Generalized_map<2, Attributes>;
    using Dart_handle_original       = typename Mesh_original::Dart_handle;
    using Dart_const_handle_original = typename Mesh_original::Dart_const_handle;
    using Copy_to_origin_map         = boost::unordered_map<typename Generic_map::Dart_handle, 
                                                            Dart_handle_original>;
    using Origin_to_copy_map         = boost::unordered_map<Dart_handle_original, 
                                                            typename Generic_map::Dart_handle>;
    
    static typename Generic_map::Dart_handle opposite(Generic_map& amap, typename Generic_map::Dart_handle dh)
    {
      return amap.template alpha<2,0>(dh);
    }

    static void copy(Generic_map& target, Mesh_original& source,
                     Origin_to_copy_map& origin_to_copy, Copy_to_origin_map& copy_to_origin)
    {
      target.copy(source, &origin_to_copy, &copy_to_origin);
    }
    
    static void set_weights(Generic_map& target, Mesh_original& source, 
                            Origin_to_copy_map& origin_to_copy, const Weight& wf)
    {
      for (auto it = source.darts().begin(), itend = source.darts().end(); it != itend; ++it)
        target.template info<1>(origin_to_copy[it]) = wf(it);
    }
  };

  struct SNC_for_combinatorial_map {
    using Generic_map                = CGAL::Combinatorial_map<2, Attributes>;
    using Dart_handle_original       = typename Mesh_original::Dart_handle;
    using Dart_const_handle_original = typename Mesh_original::Dart_const_handle;
    using Copy_to_origin_map         = boost::unordered_map<typename Generic_map::Dart_handle, 
                                                            Dart_handle_original>;
    using Origin_to_copy_map         = boost::unordered_map<Dart_handle_original, 
                                                            typename Generic_map::Dart_handle>;
    
    static typename Generic_map::Dart_handle opposite(Generic_map& amap, typename Generic_map::Dart_handle dh)
    {
      return amap.opposite(dh);
    }

    static void copy(Generic_map& target, Mesh_original& source,
                     Origin_to_copy_map& origin_to_copy, Copy_to_origin_map& copy_to_origin)
    {
      target.copy(source, &origin_to_copy, &copy_to_origin);
    }

    static void set_weights(Generic_map& target, Mesh_original& source, 
                            Origin_to_copy_map& origin_to_copy, const Weight& wf)
    {
      for (auto it = source.darts().begin(), itend = source.darts().end(); it != itend; ++it)
        target.template info<1>(origin_to_copy[it]) = wf(it);
    }
  };

  template <class>
  struct Generic_map_selector {
    using Generic_map                = CGAL::Combinatorial_map<2, Attributes>;
    using Dart_handle_original       = typename boost::graph_traits<Mesh_original>::halfedge_descriptor;
    using Dart_const_handle_original = typename boost::graph_traits<Mesh_original>::halfedge_descriptor;
    using Copy_to_origin_map         = boost::unordered_map<typename Generic_map::Dart_handle, 
                                                            Dart_handle_original>;
    using Origin_to_copy_map         = boost::unordered_map<Dart_handle_original, 
                                                            typename Generic_map::Dart_handle>;
    
    static typename Generic_map::Dart_handle opposite(Generic_map& amap, typename Generic_map::Dart_handle dh)
    {
      return amap.opposite(dh);
    }

    static void copy(Generic_map& target, Mesh_original& source,
                     Origin_to_copy_map& origin_to_copy, Copy_to_origin_map& copy_to_origin)
    {
      target.import_from_halfedge_graph(source, &origin_to_copy, &copy_to_origin);
    }
    
    static void set_weights(Generic_map& target, Mesh_original& source, 
                            Origin_to_copy_map& origin_to_copy, const Weight& wf)
    {
      for (typename boost::graph_traits<Mesh_original>::halfedge_iterator it = source.halfedges_begin(), itend = source.halfedges_end(); it != itend; ++it)
        target.template info<1>(origin_to_copy[*it]) = wf(*it);
    }
  };

  template <unsigned int d, class Items, class Alloc, class Storage>
  struct Generic_map_selector< CGAL::Generalized_map<d, Items, Alloc, Storage> > : SNC_for_generalized_map {};

  template <unsigned int d, class Refs, class Items, class Alloc, class Storage>
  struct Generic_map_selector< CGAL::Generalized_map_base
                                <d, Refs, Items, Alloc, Storage> > : SNC_for_generalized_map {};

  template <unsigned int d, unsigned int d2, class Traits, class Items,
            class Alloc, template<unsigned int,class,class,class,class> class Map, class Storage>
  struct Generic_map_selector< CGAL::Linear_cell_complex_for_generalized_map
                              <d, d2, Traits, Items, Alloc, Map, Storage> > : SNC_for_generalized_map {};

  template <unsigned int d, class Items, class Alloc, class Storage>
  struct Generic_map_selector< CGAL::Combinatorial_map<d, Items, Alloc, Storage> > : SNC_for_combinatorial_map {};

  template <unsigned int d, class Refs, class Items, class Alloc, class Storage>
  struct Generic_map_selector< CGAL::Combinatorial_map_base
                                <d, Refs, Items, Alloc, Storage> > : SNC_for_combinatorial_map {};

  template <unsigned int d, unsigned int d2, class Traits, class Items,
            class Alloc, template<unsigned int,class,class,class,class> class Map, class Storage>
  struct Generic_map_selector< CGAL::Linear_cell_complex_for_combinatorial_map
                              <d, d2, Traits, Items, Alloc, Map, Storage> > : SNC_for_combinatorial_map {};

  using Gmap                       = typename Generic_map_selector<Mesh_original>::Generic_map;
  using Gmap_wrapper               = Generic_map_selector<Mesh_original>;
  using Dart_handle_original       = typename Gmap_wrapper::Dart_handle_original;
  using Dart_const_handle_original = typename Gmap_wrapper::Dart_const_handle_original;
  using Dart_handle                = typename Gmap::Dart_handle;
  using Attribute_handle_0         = typename Gmap::template Attribute_handle<0>::type;
  using Attribute_handle_1         = typename Gmap::template Attribute_handle<1>::type;
  using size_type                  = typename Gmap::size_type;
  using Dart_container             = std::vector<Dart_handle>;
  using Path                       = std::vector<Dart_handle_original>; // Consider: CGAL::Path_on_surface<Gmap>;
  
  Shortest_noncontractible_cycle(Mesh_original& gmap, const Weight& wf = Weight())
  {
    Gmap_wrapper::copy(m_gmap, gmap, m_origin_to_copy, m_copy_to_origin);
    // m_gmap.display_characteristics(std::cerr);
    // std::cerr << '\n';
    // Initialize 2-attributes
    for (auto it = m_gmap.darts().begin(), itend = m_gmap.darts().end(); it != itend; ++it) {
      if (m_gmap.template attribute<2>(it)==NULL)
      { m_gmap.template set_attribute<2>(it, m_gmap.template create_attribute<2>()); }
    }
    // Remove all boundary by adding faces
    m_gmap.template close<2>();
    for (auto it = m_gmap.darts().begin(), itend = m_gmap.darts().end(); it != itend; ++it) {
      if (m_gmap.template attribute<1>(it)==NULL)
      { m_gmap.template set_attribute<1>(it, m_gmap.template create_attribute<1>()); }
      if (m_gmap.template attribute<0>(it)==NULL)
      { m_gmap.template set_attribute<0>(it, m_gmap.template create_attribute<0>()); }
    }
    // Initialize 1-attributes
    Gmap_wrapper::set_weights(m_gmap, gmap, m_origin_to_copy, wf);
    // Count number of vertices
    int cnt = 0;
    for (auto it = m_gmap.template one_dart_per_cell<0>().begin(), 
              itend = m_gmap.template one_dart_per_cell<0>().end(); it != itend; ++it)
      ++cnt;
    for (auto it = m_gmap.template one_dart_per_cell<2>().begin(), 
              itend = m_gmap.template one_dart_per_cell<2>().end(); it != itend; ++it)
      m_face_list.push_back(it);
    m_spanning_tree.reserve(cnt - 1);
    m_distance_from_root.reserve(cnt);
    m_trace_index.reserve(cnt - 1);
    // m_gmap.display_characteristics(std::cerr);
    // std::cerr << '\n';
  }
  
  bool find_cycle(Dart_handle_original root_vertex, Path& cycle, Distance_type* length = NULL)
  {
    Dart_handle root = m_origin_to_copy[root_vertex];
    return this->find_cycle(root, cycle, length);
  }

  void edge_width(Path& cycle, Distance_type* length = NULL)
  {
    cycle.clear();
    bool first_check = true;
    Distance_type min_length = 0;
    for (Attribute_handle_0 att_it = m_gmap.template attributes<0>().begin(), att_itend = m_gmap.template attributes<0>().end(); att_it != att_itend; ++att_it) {
      Dart_handle it = att_it->dart();
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
                          std::vector<Distance_type>& distance_from_root, std::vector<int>& trace_index)
  {
    if (std::is_same<Weight_, void>::value)
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
                          std::vector<Distance_type>& distance_from_root, std::vector<int>& trace_index)
  {
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
      Dart_handle u = (u_index == 0) ? root : m_gmap.next(spanning_tree[u_index - 1]);
      CGAL_assertion(u_index == m_gmap.template info<0>(u));
      bool first_run = true;
      for (Dart_handle it = u; first_run || it != u; it = m_gmap.next(Gmap_wrapper::opposite(m_gmap, it))) {
        first_run = false;
        Dart_handle v = m_gmap.next(it);
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
        // std::cerr << '.';
      }
    }
    m_gmap.free_mark(vertex_visited);
    // End Dijkstra
  }

  /// Create a spanning tree using BFS
  void find_BFS_tree(Dart_handle root, Dart_container& spanning_tree,
                     std::vector<Distance_type>& distance_from_root, std::vector<int>& trace_index)
  {
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
      Dart_handle u = (u_index == 0) ? root : m_gmap.next(spanning_tree[u_index - 1]);
      CGAL_assertion(u_index == m_gmap.template info<0>(u));
      bool first_run = true;
      for (Dart_handle it = u; first_run || it != u; it = m_gmap.next(Gmap_wrapper::opposite(m_gmap, it))) {
        first_run = false;
        Dart_handle v = m_gmap.next(it);
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
  bool is_degree_one_face(Dart_handle dh_face, Dart_handle& dh_only_edge, size_type edge_deleted)
  {
    Dart_handle dh_edge = NULL;
    bool first_run = true;
    for (Dart_handle dh = dh_face; first_run || dh != dh_face; dh = m_gmap.next(dh)) {
      first_run = false;
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
  void find_noncon_edges(const Dart_container& spanning_tree, Dart_container& noncon_edges)
  {
    noncon_edges.clear();
    size_type face_deleted, edge_deleted;
    try {
      face_deleted = m_gmap.get_new_mark();
      edge_deleted = m_gmap.get_new_mark();
    } catch (typename Gmap::Exception_no_more_available_mark) {
      std::cerr << "No more free mark, exit." << std::endl;
      exit(-1);
    }
    for (Dart_handle dh_face : m_face_list) {
      if (m_gmap.template attribute<2>(dh_face) == NULL) {
        bool first_run = true;
        for (Dart_handle it = dh_face; first_run || it != dh_face; it = m_gmap.next(it)) {
          first_run = false;
          if (m_gmap.is_marked(it, edge_deleted)) continue;
          m_gmap.template mark_cell<1>(it, edge_deleted);
        }
        m_gmap.template mark_cell<2>(dh_face, face_deleted);
      }
      // std::cerr << '!';
    }
    for (auto dh : spanning_tree) {
      if (m_gmap.is_marked(dh, edge_deleted)) continue;
      m_gmap.template mark_cell<1>(dh, edge_deleted);
    }
    std::queue<Dart_handle> degree_one_faces;
    // Add to queue the degree-1 faces
    for (Dart_handle it : m_face_list) {
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
      Dart_handle dh_adj_face = Gmap_wrapper::opposite(m_gmap, dh_face);
      if (m_gmap.is_marked(dh_adj_face, face_deleted)) continue;
      Dart_handle dh_only_edge = NULL;
      if (is_degree_one_face(dh_adj_face, dh_only_edge, edge_deleted))
        degree_one_faces.push(dh_only_edge);
    }
    for (Attribute_handle_1 att_it = m_gmap.template attributes<1>().begin(), att_itend = m_gmap.template attributes<1>().end(); att_it != att_itend; ++att_it) {
      Dart_handle it = att_it->dart();
      if (m_gmap.template info<0>(it) >= 0 && !m_gmap.is_marked(it, edge_deleted)) {
        noncon_edges.push_back(it);
      }
    }
    m_gmap.free_mark(edge_deleted);
    m_gmap.free_mark(face_deleted);
  }

  void add_to_cycle(Dart_handle dh, Path& cycle)
  {
    CGAL_assertion(dh != NULL);
    if (m_gmap.template attribute<2>(dh) == NULL)
      dh = Gmap_wrapper::opposite(m_gmap, dh);
    CGAL_assertion(m_gmap.template attribute<2>(dh) != NULL);
    cycle.push_back(m_copy_to_origin[dh]);
  }

  bool find_cycle(Dart_handle root, Path& cycle, Distance_type* length = NULL, const Distance_type* max_length = NULL) {
    m_spanning_tree.clear();
    m_distance_from_root.clear();
    m_trace_index.clear();
    for (Attribute_handle_0 att_it = m_gmap.template attributes<0>().begin(), att_itend = m_gmap.template attributes<0>().end(); att_it != att_itend; ++att_it) {
      Dart_handle it = att_it->dart();
      m_gmap.template info<0>(it) = -1;
    }
    find_spanning_tree(root, m_spanning_tree, m_distance_from_root, m_trace_index);
    find_noncon_edges(m_spanning_tree, m_noncon_edges);
    // std::cerr << "Done find_noncon_edges. noncon_edges.size() = " << m_noncon_edges.size() << '\n';

    bool first_check = true;
    Distance_type min_distance = 0;
    Dart_handle min_noncon_edge;
    int min_a = -1, min_b = -1;
    for (auto dh : m_noncon_edges) {
      Dart_handle a = dh, b = m_gmap.next(dh);
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
      add_to_cycle(Gmap_wrapper::opposite(m_gmap, m_spanning_tree[ind]), cycle);
    // CGAL_assertion(cycle.is_closed());

    return true;
  }

  Gmap m_gmap;
  typename Gmap_wrapper::Origin_to_copy_map m_origin_to_copy;
  typename Gmap_wrapper::Copy_to_origin_map m_copy_to_origin;
  unsigned int m_nb_of_vertices = 0;
  Dart_container m_spanning_tree, m_noncon_edges, m_face_list;
  std::vector<Distance_type> m_distance_from_root;
  std::vector<int> m_trace_index;
};

} // namespace Surface_mesh_topology
} // namespace CGAL

#endif
