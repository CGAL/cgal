// Copyright (c) 2020 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Thien Hoang <thienvhoang99@gmail.com>
//
#ifndef CGAL_SHORTEST_NONCONTRACTIBLE_CYCLE_H
#define CGAL_SHORTEST_NONCONTRACTIBLE_CYCLE_H

#include <CGAL/license/Surface_mesh_topology.h>

#include <CGAL/Path_on_surface.h>
#include <CGAL/Timer.h>
#include <CGAL/Cell_attribute.h>
#include <CGAL/Surface_mesh_topology/internal/Generic_map_selector.h>
#include <queue>
#include <tuple>
#include <unordered_map>

namespace CGAL {
namespace Surface_mesh_topology {
namespace internal {

template <class Mesh_>
class Shortest_noncontractible_cycle
{
public:
  using Self=Shortest_noncontractible_cycle<Mesh_>;
  using Mesh=Mesh_;

  using Original_map_wrapper=internal::Generic_map_selector<Mesh>;
  using Original_dart_const_handle=typename Original_map_wrapper::Dart_handle_original; // TODO SOLVE const problem with copy

  using Local_map=typename Original_map_wrapper::Generic_map;
  using Dart_handle      =typename Local_map::Dart_handle;
  using Dart_const_handle=typename Local_map::Dart_const_handle;
  using size_type        = typename Local_map::size_type;

  using Dart_container=std::vector<Dart_handle>;
  using Path          =CGAL::Surface_mesh_topology::Path_on_surface<Mesh>;

  // Associations between original darts and their copy.
  using Origin_to_copy=boost::unordered_map<Original_dart_const_handle, Dart_handle>;
  using Copy_to_origin=boost::unordered_map<Dart_handle, Original_dart_const_handle>;

  struct Default_weight_functor
  {
    using Weight_t=unsigned int;
    template <class T>
    Weight_t operator() (T) const { return 1; }
  };

  /// @return the reduced map
  const Local_map& get_local_map() const
  { return m_local_map; }

  /// @return the reduced map
  Local_map& get_local_map()
  { return m_local_map; }

  Shortest_noncontractible_cycle(const Mesh& amesh, bool display_time=false) :
    m_cycle(amesh)
  {
    CGAL::Timer t;
    if (display_time)
    { t.start(); }
    
    Original_map_wrapper::copy(m_local_map, const_cast<Mesh&>(amesh),
                               m_origin_to_copy, m_copy_to_origin);

    // Initialize m_is_hole
    try
    {
      m_is_hole=get_local_map().get_new_mark();
    }
    catch (typename Local_map::Exception_no_more_available_mark)
    {
      std::cerr<<"No more free mark, exit."<<std::endl;
      exit(-1);
    }
    get_local_map().negate_mark(m_is_hole);
    // Remove all boundary by adding faces
    get_local_map().template close<2>();
    get_local_map().negate_mark(m_is_hole);

    for (auto it=get_local_map().darts().begin(), itend=get_local_map().darts().end(); it!=itend; ++it)
    {
      if (get_local_map().template attribute<0>(it)==nullptr)
      { get_local_map().template set_attribute<0>(it, get_local_map().template create_attribute<0>()); }
      // if (get_local_map().template attribute<1>(it)==NULL) // For debug purpose only
      // { get_local_map().template set_attribute<1>(it, get_local_map().template create_attribute<1>()); }
    }
    // std::cerr << '\n';
    for (auto it=get_local_map().template one_dart_per_cell<2>().begin(),
           itend=get_local_map().template one_dart_per_cell<2>().end(); it!=itend; ++it)
    { m_face_list.push_back(it); }
    for (auto it = get_local_map().template one_dart_per_cell<1>().begin(),
              itend = get_local_map().template one_dart_per_cell<1>().end(); it != itend; ++it)
    {
      // get_local_map().template info<1>(it) = m_edge_list.size(); // For debug purpose only
      m_edge_list.push_back(it);
    }
    // get_local_map().display_characteristics(std::cerr);
    // std::cerr << '\n';

    if (display_time)
    {
      t.stop();
      std::cout<<" [TIME] Shortest_noncontractible_cycle constructor: "<<t.time()<<" seconds."<<std::endl;
    }
  }

  ~Shortest_noncontractible_cycle()
  {
    // std::cerr << "Destructor...\n"; // For testing unique_ptr
    get_local_map().free_mark(m_is_hole);
  }
  
  template <class WeightFunctor=Default_weight_functor>
  Path compute_cycle(Original_dart_const_handle root_vertex,
                     typename WeightFunctor::Weight_t* length,
                     const WeightFunctor& wf,
                     bool display_time=false)
  {
    CGAL::Timer t;
    if (display_time)
    { t.start(); }

    m_cycle.clear();
    Dart_handle root = m_origin_to_copy[root_vertex];
    this->compute_cycle(root, m_cycle, length, NULL, wf);

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] compute_cycle: "<<t.time()<<" seconds."<<std::endl;
    }

    return m_cycle;
  }

  template <class WeightFunctor=Default_weight_functor>
  Path compute_cycle(Original_dart_const_handle root_vertex,
                     typename WeightFunctor::Weight_t* length,
                     bool display_time=false)
  { return compute_cycle(root_vertex, length, Default_weight_functor(), display_time); }

  template <class WeightFunctor=Default_weight_functor>
  Path compute_cycle(Original_dart_const_handle root_vertex,
                     bool display_time=false)
  { return compute_cycle(root_vertex, nullptr, display_time); }


  template <class WeightFunctor=Default_weight_functor>
  Path compute_edgewidth(typename WeightFunctor::Weight_t* length,
                         const WeightFunctor& wf,
                         bool display_time=false)
  {
    CGAL::Timer t;
    if (display_time)
    { t.start(); }

    m_cycle.clear();
    bool first_check = true;
    typename WeightFunctor::Weight_t min_length=0;
    for (auto it=get_local_map().template one_dart_per_cell<0>().begin(),
           itend=get_local_map().template one_dart_per_cell<0>().end(); it!=itend; ++it)
    {
      typename WeightFunctor::Weight_t temp_length;
      if (first_check)
      {
        if (!compute_cycle(it, m_cycle, &temp_length, NULL, wf)) continue;
        min_length = temp_length;
        first_check = false;
      }
      else
      {
        if (compute_cycle(it, m_cycle, &temp_length, &min_length, wf))
        { min_length = temp_length; }
      }
    }
    if (length!=nullptr) { *length=min_length; }

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] compute_edgewidth: "<<t.time()<<" seconds."<<std::endl;
    }

    return m_cycle;
  }

  template <class WeightFunctor=Default_weight_functor>
  Path compute_edgewidth(typename WeightFunctor::Weight_t* length,
                         bool display_time=false)
  { return compute_edgewidth(length, Default_weight_functor(), display_time); }
  
  template <class WeightFunctor=Default_weight_functor>
  Path compute_edgewidth(bool display_time=false)
  { return compute_edgewidth(nullptr, display_time); }
  
protected:
  template <class WeightFunctor, class Distance_>
  void compute_spanning_tree(Dart_handle root, Dart_container& spanning_tree,
                             std::vector<Distance_>& distance_from_root,
                             std::vector<int>& trace_index,
                             const WeightFunctor& wf = Default_weight_functor())
  {
    if (std::is_same<WeightFunctor, Default_weight_functor>::value)
    { compute_BFS_tree(root, spanning_tree, distance_from_root, trace_index, wf); }
    else
    { compute_Dijkstra_tree(root, spanning_tree, distance_from_root, trace_index, wf); }
  }

  template <class Distance_>
  struct Dijkstra_comparator
  {
    using Distance_type=Distance_;
    Dijkstra_comparator(const std::vector<Distance_type>& distance_from_root) :
      m_distance(distance_from_root) {}
    bool operator()(const int x, const int y) const
    { return m_distance[x]>m_distance[y]; }
  private:
    const std::vector<Distance_type>& m_distance;
  };

  /// Create a spanning tree using Dijkstra
  template <class WeightFunctor, class Distance_>
  void compute_Dijkstra_tree(Dart_handle root, Dart_container& spanning_tree,
                             std::vector<Distance_>& distance_from_root,
                             std::vector<int>& trace_index,
                             const WeightFunctor& wf=Default_weight_functor())
  {
    // Preparation
    Dijkstra_comparator<Distance_> dc(distance_from_root);
    std::priority_queue<int, std::vector<int>, Dijkstra_comparator<Distance_> > pq(dc);
    int vertex_index=0;
    size_type vertex_visited;
    try
    {
      vertex_visited=get_local_map().get_new_mark();
    }
    catch (typename Local_map::Exception_no_more_available_mark)
    {
      std::cerr<<"No more free mark, exit."<<std::endl;
      exit(-1);
    }
    // Begin Dijkstra
    pq.push(0);
    get_local_map().template info<0>(root) = vertex_index;
    get_local_map().template mark_cell<0>(root, vertex_visited);
    distance_from_root.push_back(0);

    while (pq.size())
    {
      int u_index = pq.top();
      pq.pop();
      Dart_handle u = (u_index == 0) ? root : get_local_map().next(spanning_tree[u_index - 1]);
      CGAL_assertion(u_index == get_local_map().template info<0>(u));
      bool first_run = true;
      for (Dart_handle it = u; first_run || it != u; it = get_local_map().next(get_local_map().opposite2(it)))
      {
        first_run = false;
        Dart_handle v = get_local_map().next(it);
        Distance_ w = wf(m_copy_to_origin[nonhole_dart_of_same_edge(it)]);
        if (!get_local_map().is_marked(v, vertex_visited))
        {
          int v_index = ++vertex_index;
          CGAL_assertion(v_index == distance_from_root.size());
          distance_from_root.push_back(distance_from_root[u_index] + w);
          spanning_tree.push_back(it);
          trace_index.push_back(u_index - 1);
          get_local_map().template info<0>(v) = v_index;
          get_local_map().template mark_cell<0>(v, vertex_visited);
          pq.push(v_index);
        }
        else
        {
          int v_index = get_local_map().template info<0>(v);
          if (distance_from_root[v_index] > distance_from_root[u_index] + w)
          {
            CGAL_assertion(v_index > 0);
            distance_from_root[v_index] = distance_from_root[u_index] + w;
            spanning_tree[v_index - 1] = it;
            trace_index[v_index - 1] = u_index - 1;
            pq.push(v_index);
          }
        }
      }
    }
    get_local_map().free_mark(vertex_visited);
    // End Dijkstra
  }

  /// Create a spanning tree using BFS
  template <class WeightFunctor, class Distance_>
  void compute_BFS_tree(Dart_handle root, Dart_container& spanning_tree,
                        std::vector<Distance_>& distance_from_root,
                        std::vector<int>& trace_index,
                        const WeightFunctor& wf = Default_weight_functor())
  {
    // Preparation
    std::queue<int> q;
    int vertex_index = 0;
    size_type vertex_visited;
    try
    { vertex_visited=get_local_map().get_new_mark(); }
    catch (typename Local_map::Exception_no_more_available_mark)
    {
      std::cerr << "No more free mark, exit." << std::endl;
      exit(-1);
    }
    // Begin BFS
    q.push(0);
    get_local_map().template info<0>(root)=vertex_index;
    get_local_map().template mark_cell<0>(root, vertex_visited);
    distance_from_root.push_back(0);
    while (q.size())
    {
      int u_index = q.front();
      q.pop();
      Dart_handle u = (u_index == 0) ? root : get_local_map().next(spanning_tree[u_index - 1]);
      CGAL_assertion(u_index == get_local_map().template info<0>(u));
      bool first_run = true;
      for (Dart_handle it = u; first_run || it != u; it = get_local_map().next(get_local_map().opposite2(it)))
      {
        first_run = false;
        Dart_handle v = get_local_map().next(it);
        if (!get_local_map().is_marked(v, vertex_visited))
        {
          int v_index = ++vertex_index;
          distance_from_root.push_back(1 + distance_from_root[u_index]);
          spanning_tree.push_back(it);
          // `it` will lead to v
          q.push(v_index);
          trace_index.push_back(u_index-1);
          get_local_map().template info<0>(v) = v_index;
          get_local_map().template mark_cell<0>(v, vertex_visited);
        }
      }
    }
    get_local_map().free_mark(vertex_visited);
    // End BFS
  }

  /// Check if there is only one unmarked edge around a face.
  /// If there is, let dh_adjacent_edge = the edge separating it and its only adjacent face.
  bool is_degree_one_face(Dart_handle dh_face, Dart_handle& dh_only_edge, size_type edge_deleted)
  {
    Dart_handle dh_edge = NULL;
    bool first_run = true;
    for (Dart_handle dh = dh_face; first_run || dh != dh_face; dh = get_local_map().next(dh))
    {
      first_run = false;
      if (!get_local_map().is_marked(dh, edge_deleted))
      {
        if (dh_edge!=NULL) return false;
        dh_edge=dh;
      }
    }
    if (dh_edge == NULL) return false;
    dh_only_edge = dh_edge;
    return true;
  }

  /// Find E_nc
  void compute_noncon_edges(const Dart_container& spanning_tree, Dart_container& noncon_edges)
  {
    noncon_edges.clear();
    size_type face_deleted, edge_deleted;
    try
    {
      face_deleted = get_local_map().get_new_mark();
      edge_deleted = get_local_map().get_new_mark();
    }
    catch (typename Local_map::Exception_no_more_available_mark)
    {
      std::cerr << "No more free mark, exit." << std::endl;
      exit(-1);
    }
    for (Dart_handle dh_face : m_face_list)
    {
      if (get_local_map().is_marked(dh_face, m_is_hole))
      { get_local_map().template mark_cell<2>(dh_face, face_deleted); }
    }
    
    for (auto dh : spanning_tree)
    {
      if (get_local_map().is_marked(dh, edge_deleted)) continue; // TODO negate
      get_local_map().template mark_cell<1>(dh, edge_deleted);
    }
    std::queue<Dart_handle> degree_one_faces;
    // Add to queue the degree-1 faces
    for (Dart_handle it : m_face_list)
    {
      if (get_local_map().is_marked(it, face_deleted)) continue; // TODO negate
      Dart_handle dh_only_edge = NULL;
      if (is_degree_one_face(it, dh_only_edge, edge_deleted)) 
      { degree_one_faces.push(dh_only_edge); }
    }
    // Remove the degree-1 faces
    while (degree_one_faces.size())
    {
      Dart_handle dh_face=degree_one_faces.front();
      degree_one_faces.pop();
      if (!get_local_map().is_marked(dh_face, face_deleted))
        get_local_map().template mark_cell<2>(dh_face, face_deleted);
      if (!get_local_map().is_marked(dh_face, edge_deleted))
        get_local_map().template mark_cell<1>(dh_face, edge_deleted);
      Dart_handle dh_adj_face=get_local_map().opposite2(dh_face);
      if (get_local_map().is_marked(dh_adj_face, face_deleted)) continue;
      Dart_handle dh_only_edge=NULL;
      if (is_degree_one_face(dh_adj_face, dh_only_edge, edge_deleted))
        degree_one_faces.push(dh_only_edge);
    }
    for (Dart_handle it : m_edge_list)
    {
      if (get_local_map().template info<0>(it)>=0 && !get_local_map().is_marked(it, edge_deleted))
      {
        noncon_edges.push_back(it);
      }
    }
    get_local_map().free_mark(edge_deleted);
    get_local_map().free_mark(face_deleted);
  }

  Dart_handle nonhole_dart_of_same_edge(Dart_handle dh)
  {
    CGAL_assertion(dh != NULL);
    if (get_local_map().is_marked(dh, m_is_hole))
    { dh = get_local_map().opposite2(dh); }
    CGAL_assertion(!get_local_map().is_marked(dh, m_is_hole));
    return dh;
  }

  Dart_handle nonhole_dart_of_same_edge(Dart_handle dh, bool& flip)
  {
    CGAL_assertion(dh != NULL);
    if (get_local_map().is_marked(dh, m_is_hole))
    {
      dh = get_local_map().opposite2(dh);
      flip = !flip;
    }
    CGAL_assertion(!get_local_map().is_marked(dh, m_is_hole));
    return dh;
  }

  void add_to_cycle(Dart_handle dh, Path& cycle, bool flip=false)
  {
    dh = nonhole_dart_of_same_edge(dh, flip);
    Original_dart_const_handle dh_original = m_copy_to_origin[dh];
    if (cycle.can_be_pushed(dh_original, flip))
    { cycle.push_back(dh_original, flip); }
    else 
    {
      CGAL_assertion(cycle.can_be_pushed(dh_original, !flip));
      cycle.push_back(dh_original, !flip);
    }
  }

  template <class WeightFunctor>
  bool compute_cycle(Dart_handle root, Path& cycle,
                     typename WeightFunctor::Weight_t* length = NULL,
                     const typename WeightFunctor::Weight_t* max_length = NULL,
                     const WeightFunctor& wf = Default_weight_functor())
  {
    std::vector<typename WeightFunctor::Weight_t> distance_from_root;
    m_spanning_tree.clear();
    m_trace_index.clear();
    for (auto it = get_local_map().template one_dart_per_cell<0>().begin(),
                     itend = get_local_map().template one_dart_per_cell<0>().end(); it != itend; ++it)
    { get_local_map().template info<0>(it) = -1; }
    compute_spanning_tree(root, m_spanning_tree, distance_from_root, m_trace_index, wf);
    compute_noncon_edges(m_spanning_tree, m_noncon_edges);

    bool first_check = true;
    typename WeightFunctor::Weight_t min_distance = 0;
    Dart_handle min_noncon_edge;
    int min_a = -1, min_b = -1;
    for (auto dh : m_noncon_edges)
    {
      Dart_handle a = dh, b = get_local_map().next(dh);
      int index_a = get_local_map().template info<0>(a), index_b = get_local_map().template info<0>(b);
      typename WeightFunctor::Weight_t sum_distance=
        distance_from_root[index_a] + distance_from_root[index_b] 
        + wf(m_copy_to_origin[nonhole_dart_of_same_edge(dh)]);
      if (first_check || min_distance > sum_distance)
      {
        min_distance = sum_distance;
        min_noncon_edge = dh;
        min_a = index_a;
        min_b = index_b;
        first_check = false;
      }
    }
    if (first_check) return false; // no cycle found
    if (length != NULL) *length = min_distance < 0 ? 0 : min_distance;
    if (max_length != NULL && min_distance >= *max_length) return false; // abort
    cycle.clear();
    // Trace back the path from `a` to root
    for (int ind = min_a - 1; ind != -1; ind = m_trace_index[ind]) 
      add_to_cycle(m_spanning_tree[ind], cycle, true);
    // Reverse: now it is the path from root to `a`
    cycle.reverse();
    add_to_cycle(min_noncon_edge, cycle);
    // Trace back the path from `b` to root
    for (int ind = min_b - 1; ind != -1; ind = m_trace_index[ind])
      add_to_cycle(m_spanning_tree[ind], cycle, true);
    CGAL_assertion(cycle.is_closed());

    return true;
  }

protected:
  Local_map      m_local_map; /// the local map
  Origin_to_copy m_origin_to_copy;
  Copy_to_origin m_copy_to_origin;
  std::size_t    m_nb_of_vertices = 0;
  Dart_container m_spanning_tree, m_noncon_edges, m_face_list, m_edge_list;
  std::vector<int> m_trace_index;
  size_type m_is_hole;
  Path m_cycle;
};

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif
