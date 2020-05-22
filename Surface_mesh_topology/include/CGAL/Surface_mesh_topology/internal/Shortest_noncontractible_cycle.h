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
#include <CGAL/Surface_mesh_topology/internal/Edge_weight_functor.h>
#include <queue>
#include <tuple>
#include <unordered_map>
#include <boost/heap/fibonacci_heap.hpp>

namespace CGAL {
namespace Surface_mesh_topology {
namespace internal {

template<class SNC, bool Copy>
struct Get_original_dart
{
  static typename SNC::Original_dart_const_handle
  run(SNC* snc, typename SNC::Dart_handle dh)
  { return snc->m_copy_to_origin[dh]; }
};

template<class SNC>
struct Get_original_dart<SNC, false>
{
  static typename SNC::Original_dart_const_handle
  run(SNC* /* snc */, typename SNC::Dart_handle dh)
  { return dh; }
};

template <class Mesh_, bool Copy>
class Shortest_noncontractible_cycle
{
public:
  using Self=Shortest_noncontractible_cycle<Mesh_, Copy>;
  using Mesh=Mesh_;

  friend struct Get_original_dart<Self, Copy>;

  using Original_map_wrapper=internal::Generic_map_selector<Mesh, Items_for_shortest_noncontractible_cycle>;
  using Original_dart_const_handle=typename Original_map_wrapper::Dart_const_handle_original;

  using Local_map        =typename Original_map_wrapper::Generic_map;
  using Dart_handle      =typename Local_map::Dart_handle;
  using Dart_const_handle=typename Local_map::Dart_const_handle;
  using size_type        = typename Local_map::size_type;

  using Dart_container=std::vector<Dart_handle>;
  using Path          =CGAL::Surface_mesh_topology::Path_on_surface<Mesh>;

  // Associations between original darts and their copy.
  using Origin_to_copy=boost::unordered_map<Original_dart_const_handle, Dart_handle>;
  using Copy_to_origin=boost::unordered_map<Dart_handle, Original_dart_const_handle>;

  /// @return the local map
  const Local_map& get_local_map() const
  { return *m_local_map; }

  /// @return the local map
  Local_map& get_local_map()
  { return *m_local_map; }

  Shortest_noncontractible_cycle(const Mesh& amesh, bool display_time=false) :
    m_cycle(amesh)
  {
    CGAL::Timer t;
    if (display_time)
    { t.start(); }

    m_local_map=new Local_map;

    // Initialize m_is_perforated
    m_is_perforated=get_local_map().get_new_mark();

    Original_map_wrapper::copy(*m_local_map, amesh,
                               m_origin_to_copy, m_copy_to_origin, m_is_perforated);

    get_local_map().negate_mark(m_is_perforated);
    // Remove all boundary by adding faces, marked with m_is_perforated
    get_local_map().template close<2>();
    get_local_map().negate_mark(m_is_perforated);

    create_vertex_info();

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] Shortest_noncontractible_cycle constructor: "
               <<t.time()<<" seconds."<<std::endl;
    }
  }

  // Here amesh is already closed, and have 0-attributes with int.
  // Thus the mesh is not copied.
  Shortest_noncontractible_cycle(Mesh* amesh, size_type perforated_mark,
                                 bool /*display_time*/=false) :
    m_local_map(amesh),
    m_is_perforated(perforated_mark),
    m_cycle(*amesh)
  {}

  ~Shortest_noncontractible_cycle()
  {
    if (Copy)
    {
      delete m_local_map;
    }
  }

  template <class WeightFunctor>
  Path compute_cycle(Original_dart_const_handle root_vertex,
                     typename WeightFunctor::Weight_t* length,
                     const WeightFunctor& wf,
                     bool display_time=false)
  {
    CGAL::Timer t;
    if (display_time)
    { t.start(); }

    m_cycle.clear();
    Dart_handle root=m_origin_to_copy[root_vertex];
    this->compute_cycle(root, m_cycle, length, nullptr, wf);

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] compute_cycle: "<<t.time()<<" seconds."<<std::endl;
    }

    return m_cycle;
  }

  template <class WeightFunctor=Unit_weight_functor>
  Path compute_cycle(Original_dart_const_handle root_vertex,
                     typename WeightFunctor::Weight_t* length,
                     bool display_time=false)
  { return compute_cycle(root_vertex, length, WeightFunctor(), display_time); }

  template <class WeightFunctor=Unit_weight_functor>
  Path compute_cycle(Original_dart_const_handle root_vertex,
                     bool display_time=false)
  { return compute_cycle(root_vertex, nullptr, display_time); }

  template <class WeightFunctor>
  Path compute_shortest_non_contractible_cycle(typename WeightFunctor::Weight_t* length,
                                               const WeightFunctor& wf,
                                               bool display_time=false)
  {
    CGAL::Timer t;
    if (display_time)
    { t.start(); }

    m_cycle.clear();
    bool first_check = true;
    typename WeightFunctor::Weight_t min_length=0;
    for (auto it=get_local_map().template attributes<0>().begin(),
         itend=get_local_map().template attributes<0>().end(); it!=itend; ++it)
    {
      Dart_handle dh=get_local_map().template dart_of_attribute<0>(it);
      typename WeightFunctor::Weight_t temp_length;
      if (first_check)
      {
        if (compute_cycle(dh, m_cycle, &temp_length, nullptr, wf))
        {
          min_length = temp_length;
          first_check = false;
        }
      }
      else
      {
        if (compute_cycle(dh, m_cycle, &temp_length, &min_length, wf))
        { min_length = temp_length; }
      }
    }
    if (length!=nullptr) { *length=min_length; }

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] compute_shortest_non_contractible_cycle: "<<t.time()<<" seconds."<<std::endl;
    }

    return m_cycle;
  }

  template <class WeightFunctor=Unit_weight_functor>
  Path compute_shortest_non_contractible_cycle(typename WeightFunctor::Weight_t* length,
                                              bool display_time=false)
  { return compute_shortest_non_contractible_cycle(length, WeightFunctor(), display_time); }

  Path compute_edge_width(bool display_time=false)
  { return compute_shortest_non_contractible_cycle(nullptr, display_time); }

protected:
  int vertex_info(Dart_handle dh) const
  { return get_local_map().template info<0>(dh); }
  int& vertex_info(Dart_handle dh)
  { return get_local_map().template info<0>(dh); }

  void create_vertex_info()
  {
    for (auto it=get_local_map().darts().begin(),
         itend=get_local_map().darts().end(); it!=itend; ++it)
    {
      if (get_local_map().template attribute<0>(it)==nullptr)
      { get_local_map().template set_attribute<0>
          (it, get_local_map().template create_attribute<0>(-1)); }
    }
  }

  void initialize_vertex_info()
  {
    for (auto it=get_local_map().template attributes<0>().begin(),
         itend = get_local_map().template attributes<0>().end(); it != itend; ++it)
    { get_local_map().template info_of_attribute<0>(it)=-1; }
  }

  template <class WeightFunctor, class Distance_>
  void compute_spanning_tree(Dart_handle root, Dart_container& spanning_tree,
                             std::vector<Distance_>& distance_from_root,
                             std::vector<int>& trace_index,
                             const WeightFunctor& wf)
  {
    if (std::is_same<WeightFunctor, Unit_weight_functor>::value)
    { compute_BFS_tree(root, spanning_tree, distance_from_root, trace_index); }
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
                             const WeightFunctor& wf)
  {
    // Preparation
    Dijkstra_comparator<Distance_> dc(distance_from_root);
    //std::priority_queue<int, std::vector<int>, Dijkstra_comparator<Distance_> > pq(dc);
    typedef boost::heap::fibonacci_heap<int, boost::heap::compare<Dijkstra_comparator<Distance_>>> Heap;
    typedef typename Heap::handle_type heap_handle;

    Heap pq(dc);
    std::unordered_map<int, heap_handle> inqueue;

    int vertex_index=0;
    size_type vertex_visited=get_local_map().get_new_mark();
    // Begin Dijkstra
    distance_from_root.reserve(get_local_map().template attributes<0>().size());
    distance_from_root.push_back(0);
    inqueue[0]=pq.push(0);
    vertex_info(root)=vertex_index;
    get_local_map().template mark_cell<0>(root, vertex_visited);

    while (!pq.empty())
    {
      int u_index=pq.top();
      pq.pop();

      inqueue.erase(u_index);

      Dart_handle u=(u_index==0)?root:get_local_map().next(spanning_tree[u_index-1]);
      CGAL_assertion(u_index==vertex_info(u));
      Dart_handle it=u;
      do
      {
        Dart_handle v=get_local_map().next(it);
        Distance_ w=wf(Get_original_dart<Self, Copy>::run(this, nonhole_dart_of_same_edge(it)));
        if (!get_local_map().is_marked(v, vertex_visited))
        {
          int v_index=++vertex_index;
          CGAL_assertion(v_index==static_cast<int>(distance_from_root.size()));
          CGAL_assertion(inqueue.count(v_index)==0);
          distance_from_root.push_back(distance_from_root[u_index]+w);
          spanning_tree.push_back(it);
          trace_index.push_back(u_index-1);
          vertex_info(v)=v_index;
          get_local_map().template mark_cell<0>(v, vertex_visited);
          inqueue[v_index]=pq.push(v_index);
        }
        else
        {
          int v_index=vertex_info(v);
          if (distance_from_root[v_index]>distance_from_root[u_index]+w)
          {
            CGAL_assertion(v_index>0);
            distance_from_root[v_index]=distance_from_root[u_index]+w;
            spanning_tree[v_index-1]=it;
            trace_index[v_index-1]=u_index-1;

            auto it_inqueue=inqueue.find(v_index);
            if (it_inqueue==inqueue.end())
            { inqueue[v_index]=pq.push(v_index); }
            else
            { pq.decrease(it_inqueue->second, v_index); }
          }
        }
        it=get_local_map().next(get_local_map().opposite2(it));
      }
      while(it!=u);
    }
    get_local_map().free_mark(vertex_visited);
    // End Dijkstra
  }

  /// Create a spanning tree using BFS
  template <class Distance_>
  void compute_BFS_tree(Dart_handle root, Dart_container& spanning_tree,
                        std::vector<Distance_>& distance_from_root,
                        std::vector<int>& trace_index)
  {
    // Preparation
    std::queue<int> q;
    int vertex_index=0;
    size_type vertex_visited=get_local_map().get_new_mark();
    // Begin BFS
    q.push(0);
    vertex_info(root)=vertex_index;
    get_local_map().template mark_cell<0>(root, vertex_visited);
    distance_from_root.reserve(get_local_map().template attributes<0>().size());
    distance_from_root.push_back(0);
    while (!q.empty())
    {
      int u_index=q.front();
      q.pop();
      Dart_handle u=(u_index==0)?root:get_local_map().next(spanning_tree[u_index-1]);
      CGAL_assertion(u_index==vertex_info(u));
      Dart_handle it=u;
      do
      {
        Dart_handle v=get_local_map().next(it);
        if (!get_local_map().is_marked(v, vertex_visited))
        {
          int v_index=++vertex_index;
          CGAL_assertion(v_index==static_cast<int>(distance_from_root.size()));
          distance_from_root.push_back(1+distance_from_root[u_index]);
          spanning_tree.push_back(it);
          // `it` will lead to v
          q.push(v_index);
          trace_index.push_back(u_index-1);
          vertex_info(v)=v_index;
          get_local_map().template mark_cell<0>(v, vertex_visited);
        }
        it=get_local_map().next(get_local_map().opposite2(it));
      }
      while(it!=u);
    }
    get_local_map().free_mark(vertex_visited);
    // End BFS
  }

  /// Check if there is only one unmarked edge around a face.
  /// If there is, let dh_only_edge=the edge separating it and its only adjacent face.
  bool is_degree_one_face(Dart_handle dh_face, Dart_handle& dh_only_edge,
                          size_type edge_deleted)
  {
    Dart_handle dh=dh_face;
    dh_only_edge=nullptr;
    do
    {
      if (!get_local_map().is_marked(dh, edge_deleted))
      {
        if (dh_only_edge!=nullptr)
        { dh_only_edge=nullptr; return false; }
        dh_only_edge=dh;
      }
      dh=get_local_map().next(dh);
    }
    while(dh!=dh_face);
    return (dh_only_edge!=nullptr);
  }

  /// Find E_nc
  void compute_noncon_edges(const Dart_container& spanning_tree,
                            Dart_container& noncon_edges)
  {
    noncon_edges.clear();
    size_type face_deleted=get_local_map().get_new_mark();
    size_type edge_deleted=get_local_map().get_new_mark();
    size_type tested=get_local_map().get_new_mark();
    for (auto dh : spanning_tree)
    {
      if (!get_local_map().is_marked(dh, edge_deleted))
      { get_local_map().template mark_cell<1>(dh, edge_deleted); }
    }
    std::queue<Dart_handle> degree_one_faces;
    // Add to queue the degree-1 faces
    for (auto it=get_local_map().darts().begin(), itend=get_local_map().darts().end();
         it!=itend; ++it)
    {
      if (get_local_map().is_marked(it, m_is_perforated))
      {
        if (!get_local_map().is_marked(it, face_deleted))
        { get_local_map().template mark_cell<2>(it, face_deleted); }
      }
      else if (!get_local_map().is_marked(it, tested))
      {
        Dart_handle dh=it, dh_only_edge=nullptr;
        bool degree_one=true;
        do // Here we do not use is_degree_one_face method because we want to
        {  // mark tested all darts of the face in the same loop.
          get_local_map().template mark_cell<1,1>(dh, tested); // For CMap and GMap
          if (degree_one)
          {
            if (!get_local_map().is_marked(dh, edge_deleted))
            {
              if (dh_only_edge!=nullptr) { degree_one=false; }
              else { dh_only_edge=dh; }
            }
          }
          dh=get_local_map().next(dh);
        }
        while(dh!=it);

        if (degree_one && dh_only_edge!=nullptr)
        {
          degree_one_faces.push(dh_only_edge);
          CGAL_assertion(!get_local_map().is_marked(dh_only_edge, edge_deleted));
          CGAL_assertion(!get_local_map().is_marked(dh_only_edge, face_deleted));
        }
      }
    }

    // Remove the degree-1 faces
    while (!degree_one_faces.empty())
    {
      Dart_handle dh_face=degree_one_faces.front();
      degree_one_faces.pop();
      if (!get_local_map().is_marked(dh_face, edge_deleted))
      { get_local_map().template mark_cell<1>(dh_face, edge_deleted); }
      if (!get_local_map().is_marked(dh_face, face_deleted))
      { get_local_map().template mark_cell<2>(dh_face, face_deleted); }
      Dart_handle dh_adj_face=get_local_map().opposite2(dh_face);
      if (!get_local_map().is_marked(dh_adj_face, face_deleted))
      {
        Dart_handle dh_only_edge=nullptr;
        if (is_degree_one_face(dh_adj_face, dh_only_edge, edge_deleted))
        {
          degree_one_faces.push(dh_only_edge);
          CGAL_assertion(!get_local_map().is_marked(dh_only_edge, edge_deleted));
          CGAL_assertion(!get_local_map().is_marked(dh_only_edge, face_deleted));
        }
      }
    }
    for (auto it=get_local_map().darts().begin(),
         itend=get_local_map().darts().end(); it!=itend; ++it)
    {
      if (it<get_local_map().opposite2(it) && vertex_info(it)>=0 &&
          !get_local_map().is_marked(it, edge_deleted))
      { noncon_edges.push_back(it); }

      get_local_map().unmark(it, tested);
      get_local_map().unmark(it, edge_deleted);
      get_local_map().unmark(it, face_deleted);
    }

    get_local_map().free_mark(tested);
    get_local_map().free_mark(edge_deleted);
    get_local_map().free_mark(face_deleted);
  }

  Dart_handle nonhole_dart_of_same_edge(Dart_handle dh)
  {
    CGAL_assertion(dh!=nullptr);
    if (get_local_map().is_marked(dh, m_is_perforated))
    { dh = get_local_map().opposite2(dh); }
    CGAL_assertion(!get_local_map().is_marked(dh, m_is_perforated));
    return dh;
  }

  Dart_handle nonhole_dart_of_same_edge(Dart_handle dh, bool& flip)
  {
    CGAL_assertion(dh!=nullptr);
    if (get_local_map().is_marked(dh, m_is_perforated))
    {
      dh=get_local_map().opposite2(dh);
      flip=!flip;
    }
    CGAL_assertion(!get_local_map().is_marked(dh, m_is_perforated));
    return dh;
  }

  void add_to_cycle(Dart_handle dh, Path& cycle, bool flip)
  {
    dh=nonhole_dart_of_same_edge(dh, flip);
    Original_dart_const_handle
        dh_original=Get_original_dart<Self, Copy>::run(this, dh);
    if (cycle.can_be_pushed(dh_original, flip))
    { cycle.push_back(dh_original, flip, false); }
    else
    {
      CGAL_assertion(cycle.can_be_pushed(dh_original, !flip));
      cycle.push_back(dh_original, !flip, false);
    }
  }

  template <class WeightFunctor>
  bool compute_cycle(Dart_handle root, Path& cycle,
                     typename WeightFunctor::Weight_t* length,
                     const typename WeightFunctor::Weight_t* max_length,
                     const WeightFunctor& wf)
  {
    std::vector<typename WeightFunctor::Weight_t> distance_from_root;
    m_spanning_tree.clear();
    m_trace_index.clear();
    initialize_vertex_info();
    compute_spanning_tree(root, m_spanning_tree, distance_from_root, m_trace_index, wf);
    compute_noncon_edges(m_spanning_tree, m_noncon_edges);

    bool first_check=true;
    typename WeightFunctor::Weight_t min_distance=0;
    Dart_handle min_noncon_edge;
    int min_a=-1, min_b=-1;
    for (auto dh : m_noncon_edges)
    {
      Dart_handle a=dh, b=get_local_map().next(dh);
      int index_a=vertex_info(a), index_b=vertex_info(b);
      typename WeightFunctor::Weight_t sum_distance=
        distance_from_root[index_a]+distance_from_root[index_b]
        +wf(Get_original_dart<Self, Copy>::run(this, nonhole_dart_of_same_edge(dh)));
      if (first_check || sum_distance<min_distance)
      {
        min_distance=sum_distance;
        min_noncon_edge=dh;
        min_a=index_a;
        min_b=index_b;
        first_check=false;
      }
    }

    if (first_check) { return false; } // no cycle found

    if (max_length!=nullptr && min_distance>=*max_length)
    { return false; } // Here the cycle is too long

    cycle.clear();
    // Trace back the path from `a` to root
    for (int ind=min_a-1; ind!=-1; ind=m_trace_index[ind])
    {  add_to_cycle(m_spanning_tree[ind], cycle, true); }
    // Reverse: now it is the path from root to `a`
    cycle.reverse();
    add_to_cycle(min_noncon_edge, cycle, false);
    // Trace back the path from `b` to root
    for (int ind=min_b-1; ind!=-1; ind=m_trace_index[ind])
    { add_to_cycle(m_spanning_tree[ind], cycle, true); }

    cycle.update_is_closed();
    CGAL_assertion(cycle.is_closed());
    if (length!=nullptr) { *length=(min_distance<0?0:min_distance); }

    return true;
  }

protected:
  Local_map*       m_local_map; /// the local map
  size_type        m_is_perforated;   /// mark for perforated darts
  Origin_to_copy   m_origin_to_copy; /// associative array from original darts to copies
  Copy_to_origin   m_copy_to_origin; /// associative array from copies to original darts
  Dart_container   m_spanning_tree, m_noncon_edges; // Darts in the spanning tree and in E_nc
  std::vector<int> m_trace_index;
  Path             m_cycle;
};

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif
