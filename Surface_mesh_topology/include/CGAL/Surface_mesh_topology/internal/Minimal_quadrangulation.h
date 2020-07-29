// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_MINIMAL_QUADRANGULATION_H
#define CGAL_MINIMAL_QUADRANGULATION_H 1

// Should be defined before to include Path_on_surface_with_rle.h
// If nothing is defined, use V1
// #define CGAL_PWRLE_TURN_V1  // Compute turns by turning (method of CMap)
// #define CGAL_PWRLE_TURN_V2  // Compute turns by using an id of darts, given by an hash-table (built and given by Minimal_quadrangulation)
#define CGAL_PWRLE_TURN_V3  // Compute turns by using an id of darts, associated in Info of Darts (build by Minimal_quadrangulation)

#include <CGAL/license/Surface_mesh_topology.h>

#include <CGAL/Face_graph_wrapper.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Surface_mesh_topology/internal/Path_on_surface_with_rle.h>
#include <CGAL/Surface_mesh_topology/internal/Path_generators.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <CGAL/Cell_attribute.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>
#include <unordered_map>
#include <queue>
#include <tuple>
#include <iostream>

namespace CGAL {
namespace Surface_mesh_topology {
namespace internal {

struct Minimal_quadrangulation_local_map_items
{
  template <class CMap>
  struct Dart_wrapper
  {
#ifdef CGAL_PWRLE_TURN_V3
    typedef std::size_t Dart_info;
#endif // CGAL_PWRLE_TURN_V3
    typedef CGAL::Cell_attribute<CMap, int32_t> Vertex_attribute;
    typedef std::tuple<Vertex_attribute, void, void> Attributes;
  };
};

template<typename Mesh_>
class Minimal_quadrangulation
{
public:
  typedef Minimal_quadrangulation<Mesh_> Self;
  typedef Mesh_                          Mesh;

  typedef typename Get_map<Mesh, Mesh>::type       Original_map; // Mesh seen as a 2-map
  typedef typename Original_map::Dart_const_handle Original_dart_const_handle;
  typedef typename Original_map::size_type         Original_size_type;

  typedef CGAL::Combinatorial_map<2, Minimal_quadrangulation_local_map_items> Local_map;
  typedef typename Local_map::Dart_handle             Dart_handle;
  typedef typename Local_map::Dart_const_handle       Dart_const_handle;
  typedef typename Local_map::size_type               size_type;

  // Associate each dart of the original map, not removed, a pair of darts in
  // the reduced map.
  typedef std::unordered_map<Original_dart_const_handle,
  std::pair<Dart_handle, Dart_handle> > TPaths;

  // Associations between original darts and their copy.
  // Use only locally during the construction of the minimal quadrangulation.
  typedef std::unordered_map<Original_dart_const_handle, Dart_handle> Origin_to_copy;
  typedef std::unordered_map<Dart_handle, Original_dart_const_handle> Copy_to_origin;

#ifdef CGAL_PWRLE_TURN_V2
  typedef std::unordered_map<Original_dart_const_handle, std::size_t> TDartIds;
#endif //CGAL_PWRLE_TURN_V2

  /// @return the original map
  const Original_map& get_original_map() const
  { return m_original_map; }

  /// @return the reduced map
  const Local_map& get_local_map() const
  { return m_local_map; }

  /// @return the reduced map
  Local_map& get_local_map()
  { return m_local_map; }

  /// Constructor taking a mesh as parameter.
  Minimal_quadrangulation(const Mesh& amesh, bool display_time=false) :
    m_original_map(amesh)
  {
    if (!get_original_map().is_without_boundary(1))
    {
      std::cerr<<"ERROR: the given amap has 1-boundaries; "
               <<"such a surface is not possible to process here."
               <<std::endl;
      return;
    }

    CGAL::Timer t, t2;
    if (display_time)
    { t.start(); t2.start(); }

    // The mapping between darts of the original map into the copied map.
    Origin_to_copy origin_to_copy;

    // The mapping between darts of the copy into darts of the original map.
    Copy_to_origin copy_to_origin;

    // We reserve three marks: m_mart_T and m_mark_L to mark darts in
    // get_original_map() that belong to T or to L (i.e. edges contracted or
    // removed), m_mark_perforated to mark perforated faces.
    m_mark_T=get_original_map().get_new_mark();
    m_mark_L=get_original_map().get_new_mark();
    m_mark_perforated=get_local_map().get_new_mark();

    // 1) We create m_local_map as a copy of amap, by contracting all
    //    non loop.
    surface_simplification_in_one_vertex(origin_to_copy, copy_to_origin);

    if (display_time)
    {
      t2.stop();
      std::cout<<"[TIME] Simplification in one vertex: "
               <<t2.time()<<" seconds"<<std::endl;

      t2.reset(); t2.start();
    }

#ifdef CGAL_TRACE_CMAP_TOOLS
    std::cout<<"Creation of reduced map with non loops contracted: ";
    get_local_map().display_characteristics(std::cout);
    std::cout<<", valid="<<get_local_map().is_valid()<<std::endl;
#endif

    // 2) We simplify m_local_map in a surface with only one face
    //    (or maybe more if we have boundaries!)
    surface_simplification_in_one_face(origin_to_copy, copy_to_origin);

    // 3) We remove all non perforated loops
    remove_non_perforated_loops(copy_to_origin);

    if (display_time)
    {
      t2.stop();
      std::cout<<"[TIME] Simplification in one face: "
               <<t2.time()<<" seconds"<<std::endl;
      t2.reset(); t2.start();
    }

#ifdef CGAL_TRACE_CMAP_TOOLS
    std::cout<<"All faces merges: ";
    get_local_map().display_characteristics(std::cout);
    std::cout<<", valid="<<get_local_map().is_valid()<<std::endl;
#endif

    if (!get_local_map().is_empty()) // is_empty if the surface is a sphere
    {
      // 4) We quadrangulate the face, except for the torus surfaces
      if (!local_map_is_a_torus())
      {
        surface_quadrangulate();

        if (display_time)
        {
          t2.stop();
          std::cout<<"[TIME] Face quadrangulation: "
                  <<t2.time()<<" seconds"<<std::endl;
          t2.reset(); t2.start();
        }
      }

#if defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
      // Now we label all the darts of the reduced map, to allow the computation
      // of turns in constant time: only for methods V2 and V3
      initialize_vertices_attributes();
      initialize_ids();
#endif // defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
    }

#ifdef CGAL_TRACE_CMAP_TOOLS
    std::cout<<"After quadrangulation: ";
    get_local_map().display_characteristics(std::cout);
    std::cout<<", valid="<<get_local_map().is_valid()<<std::endl;

    std::cout<<"Paths are all valid ? "<<(are_paths_valid()?"YES":"NO")
             <<std::endl;
    auto marktemp=get_local_map().get_new_mark();
    Dart_handle dh2=nullptr;
    for (auto it=get_local_map().darts().begin();
         it!=get_local_map().darts().end(); ++it)
    {
      if (!get_local_map().is_marked(it, marktemp))
      {
        std::cout<<"Degree="<<CGAL::template
                   degree<Local_map, 0>(get_local_map(), it)
            <<std::endl;
        std::cout<<"Co-degree="<<CGAL::template
                   codegree<Local_map, 2>(get_local_map(), it)
            <<std::endl;
        dh2=it;
        do
        {
          get_local_map().mark(dh2, marktemp);
          std::cout<<get_local_map().darts().index(dh2)<<"   "
                   <<get_local_map().darts().
                     index(get_local_map().template beta<0>(dh2))
                   <<std::endl;
          dh2=get_local_map().template beta<0,2>(dh2);
        }
        while(dh2!=it);
      }
    }
    get_local_map().free_mark(marktemp);
    get_local_map().display_darts(std::cout);
#endif

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] Total time for computation of reduced map: "
               <<t.time()<<" seconds"<<std::endl;
    }

    CGAL_assertion(perforated_faces_correctly_marked());
    CGAL_assertion(local_map_is_a_torus() || are_paths_valid()); // Because torus is a special case
  }

  ~Minimal_quadrangulation()
  {
    get_original_map().free_mark(m_mark_T);
    get_original_map().free_mark(m_mark_L);
    get_local_map().free_mark(m_mark_perforated);
  }

  /// @return true iff 'p' is contractible.
  bool is_contractible(const Path_on_surface<Mesh>& p,
                       bool display_time=false) const
  {
    if (p.is_empty())
    { return true; }

    if (!p.is_closed())
    {
      std::cerr<<"Error: is_contractible requires a closed path."<<std::endl;
      return false;
    }

    if (get_local_map().is_empty())
    { return true; } // A closed path on a sphere is always contractible.

    CGAL::Timer t;
    if (display_time)
    { t.start(); }

    bool res=false;
    if (local_map_is_a_torus())
    { // Case of torus
      Path_on_surface<Local_map>
        pt=transform_original_path_into_quad_surface_for_torus(p);

      int a, b;
      count_edges_of_path_on_torus(pt, a, b);
      res=(a==0 && b==0);
    }
    else
    {
      internal::Path_on_surface_with_rle<Self>
        pt=transform_original_path_into_quad_surface_with_rle(p);

      pt.canonize();
      res=pt.is_empty();
    }

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] is_contractible: "<<t.time()<<" seconds"
               <<std::endl;
    }

    return res;
  }

  /// @return true iff 'p1' and 'p2' are freely homotopic.
  bool are_freely_homotopic(const Path_on_surface<Mesh>& p1,
                            const Path_on_surface<Mesh>& p2,
                            bool display_time=false) const
  {
    if (p1.is_empty() && p2.is_empty()) { return true; }

    if ((!p1.is_empty() && !p1.is_closed()) ||
        (!p2.is_empty() && !p2.is_closed()))
    {
      std::cerr<<"Error: are_freely_homotopic requires two closed paths."
               <<std::endl;
      return false;
    }

    // Here we have two non empty closed paths.
    if (get_local_map().is_empty())
    { return true; } // Two closed paths on a sphere are always homotopic.

    CGAL::Timer t;
    if (display_time)
    { t.start(); }

    bool res=false;
    if (local_map_is_a_torus())
    { // Case of torus
      Path_on_surface<Local_map>
        pt1=transform_original_path_into_quad_surface_for_torus(p1);
      Path_on_surface<Local_map>
        pt2=transform_original_path_into_quad_surface_for_torus(p2);

      int a1, a2, b1, b2;
      count_edges_of_path_on_torus(pt1, a1, b1);
      count_edges_of_path_on_torus(pt2, a2, b2);
      res=(a1==a2 && b1==b2);
    }
    else
    {
      internal::Path_on_surface_with_rle<Self>
        pt1=transform_original_path_into_quad_surface_with_rle(p1);
      internal::Path_on_surface_with_rle<Self>
        pt2=transform_original_path_into_quad_surface_with_rle(p2);
      pt1.canonize();
      pt2.canonize();
      res=(pt1==pt2); // Do here to be counted in the computation time

#ifdef CGAL_TRACE_PATH_TESTS
      std::cout<<"Length of reduced paths: "<<pt1.length()<<" and "
               <<pt2.length()<<std::endl;
#endif
    }

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] are_freely_homotopic: "<<t.time()<<" seconds"<<std::endl;
    }

    return res;
  }

  /// @return true iff 'p1' and 'p2' are base point freely homotopic.
  bool are_base_point_homotopic(const Path_on_surface<Mesh>& p1,
                                const Path_on_surface<Mesh>& p2,
                                bool display_time=false) const
  {
    if (p1.is_empty() && p2.is_empty()) { return true; }
    if (p1.is_empty() || p2.is_empty()) { return false; }

    if (!get_original_map().template belong_to_same_cell<0>
        (p1.front_flip()?get_original_map().template beta<1>(p1.front()):p1.front(),
         p2.front_flip()?get_original_map().template beta<1>(p2.front()):p2.front()) ||
        !get_original_map().template belong_to_same_cell<0>
        (p1.back_flip()?p1.back():get_original_map().template beta<1>(p1.back()),
         p2.back_flip()?p2.back():get_original_map().template beta<1>(p2.back())))
    {
      std::cerr<<"Error: are_base_point_homotopic requires two paths that"
               <<" share the same vertices as extremities."<<std::endl;
      return false;
    }

    if (get_local_map().is_empty())
    { return true; } // Two paths on a sphere are always base_point_homotopic.

    CGAL::Timer t;
    if (display_time)
    { t.start(); }

    Path_on_surface<Mesh> path=p1;
    Path_on_surface<Mesh> path2=p2; path2.reverse();
    path+=path2;

    bool res=is_contractible(path);

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] are_base_point_homotopic: "<<t.time()<<" seconds."
               <<std::endl;
    }

    return res;
  }

  /// @return the positive turn between the two given darts.
  /// i.e. turning right if the face bounded by a dart lays on its right
  /// and left otherwise
  std::size_t positive_turn(Dart_const_handle dh1,
                            Dart_const_handle dh2) const
  {
#if defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
    return compute_positive_turn_given_ids
        (get_local_map().template beta<2>(dh1), dh2);
#else // defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)

    CGAL_assertion((!get_local_map().template is_free<1>(dh1)));
    CGAL_assertion((!get_local_map().template is_free<2>(dh1)));

    if (dh2==get_local_map().template beta<2>(dh1) &&
        dh2==get_local_map().template beta<1>(dh1))
    { return 0; }

    if (get_local_map().is_marked(dh1, m_mark_hole))
    { return (std::numeric_limits<std::size_t>::max)(); }

    Dart_const_handle ddh1=dh1;
    std::size_t res=1;
    while (get_local_map().template beta<1>(ddh1)!=dh2)
    {
      CGAL_assertion(!get_local_map().template is_free<2>
                     (get_local_map().template beta<1>(ddh1)));

      ++res;
      ddh1=get_local_map().template beta<1, 2>(ddh1);
      if (get_local_map().is_marked(ddh1, m_mark_hole))
      { return (std::numeric_limits<std::size_t>::max)(); }

      CGAL_assertion(!get_local_map().template is_free<1>(ddh1));
      CGAL_assertion(get_local_map().template beta<1>(ddh1)==dh2 || ddh1!=dh1);
    }
    return res;
#endif // defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
  }

  /// @return the negative turn between the two given darts.
  /// i.e. turning right if the face bounded by a dart lays on its left
  /// and left otherwise
  std::size_t negative_turn(Dart_const_handle dh1,
                            Dart_const_handle dh2) const
  {
#if defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
    return compute_negative_turn_given_ids
        (get_local_map().template beta<2>(dh1), dh2);
#else // defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)

    CGAL_assertion((!get_local_map().template is_free<1>(dh1)));
    CGAL_assertion((!get_local_map().template is_free<2>(dh1)));

    if (dh2==get_local_map().template beta<2>(dh1) &&
        dh1==get_local_map().template beta<0>(dh2))
    { return 0; }

    if (get_local_map().is_marked
        (get_local_map().template beta <2>(dh1), m_mark_hole))
    { return (std::numeric_limits<std::size_t>::max)(); }

    dh1=get_local_map().template beta<2>(dh1);
    dh2=get_local_map().template beta<2>(dh2);
    Dart_const_handle ddh1=dh1;
    std::size_t res=1;
    while (get_local_map().template beta<0>(ddh1)!=dh2)
    {
      CGAL_assertion(!get_local_map().template is_free<2>
                     (get_local_map().template beta<0>(ddh1)));

      ++res;
      ddh1=get_local_map().template beta<0, 2>(ddh1);
      if (get_local_map().is_marked(ddh1, m_mark_hole))
      { return (std::numeric_limits<std::size_t>::max)(); }

      CGAL_assertion(!get_local_map().template is_free<0>(ddh1));
      CGAL_assertion(get_local_map().template beta<0>(ddh1)==dh2 || ddh1!=dh1);
    }
    return res;
#endif // defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
  }

protected:

  /// @return true iff the dart dh is contracted, i.e. if it belongs to T
  bool is_contracted(Original_dart_const_handle dh) const
  { return get_original_map().is_marked(dh, m_mark_T); }

  void count_edges_of_path_on_torus
  (const Path_on_surface<Local_map>& path,
   int& a, int& b) const
  {
    CGAL_assertion(local_map_is_a_torus());

    Dart_const_handle dha=get_local_map().darts().begin();
    Dart_const_handle dhb=get_local_map().template beta<1>(dha);

    a=0; b=0;
    for (std::size_t i=0; i<path.length(); ++i)
    {
      if(path.get_ith_flip(i))
      {
        if (path[i]==dha) { --a; }
        else if (path[i]==get_local_map().template beta<2>(dha)) { ++a; }
        else if (path[i]==dhb) { --b; }
        else if (path[i]==get_local_map().template beta<2>(dhb)) { ++b; }
      }
      else
      {
        if (path[i]==dha) { ++a; }
        else if (path[i]==get_local_map().template beta<2>(dha)) { --a; }
        else if (path[i]==dhb) { ++b; }
        else if (path[i]==get_local_map().template beta<2>(dhb)) { --b; }
      }
    }
  }

  Path_on_surface<Local_map>
  transform_original_path_into_quad_surface_for_torus
  (const Path_on_surface<Mesh>& path) const
  {
    CGAL_assertion(local_map_is_a_torus());

    Path_on_surface<Local_map> res(get_local_map());
    if (path.is_empty()) return res;

    Dart_const_handle cur;
    for (std::size_t i=0; i<path.length(); ++i)
    {
      if (!is_contracted(path[i]))// here flip doesn't matter
      {
        cur=get_first_dart_of_the_path(path[i], path.get_ith_flip(i), false);
        CGAL_assertion(!get_local_map().is_marked(cur, m_mark_perforated));
        while(cur!=get_second_dart_of_the_path(path[i], path.get_ith_flip(i), false))
        {
          res.push_back(cur, false, false);
          cur=get_local_map().template beta<1>(cur);
        }
      }
    }
    res.update_is_closed();
    CGAL_assertion(res.is_empty() || res.is_closed());
    CGAL_assertion(res.is_valid());
    return res;
  }

  internal::Path_on_surface_with_rle<Self>
  transform_original_path_into_quad_surface_with_rle
  (const Path_on_surface<Mesh>& path) const
  {
    internal::Path_on_surface_with_rle<Self>
      res(*this
#ifdef CGAL_PWRLE_TURN_V2
          , m_dart_ids
#endif //CGAL_PWRLE_TURN_V2
          );

    if (path.is_empty()) return res;

    for (std::size_t i=0; i<path.length(); ++i)
    {
      if (!is_contracted(path[i]))
      {
        res.push_back(get_first_dart_of_the_path
                      (path[i], path.get_ith_flip(i)), true);
        if (!edge_path_has_only_one_dart(path[i]))
        { res.push_back(get_second_dart_of_the_path
                        (path[i], path.get_ith_flip(i)), true); }
      }
    }
    res.update_is_closed();
    if (!res.is_empty())
    { res.merge_last_flat_with_next_if_possible(); }
    CGAL_assertion(res.is_closed() || res.is_empty());
    CGAL_assertion(res.is_valid());
    return res;
  }

  /// Mark the edge containing adart in the original map.
  void mark_original_edge(Original_dart_const_handle adart,
                          Original_size_type amark)
  {
    get_original_map().mark(adart, amark);
    if (!get_original_map().template is_free<2>(adart))
    { get_original_map().mark(get_original_map().template beta<2>(adart),
                              amark); }
  }

  /// Mark the edge containing adart in the reduced map (which is 2-closed)
  void mark_reduced_edge(Dart_const_handle adart, size_type amark)
  {
    get_local_map().mark(adart, amark);
    get_local_map().mark(get_local_map().template beta<2>(adart), amark);
  }

  /// Erase the edge given by adart (which belongs to the map m_local_map, i.e. the
  /// copy) from the associative array copy_to_origin, and erase the
  /// corresponding edge (which belongs to the map get_original_map()) from the
  /// array origin_to_copy
  void erase_edge_from_associative_arrays(Dart_handle adart,
                                          Origin_to_copy& origin_to_copy,
                                          Copy_to_origin & copy_to_origin)
  {
    if (!get_original_map().template is_free<2>(copy_to_origin[adart]))
    {
      origin_to_copy.erase(get_original_map().template beta<2>
                           (copy_to_origin[adart]));
      copy_to_origin.erase(get_local_map().template beta<2>(adart));
    }

    origin_to_copy.erase(copy_to_origin[adart]);
    copy_to_origin.erase(adart);
  }

  Original_dart_const_handle prev_in_boundary(Original_dart_const_handle d)
  {
    CGAL_assertion(get_original_map().template is_free<2>(d));
    Original_dart_const_handle res=get_original_map().template beta<0>(d);
    while(!get_original_map().template is_free<2>(res))
    { res=get_original_map().template beta<2, 0>(res); }
    return res;
  }

  /// Step 1) Transform m_local_map into an equivalent surface having only one
  /// vertex. All edges contracted during this step belong to the spanning
  /// tree T, and thus corresponding edges in get_original_map() are marked.

  /// Marks all darts belonging to T using a BFS
  void compute_T()
  {
    Original_dart_const_handle dh;
    auto grey=get_original_map().get_new_mark();
    std::queue<Original_dart_const_handle> queue;
    get_original_map().template mark_cell<0>
        (get_original_map().darts().begin(), grey);
    queue.push(get_original_map().darts().begin());

    while (!queue.empty())
    {
      dh=queue.front();
      queue.pop();
      for (auto it=get_original_map().template darts_of_cell<0>(dh).begin(),
             itend=get_original_map().template darts_of_cell<0>(dh).end();
               it!=itend; ++it)
      {
        if (!get_original_map().is_marked
            (get_original_map().template beta<1>(it), grey))
        {
          mark_original_edge(it, m_mark_T);
          get_original_map().template mark_cell<0>
              (get_original_map().template beta<1>(it), grey);
          queue.push(get_original_map().template beta<1>(it));
        }
      }
    }

    CGAL_assertion(get_original_map().is_whole_map_marked(grey));
    get_original_map().free_mark(grey);
  }

  /// Step 1) Create the copy of the input mesh, while contracting all edges of
  /// THE spanning tree. Initializes also m_path, the pair of darts of the edge
  /// in m_copy. This pair of edges will be updated later
  /// (in surface_simplification_in_one_face() and in surface_quadrangulate()).
  void surface_simplification_in_one_vertex
  (Origin_to_copy& origin_to_copy,Copy_to_origin& copy_to_origin)
  {
    m_paths.clear();
    compute_T();

    Dart_handle d1, d2;
    for (typename Original_map::Dart_range::const_iterator
         it=get_original_map().darts().begin(),
         itend=get_original_map().darts().end(); it!=itend; ++it)
    {
      if (!get_original_map().is_marked(it, m_mark_T))
      {
        if (get_original_map().template is_free<2>(it))
        {// case of a boundary
          d1=get_local_map().create_dart();
          d2=get_local_map().create_dart();
          get_local_map().template basic_link_beta_for_involution<2>(d1, d2);

          origin_to_copy[it]=d1;
          copy_to_origin[d1]=it;
          get_local_map().mark(d2, m_mark_perforated);
          if (get_original_map().is_perforated(it))
          { get_local_map().mark(d1, m_mark_perforated); }

          m_paths[it]=std::make_pair(d1, nullptr); // Initialize m_paths
        }
        else if (Original_dart_const_handle(it)<
                 get_original_map().template beta<2>(it))
        {
          d1=get_local_map().create_dart();
          d2=get_local_map().create_dart();
          get_local_map().template basic_link_beta_for_involution<2>(d1, d2);

          origin_to_copy[it]=d1;
          origin_to_copy[get_original_map().template beta<2>(it)]=d2;
          copy_to_origin[d1]=it;
          copy_to_origin[d2]=get_original_map().template beta<2>(it);

          if (get_original_map().is_perforated(it))
          { get_local_map().mark(d1, m_mark_perforated); }
          if (get_original_map().is_perforated
              (get_original_map().template beta<2>(it)))
          { get_local_map().mark(d2, m_mark_perforated); }

          m_paths[it]=std::make_pair(d1, nullptr); // Initialize m_paths
        }
      }
    }

    /// Now we only need to do the basic_link_beta_1
    Original_dart_const_handle dd1;
    for (typename Original_map::Dart_range::const_iterator
         it=get_original_map().darts().begin(),
         itend=get_original_map().darts().end(); it!=itend; ++it)
    {
      if (!get_original_map().is_marked(it, m_mark_T))
      {
        dd1=get_original_map().template beta<1>(it);
        while(get_original_map().is_marked(dd1, m_mark_T))
        { dd1=get_original_map().template beta<1>(dd1); }
        get_local_map().basic_link_beta_1(origin_to_copy[it],
                                            origin_to_copy[dd1]); // let's link both

        if (get_original_map().template is_free<2>(it))
        {
          dd1=prev_in_boundary(it);
          while(get_original_map().is_marked(dd1, m_mark_T))
          { dd1=prev_in_boundary(dd1); }
          get_local_map().basic_link_beta_1
              (get_local_map().template beta<2>(origin_to_copy[it]),
               get_local_map().template beta<2>(origin_to_copy[dd1]));
        }
      }
    }
  }

  /// Step 2) Transform the 2-map into an equivalent surface having only
  /// one face. All edges removed during this step belong to the
  /// dual spanning tree L (spanning tree of the dual 2-map).

  /// Marks all darts belonging to L using a BFS
  void compute_L(size_type toremove, Copy_to_origin& copy_to_origin)
  {
    Dart_handle dh;
    Dart_handle ddh;
    auto grey=get_local_map().get_new_mark();
    std::queue<Dart_handle> queue;

    for (auto it=get_local_map().darts().begin(),
         itend=get_local_map().darts().end(); it!=itend; ++it)
    {
      if (!get_local_map().is_marked(it, grey))
      {
        get_local_map().template mark_cell<2>(it, grey);
        if (!get_local_map().is_marked(it, m_mark_perforated))
        {
          queue.push(it);

          while (!queue.empty())
          {
            dh=queue.front(); // face(dh) is not perforated (and thus all its darts)
            queue.pop();
            ddh=dh;
            do
            {
              if (!get_local_map().is_marked
                  (get_local_map().template beta<2>(ddh), grey) &&
                  !get_local_map().is_marked
                  (get_local_map().template beta<2>(ddh), m_mark_perforated))
              {
                mark_reduced_edge(ddh, toremove);
                mark_original_edge(copy_to_origin[ddh], m_mark_L);
                get_local_map().template mark_cell<2>
                    (get_local_map().template beta<2>(ddh), grey);
                queue.push(get_local_map().template beta<2>(ddh));
              }
              ddh=get_local_map().template beta<1>(ddh);
            }
            while (dh!=ddh);
          }
        }
      }
    }

    CGAL_assertion(get_local_map().is_whole_map_marked(grey));
    get_local_map().free_mark(grey);
  }

  /// Update all length two paths, before edge removal. Edges that will be
  /// removed are marked with toremove mark. Dart mark as perforated are not
  /// updates, since they are associated only with their corresponding dart
  /// in the original map (and thus with only one dart and not two).
  void update_length_two_paths_before_edge_removals(size_type toremove,
                                                    Copy_to_origin& copy_to_origin)
  {
    // std::cout<<"************************************************"<<std::endl;
    Dart_handle initdart, curdart;
    Original_dart_const_handle d1, d2;
    for (typename Local_map::Dart_range::iterator
         it=get_local_map().darts().begin();
         it!=get_local_map().darts().end(); ++it)
    {
      initdart=it;
      if (!get_local_map().is_marked(initdart, toremove) &&
          !get_local_map().is_marked(initdart, m_mark_perforated))
      { // Here we are on a border edge of a "real" face (i.e. non perforated)
        curdart=get_local_map().template beta<0, 2>(initdart);
        while(get_local_map().is_marked(curdart, toremove))
        { // Here, all edges marked to remove are between two real faces.
          CGAL_assertion(copy_to_origin.count(curdart)==1);
          set_first_dart_of_the_path(copy_to_origin.find(curdart)->second,
                                     initdart);
          curdart=get_local_map().template beta<0, 2>(curdart);
        }

        if (!get_local_map().is_marked
            (get_local_map().template beta<2>(curdart), m_mark_perforated))
        {
          d1=copy_to_origin.find(get_local_map().template beta<2>(curdart))->second;
          if (get_original_map().template is_free<2>(d1) ||
              d1<get_original_map().template beta<2>(d1))
          { m_paths[d1].second=initdart; }
        }
      }
    }

#ifndef NDEBUG
    for (auto it=m_paths.begin(), itend=m_paths.end(); it!=itend; ++it)
    {
      CGAL_assertion(!get_local_map().is_marked(it->second.first, toremove));
      CGAL_assertion(it->second.second==nullptr ||
                     !get_local_map().is_marked(it->second.second, toremove));
    }
#endif
  }

  /// Remove all loops after having merge all real faces in one.
  /// Erase from m_paths all darts linked with such a dart.
  void remove_non_perforated_loops(Copy_to_origin& copy_to_origin)
  {
    get_local_map().set_automatic_attributes_management(false);

    Original_dart_const_handle origin_dart;
    for (typename Local_map::Dart_range::iterator
         it=get_local_map().darts().begin(),
         itend=get_local_map().darts().end(); it!=itend; ++it)
    {
      if (get_local_map().is_dart_used(it) &&
          !get_local_map().is_marked(it, m_mark_perforated) &&
          get_local_map().template beta<1>(it)==it)
        { // A non perforated loop.
          origin_dart=copy_to_origin[it];
          if (!get_original_map().template is_free<2>(origin_dart) &&
              origin_dart>get_original_map().template beta<2>(origin_dart))
          { origin_dart=get_original_map().template beta<2>(origin_dart); }

          mark_original_edge(origin_dart, m_mark_T);

          if (get_local_map().is_marked(get_local_map().template beta<2>(it),
                                  m_mark_perforated))
          { get_local_map().template mark_cell<2>(it, m_mark_perforated); }

          get_local_map().template remove_cell<1>(it);
        }
      }

    get_local_map().set_automatic_attributes_management(true);

    /// We remove from m_paths all associations having a removed dart
    for (auto it=m_paths.begin(), itend=m_paths.end(); it!=itend; )
    {
      if (!get_local_map().is_dart_used(it->second.first))
      {
        mark_original_edge(it->first, m_mark_T);
        it=m_paths.erase(it);
      }
      else if (it->second.second!=nullptr &&
               !get_local_map().is_dart_used(it->second.second))
      {
        mark_original_edge(it->first, m_mark_T);
        it=m_paths.erase(it);
      }
      else { ++it; }
    }
  }

  /// Simplify the reduced map by merging all adjacent non perforated faces.
  void surface_simplification_in_one_face
  (Origin_to_copy& origin_to_copy, Copy_to_origin& copy_to_origin)
  {
    get_local_map().set_automatic_attributes_management(false);

    size_type toremove=get_local_map().get_new_mark();
    compute_L(toremove, copy_to_origin);

    if (get_local_map().number_of_marked_darts(toremove)==
        get_local_map().number_of_darts())
    { // Case of sphere; all darts are removed.
      m_paths.clear();
      get_local_map().clear();
    }
    else
    {
      // Update the length two paths before to remove the edges.
      update_length_two_paths_before_edge_removals(toremove, copy_to_origin);

      // Then remove all the edges marekd to remove.
      for (typename Local_map::Dart_range::iterator
             it=get_local_map().darts().begin(),
           itend=get_local_map().darts().end(); it!=itend; ++it)
      {
        if (get_local_map().is_dart_used(it) &&
            get_local_map().is_marked(it, toremove))
        {
          CGAL_assertion(get_local_map().is_marked
                         (get_local_map().template beta<2>(it), toremove));
          erase_edge_from_associative_arrays(it, origin_to_copy, copy_to_origin);
          get_local_map().template remove_cell<1>(it);
        }
      }
    }

    get_local_map().set_automatic_attributes_management(true);
    get_local_map().free_mark(toremove);

#ifndef NDEBUG
    for (auto it=m_paths.begin(), itend=m_paths.end(); it!=itend; ++it)
    {
      CGAL_assertion(get_local_map().is_dart_used(it->second.first));
      CGAL_assertion(it->second.second==nullptr ||
                     get_local_map().is_dart_used(it->second.second));
    }
#endif
  }

  /// Step 4) quadrangulate the surface.
  void surface_quadrangulate()
  {
    // Here the map has only one vertex and one face if we have a closed surface,
    // and maybe several faces if the surface has boundaries
    size_type oldedges=get_local_map().get_new_mark();
    get_local_map().negate_mark(oldedges); // now all edges are marked

    // 1) We insert a vertex in each face which is not perforated.
    //    New edges created by the operation are not marked oldedges.
    size_type treated=get_local_map().get_new_mark();

    for (typename Local_map::Dart_range::iterator
         it=get_local_map().darts().begin();
         it!=get_local_map().darts().end(); ++it)
    {
      if (!get_local_map().is_marked(it, treated))
      {
        get_local_map().template mark_cell<2>(it, treated);
        if (get_local_map().is_marked(it, oldedges) &&
            !get_local_map().is_marked(it, m_mark_perforated))
        {
          CGAL_assertion(get_local_map().template beta<1>(it)!=it);
          get_local_map().negate_mark(treated); // To mark new darts treated
          get_local_map().insert_cell_0_in_cell_2(it);
          get_local_map().negate_mark(treated);
        }
      }
    }
    CGAL_assertion(get_local_map().is_whole_map_marked(treated));
    get_local_map().free_mark(treated);

#ifdef NDEBUG
    for (typename Local_map::Dart_range::iterator
         it=get_local_map().darts().begin();
         it!=get_local_map().darts().end(); ++it)
    {
      if (!get_local_map().is_marked(it, m_mark_perforated))
      {
        CGAL_assertion(get_local_map().template beta<1>(it)!=it);
        CGAL_assertion((get_local_map().template beta<1,1,1>(it)==it));
      }
    }
#endif

    // 2) We update the pair of darts
    for (typename TPaths::iterator itp=m_paths.begin(), itpend=m_paths.end();
         itp!=itpend; ++itp)
    {
      std::pair<Dart_handle, Dart_handle>& p=itp->second;
      if (p.second!=nullptr)
      { // Edge between two real faces, removed during the quadrangulation
        CGAL_assertion(!get_local_map().is_marked(p.first, m_mark_perforated));
        CGAL_assertion(!get_local_map().is_marked(p.second, m_mark_perforated));
        CGAL_assertion(get_local_map().template beta<1>(p.first)!=p.first);
        CGAL_assertion(get_local_map().template beta<1>(p.second)!=p.second);
        CGAL_assertion(get_local_map().is_marked(p.first, oldedges));
        CGAL_assertion(get_local_map().is_marked(p.second, oldedges));
        CGAL_assertion((get_local_map().template beta<1,1,1>(p.first)==p.first));
        CGAL_assertion((get_local_map().template beta<1,1,1>(p.second)==p.second));
        p.first=get_local_map().template beta<0, 2>(p.first);
        p.second=get_local_map().template beta<0>(p.second);
        CGAL_assertion(!get_local_map().is_marked(p.first, oldedges));
        CGAL_assertion(!get_local_map().is_marked(p.second, oldedges));
      }
      else if (!get_local_map().is_marked(p.first, m_mark_perforated))
      { // Edge between a real face and a perforated one, not updated during update_length_two_paths_before_edge_removals
        CGAL_assertion(get_local_map().is_marked(p.first, oldedges));
        CGAL_assertion(get_local_map().template beta<1>(p.first)!=p.first);
        CGAL_assertion((get_local_map().template beta<1,1,1>(p.first)==p.first));
        p.second=get_local_map().template beta<1, 2>(p.first);
        p.first=get_local_map().template beta<0, 2>(p.first);
        CGAL_assertion(!get_local_map().is_marked(p.first, oldedges));
        CGAL_assertion(!get_local_map().is_marked(p.second, oldedges));
      }
      else if (!get_local_map().is_marked
               (get_local_map().template beta<2>(p.first), m_mark_perforated))
      { // Edge between a perforated face and a real one, not updated during update_length_two_paths_before_edge_removals
        CGAL_assertion(get_local_map().is_marked(p.first, oldedges));
        CGAL_assertion((get_local_map().template beta<2,1,1,1>(p.first)==
                        get_local_map().template beta<2>(p.first)));
        CGAL_assertion((get_local_map().template beta<2,1>(p.first)!=
            get_local_map().template beta<2>(p.first)));
        p.first=get_local_map().template beta<2, 1>(p.first);
        p.second=get_local_map().template beta<1>(p.first);
        CGAL_assertion(!get_local_map().is_marked(p.first, oldedges));
        CGAL_assertion(!get_local_map().is_marked(p.second, oldedges));
      }
    }

#ifdef NDEBUG
    for (auto it=m_paths.begin(), itend=m_paths.end(); it!=itend; ++it)
    {
      if (it->second.second==nullptr)
      {
        CGAL_assertion(get_local_map().is_marked(it->second.first, oldedges));
        CGAL_assertion(get_local_map().is_marked(it->second.first, m_mark_perforated) &&
                       get_local_map().is_marked
                       (get_local_map().template beta<2>(it->second.first), m_mark_perforated));
      }
      else
      {
        CGAL_assertion(!get_local_map().is_marked(it->second.first, oldedges));
        CGAL_assertion(!get_local_map().is_marked(it->second.second, oldedges));
      }
    }
#endif

    // 3) We remove all the old edges, and extend the perforated faces
    //    (i.e. we remove edges between real and perforated faces).
    for (typename Local_map::Dart_range::iterator
         it=get_local_map().darts().begin(),
         itend=get_local_map().darts().end(); it!=itend; ++it)
    {
      if (get_local_map().is_dart_used(it) &&
          get_local_map().is_marked(it, oldedges) &&
          !get_local_map().is_marked(it, m_mark_perforated))
      {
        get_local_map().unmark(it, oldedges);
        if (get_local_map().is_marked
            (get_local_map().template beta<2>(it), m_mark_perforated))
        { get_local_map().template mark_cell<2>(it, m_mark_perforated); }
        get_local_map().template remove_cell<1>(it);
      }
      else { get_local_map().unmark(it, oldedges); }
    }

    CGAL_assertion(get_local_map().is_whole_map_unmarked(oldedges));
    get_local_map().free_mark(oldedges);

#ifdef NDEBUG
    for (auto it=m_paths.begin(), itend=m_paths.end(); it!=itend; ++it)
    {
      CGAL_assertion(get_local_map().is_dart_used(it->second.first));
      CGAL_assertion(it->second.second==nullptr ||
                     get_local_map().is_dart_used(it->second.second));
    }
#endif
  }

#if defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
  // Initialize all vertices attributes.
  // For all vertex v, the attribute of v is set to deg(v) if there is a
  // hole around v and -deg(v) if there are no holes around v.
  void initialize_vertices_attributes()
  {
    int deg;
    bool hole_detected;
    Dart_handle dh;
    for (auto it=get_local_map().darts().begin(),
         itend=get_local_map().darts().end(); it!=itend; ++it)
    {
      if (get_local_map().template attribute<0>(it)==NULL)
      {
        // We count the degree, while testing if there is a perforated face
        // incident to the vertex.
        deg=0;
        hole_detected=false;
        dh=it;
        do
        {
          ++deg;
          if (get_local_map().is_marked(dh, m_mark_perforated))
          { hole_detected=true; }
          dh=get_local_map().template beta<2, 1>(dh);
        }
        while(dh!=it);

        // Then we set the vertex attribute to deg id there is a
        // hole and -deg otherwise
        if (!hole_detected) { deg=-deg; }

        // We create the 0-attribute
        get_local_map().template set_attribute<0>
            (it, get_local_map().template create_attribute<0>(deg));
      }
    }
  }

  void initialize_ids()
  {
    std::size_t id;
    Dart_handle dh, ddh;
    auto treated=get_local_map().get_new_mark();

    for (auto it=get_local_map().darts().begin(),
         itend=get_local_map().darts().end(); it!=itend; ++it)
    {
      if (!get_local_map().is_marked(it, treated))
      {
        id=0;
        if (get_local_map().template info<0>(it)<0)
        {// there are no holes around this vertex
         // we just set the ids from 0 to deg(v)-1
          dh=it;
          do
          {
#ifdef CGAL_PWRLE_TURN_V2
            m_dart_ids[dh]=id;
#else // CGAL_PWRLE_TURN_V2
            // this is for the turn V3
            get_local_map().info(dh)=id;
          #endif // CGAL_PWRLE_TURN_V2
            ++id;
            get_local_map().mark(dh, treated);
            dh=get_local_map().template beta<2, 1>(dh);
          }
          while(dh!=it);
        }
        else
        {// there is at least a hole around the vertex
         // we set the 0-th dart just after a hole
         // then we add 1 to the next dart if we dont cross a hole
         // and we add deg(v)+1 if we cross a hole
          dh=it;
          while(!get_local_map().is_marked(dh, m_mark_perforated))
          {
            dh=get_local_map().template beta<2, 1>(dh);
          }
          // now dh is right after a hole
          ddh=dh;
          do
          {
#ifdef CGAL_PWRLE_TURN_V2
            m_dart_ids[ddh]=id;
#else // CGAL_PWRLE_TURN_V2
            // this is for the turn V3
            get_local_map().info(ddh)=id;
          #endif // CGAL_PWRLE_TURN_V2
            if (get_local_map().is_marked(get_local_map().template beta<2>(ddh), m_mark_perforated))
            { id+=get_local_map().template info<0>(ddh)+1; }
            else
            { id+=1; }
            get_local_map().mark(ddh, treated);
            ddh=get_local_map().template beta<2, 1>(ddh);
          }
          while(ddh!=dh);
        }
      }
    }
    get_local_map().free_mark(treated);
  }

  std::size_t get_dart_id(Dart_const_handle dh) const
  {
#ifdef CGAL_PWRLE_TURN_V2
    return m_dart_ids.at(dh);
#else //  CGAL_PWRLE_TURN_V2
    // this is for the turn V3
    return get_local_map().info(dh);
#endif // CGAL_PWRLE_TURN_V2
    std::cerr<<"Error: impossible to get dart id without method V2 or V3."<<std::endl;
    return (std::numeric_limits<std::size_t>::max)();
  }

  /// @return the positive turn given two darts using their ids (unsed for CGAL_PWRLE_TURN_V2 and V3)
  std::size_t compute_positive_turn_given_ids(Dart_const_handle dh1,
                                              Dart_const_handle dh2) const
  {
    if (get_local_map().template info<0>(dh1)<0)
    {// there is no hole around dh1 and dh2
      if (get_dart_id(dh1)<=get_dart_id(dh2))
      {
        return get_dart_id(dh2)-get_dart_id(dh1);
      }
      // here we have to add the degree (i.e. substract the vertex info)
      return get_dart_id(dh2)-get_local_map().template info<0>(dh1)-get_dart_id(dh1);
    }
    // here we know there is a hole just before the dart 0 (plus maybe other ones)
    if (get_dart_id(dh1)>get_dart_id(dh2) || // we crossed dart 0, so we crossed a hole
        get_dart_id(dh2)-get_dart_id(dh1)>
        static_cast<std::size_t>(get_local_map().template info<0>(dh1))) // the gap is more than the degree, so we crossed a hole
    {// so we return an "infinite" value
      return (std::numeric_limits<std::size_t>::max)();
    }
    return get_dart_id(dh2)-get_dart_id(dh1);
  }

  /// @return the negative turn given two darts using their ids (unsed for CGAL_PWRLE_TURN_V2 and V3)
  std::size_t compute_negative_turn_given_ids(Dart_const_handle dh1,
                                              Dart_const_handle dh2) const
  {
    if (get_local_map().template info<0>(dh1)<0)
    {// there is no hole around dh1 and dh2
      if (get_dart_id(dh1)>=get_dart_id(dh2))
      {
        return get_dart_id(dh1)-get_dart_id(dh2);
      }
      // here we have to add the degree (i.e. substract the vertex info)
      return get_dart_id(dh1)-get_local_map().template info<0>(dh1)-get_dart_id(dh2);
    }
    // here we know there is a hole just before the dart 0 (plus maybe other ones)
    if (get_dart_id(dh1)<get_dart_id(dh2) || // we crossed dart 0, so we crossed a hole
        get_dart_id(dh1)-get_dart_id(dh2)>
        static_cast<std::size_t>(get_local_map().template info<0>(dh1))) // the gap is more than the degree, so we crossed a hole
    {// so we return an "infinite" value
      return (std::numeric_limits<std::size_t>::max)();
    }
    return get_dart_id(dh1)-get_dart_id(dh2);
  }
#endif // defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)

  /// @return true iff the edge containing adart is associated with a path.
  ///         (used for debug purpose because we are suppose to be able to
  ///          test this by using directly the mark m_mark_T).
  bool is_edge_has_path(Original_dart_const_handle adart) const
  {
    if (get_original_map().template is_free<2>(adart) ||
        adart<get_original_map().template beta<2>(adart))
    { return m_paths.find(adart)!=m_paths.end(); }
    return m_paths.find(get_original_map().template beta<2>(adart))!=m_paths.end();
  }

  /// @return true iff the edge containing adart is associated with a path
  ///         of only 1 dart (case of an edge bewteen two perforated faces)
  bool edge_path_has_only_one_dart(Original_dart_const_handle adart) const
  {
    return
        (get_original_map().is_perforated(adart) &&
         (get_original_map().template is_free<2>(adart) ||
          get_original_map().is_perforated
          (get_original_map().template beta<2>(adart))));
  }

  /// @return the pair of darts associated with the edge containing adart
  ///         in get_original_map().
  /// @pre the edge containing adart must not belong to T.
  std::pair<Dart_handle, Dart_handle>& get_pair_of_darts
  (Original_dart_const_handle adart)
  {
    CGAL_assertion(!is_contracted(adart));
    CGAL_assertion(is_edge_has_path(adart));

    if (get_original_map().template is_free<2>(adart) ||
        adart<get_original_map().template beta<2>(adart))
    return m_paths.find(get_original_map().template beta<2>(adart))->second;
  }

  /** @return the first dart of the length 2 path associated with adart.
   *  @param adart a dart in the original map
   *  @param flip to return the first dart of beta2(adart)
   *  @param isquadrangulation true if the minimal map is a quadrangulation,
   *                           false otherwise (special case of a torus)
   *  There is a special case for an edge between two perforated faces. In this
   *  case there is only one dart in the path.
   */
  Dart_const_handle get_first_dart_of_the_path(Original_dart_const_handle adart,
                                               bool flip=false,
                                               bool isquadrangulation=true) const
  {
    CGAL_assertion(is_edge_has_path(adart));

    if (isquadrangulation && get_original_map().is_perforated(adart) &&
        (get_original_map().template is_free<2>(adart) ||
         get_original_map().is_perforated
         (get_original_map().template beta<2>(adart))))
    {
      if (get_original_map().template is_free<2>(adart) ||
          adart<get_original_map().template beta<2>(adart))
      {
        const std::pair<Dart_handle, Dart_handle>&
          p=m_paths.find(adart)->second;
        return flip?get_local_map().template beta<2>(p.first):p.first;
      }
      else
      {
        const std::pair<Dart_handle, Dart_handle>&
            p=m_paths.find(get_original_map().template beta<2>(adart))->second;
        return flip?p.first:get_local_map().template beta<2>(p.first);
      }
    }

    CGAL_assertion(!get_original_map().is_perforated(adart) ||
                   (!get_original_map().template is_free<2>(adart) &&
                    !get_original_map().is_perforated
                    (get_original_map().template beta<2>(adart))));

    if (get_original_map().template is_free<2>(adart) ||
        adart<get_original_map().template beta<2>(adart))
    {
      const std::pair<Dart_handle, Dart_handle>&
        p=m_paths.find(adart)->second;
      return flip?
            (isquadrangulation?get_local_map().template beta<2>(p.second):p.second):
            p.first;
    }

    const std::pair<Dart_handle, Dart_handle>&
        p=m_paths.find(get_original_map().template beta<2>(adart))->second;

    return flip?
          p.first:
          (isquadrangulation?get_local_map().template beta<2>(p.second):p.second);
  }

  /** @return the second dart of the length 2 path associated with adart.
   *  @param adart a dart in the original map
   *  @param flip to return the first dart of beta2(adart)
   *  @param isquadrangulation true if the minimal map is a quadrangulation,
   *                           false otherwise (special case of a torus)
   *  @pre edge containing adart should be incident to at least one non
   *       perforated face.
   */
  Dart_const_handle get_second_dart_of_the_path(Original_dart_const_handle adart,
                                                bool flip=false,
                                                bool isquadrangulation=true) const
  {
    CGAL_assertion(is_edge_has_path(adart));
    CGAL_assertion(!get_original_map().is_perforated(adart) ||
                   (!get_original_map().template is_free<2>(adart) &&
                    !get_original_map().is_perforated
                    (get_original_map().template beta<2>(adart))));

    if (get_original_map().template is_free<2>(adart) ||
        adart<get_original_map().template beta<2>(adart))
    {
      const std::pair<Dart_handle, Dart_handle>&
          p=m_paths.find(adart)->second;
      return flip?
            (isquadrangulation?get_local_map().template beta<2>(p.first):p.first):
            p.second;
    }

    const std::pair<Dart_handle, Dart_handle>&
        p=m_paths.find(get_original_map().template beta<2>(adart))->second;
    return flip?
          p.second:
          (isquadrangulation?get_local_map().template beta<2>(p.first):p.first);
  }

  /// @return the first dart of the path, direct version without taking into
  /// account possible flip and the two case of torus/quadrangulation.
  Dart_handle get_first_dart_of_the_path_direct(Original_dart_const_handle adart) const
  {
    CGAL_assertion(is_edge_has_path(adart));

    if (get_original_map().template is_free<2>(adart) ||
        adart<get_original_map().template beta<2>(adart))
    { return m_paths.find(adart)->second.first; }
    return m_paths.find(get_original_map().template beta<2>(adart))->second.second;
  }

  void set_first_dart_of_the_path(Original_dart_const_handle adart,
                                  Dart_handle d)
  {
    CGAL_assertion(is_edge_has_path(adart));

    if (get_original_map().template is_free<2>(adart) ||
        adart<get_original_map().template beta<2>(adart))
    { m_paths.find(adart)->second.first=d; }
    else
    { m_paths.find(get_original_map().template beta<2>(adart))->second.second=d; }
  }

  void set_second_dart_of_the_path(Original_dart_const_handle adart,
                                   Dart_handle d)
  {
    CGAL_assertion(is_edge_has_path(adart));

    if (get_original_map().template is_free<2>(adart) ||
        adart<get_original_map().template beta<2>(adart))
    { m_paths.find(adart)->second.second=d; }
    else { m_paths.find(adart)->second.first=d; }
  }

  bool local_map_is_a_torus() const
  {
    if (get_local_map().number_of_darts()!=4)
    { return false; }
    return (get_local_map().number_of_marked_darts(m_mark_perforated)==0);
  }

  /// @return true iff the perforated faces are correctly marked
  /// (i.e. either fully marked or fully unmarked).
  bool perforated_faces_correctly_marked() const
  {
    for (typename Local_map::Dart_range::const_iterator
         it=get_local_map().darts().begin(),
         itend=get_local_map().darts().end(); it!=itend; ++it)
    {
      CGAL_assertion(get_local_map().is_without_boundary(1));
      if (get_local_map().is_marked(it, m_mark_perforated))
      {
        if (!get_local_map().template is_whole_cell_marked<2>
            (it, m_mark_perforated))
        {
          std::cout<<"[ERROR] perforated_faces_correctly_marked: it is marked, "
                     <<"but its face is not fully marked."<<std::endl;
          return false;
        }
      }
      else
      {
        if (!get_local_map().template is_whole_cell_unmarked<2>
            (it, m_mark_perforated))
        {
          std::cout<<"[ERROR] perforated_faces_correctly_marked: it is not marked, "
                     <<"but its face is not fully unmarked."<<std::endl;
          return false;
        }
      }
    }
    return true;
  }

  /// Test if m_paths are valid, i.e.:
  /// 1) all the darts of get_original_map() that do not belong to T are
  ///    associated with a pair of darts;
  /// 2) all the darts of m_paths belong to m_local_map;
  /// 3) the origin of the second dart of the pair is the extremity of the
  ///    first dart.
  /// 4) all the darts of m_local_map are not free (both for beta 1 and 2)
  /// 5) The two darts in a pair are different
  bool are_paths_valid() const
  {
    if (m_paths.empty()) { return true; }

    bool res=true;
    for (auto it=get_original_map().darts().begin(),
           itend=get_original_map().darts().end(); it!=itend; ++it)
    {

      if (!is_contracted(it))
      {
        if (!is_edge_has_path(it))
        {
          std::cout<<"ERROR: an edge that is not contracted "
                   <<"has no associated path."
                   <<"BTW is_marked(it)="
                   <<get_original_map().is_marked(it, m_mark_T)<<std::endl;
          res=false;
        }
      }
      else
      {
        if (is_edge_has_path(it))
        {
          std::cout<<"ERROR: an edge that is contracted"
                   <<" has an associated path."
                   <<"BTW is_marked(it)="
                   <<get_original_map().is_marked(it, m_mark_T)<<std::endl;
          res=false;
        }
      }
    }

    for (auto it=get_local_map().darts().begin(),
           itend=get_local_map().darts().end(); it!=itend; ++it)
    {
      if (get_local_map().is_free(it, 1))
      {
        std::cout<<"ERROR: a dart of the quandrangulated map is 1-free"
                 <<std::endl;
        res=false;
      }
      if (get_local_map().is_free(it, 2))
      {
        std::cout<<"ERROR: a dart of the quandrangulated map is 2-free"
                 <<std::endl;
        res=false;
      }
    }

    for (auto it=m_paths.begin(); it!=m_paths.end(); ++it)
    {
      if (!get_local_map().is_dart_used(it->second.first))
      {
        std::cout<<"ERROR: first dart in m_paths does not exist anymore in m_local_map."
                 <<std::endl;
        res=false;
      }
      else if (!get_local_map().darts().owns(it->second.first))
      {
        std::cout<<"ERROR: first dart in m_paths does not belong to m_local_map."
                 <<std::endl;
        res=false;
      }

      if (!edge_path_has_only_one_dart(it->first))
      {
        if (!get_local_map().is_dart_used(it->second.second))
        {
          std::cout<<"ERROR: second dart in m_paths does not exist anymore in m_local_map."
                  <<std::endl;
          res=false;
        }
        else if (!get_local_map().darts().owns(it->second.second))
        {
          std::cout<<"ERROR: second dart in m_paths does not belong to m_local_map."
                  <<std::endl;
          res=false;
        }
      }

      if (it->second.first==it->second.second)
      {
        std::cout<<"ERROR: two darts in the same pair are equal."
                 <<std::endl;
        res=false;
      }
    }

    for (auto it=get_original_map().darts().begin(),
           itend=get_original_map().darts().end(); it!=itend; ++it)
    {
      if (!is_contracted(it))
      {
        Dart_const_handle d1=get_first_dart_of_the_path(it);
        if (d1==nullptr)
        {
          std::cout<<"ERROR: an edge is associated with a null dart in m_paths"
                   <<" for its first dart."<<std::endl;
          res=false;
        }
        else if (!edge_path_has_only_one_dart(it))
        {
          Dart_const_handle d2=get_second_dart_of_the_path(it);
          if (d2==nullptr)
          {
            std::cout<<"ERROR: an edge is associated with a null dart in m_paths"
                     <<" for its second dart."<<std::endl;
            res=false;
          }
          else
          {
            Dart_const_handle dd1=get_local_map().other_extremity(d1);
            CGAL_assertion(dd1!=NULL);
            if (!CGAL::belong_to_same_cell<Local_map,0>
                (get_local_map(), dd1, d2))
            {
              std::cout<<"ERROR: the two darts in a path are not consecutive."
                      <<std::endl;
              res=false;
            }
          }
        }
      }
    }

    return res;
  }

protected:
  /// The original map (the mesh seen as a 2-map)
  const typename Get_map<Mesh, Mesh>::storage_type m_original_map;
  Local_map m_local_map; /// the reduced map
  TPaths m_paths; /// Pair of edges associated with each edge of get_original_map()
                  /// (except the edges that belong to the spanning tree T).
  Original_size_type m_mark_T;   /// mark each edge of get_original_map() that belong to the spanning tree T
  Original_size_type m_mark_L;   /// mark each edge of get_original_map() that belong to the dual spanning tree L
  size_type m_mark_perforated; /// mark each edge of m_local_map that bounds a hole

#ifdef CGAL_PWRLE_TURN_V2
  TDartIds m_dart_ids; /// Ids of each dart of the transformed map, between 0 and n-1 (n being the number of darts)
                       /// so that darts between 0...(n/2)-1 belong to the same vertex and
                       /// d1=beta<1, 2>(d0), d2=beta<1, 2>(d1)...
                       /// The same for darts between n/2...n-1 for the second vertex
                       /// Thanks to these ids, we can compute in constant time the positive and
                       /// negative turns between two consecutive darts
#endif // CGAL_PWRLE_TURN_V2
};

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif // CGAL_MINIMAL_QUADRANGULATION_H //
// EOF //
