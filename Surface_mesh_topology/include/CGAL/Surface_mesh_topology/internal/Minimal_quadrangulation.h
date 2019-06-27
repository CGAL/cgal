// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_MINIMAL_QUADRANGULATION_H
#define CGAL_MINIMAL_QUADRANGULATION_H 1

// Should be defined before to include Path_on_surface_with_rle.h
// If nothing is defined, use V1
// #define CGAL_PWRLE_TURN_V1  // Compute turns by turning (method of CMap)
// #define CGAL_PWRLE_TURN_V2  // Compute turns by using an id of darts, given by an hash-table (built and given by Curves_on_surface_topology)
#define CGAL_PWRLE_TURN_V3  // Compute turns by using an id of darts, associated in Info of Darts (build by Curves_on_surface_topology)

#include <CGAL/license/Surface_mesh_topology.h>

#include <CGAL/Union_find.h>
#include <CGAL/Random.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Surface_mesh_topology/internal/Path_on_surface_with_rle.h>
#include <CGAL/Surface_mesh_topology/internal/Path_generators.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <CGAL/Timer.h>
#include <boost/unordered_map.hpp>
#include <stack>
#include <iostream>

namespace CGAL {
namespace Surface_mesh_topology {
namespace internal {

struct CMap_for_minimal_quadrangulation_items
{
  template <class CMap>
  struct Dart_wrapper
  {
#ifdef CGAL_PWRLE_TURN_V3
    typedef std::size_t Dart_info;
#endif // CGAL_PWRLE_TURN_V3
    typedef CGAL::cpp11::tuple<> Attributes;
  };
};
typedef CGAL::Combinatorial_map<2, CMap_for_minimal_quadrangulation_items>
CMap_for_minimal_quadrangulation;

template<typename Mesh>
class Minimal_quadrangulation
{
public:
  typedef typename CMap_for_minimal_quadrangulation::Dart_handle
  Dart_handle;
  typedef typename CMap_for_minimal_quadrangulation::Dart_const_handle
  Dart_const_handle;
  typedef CGAL::Union_find<Dart_handle> UFTree;
  typedef typename UFTree::handle UFTree_handle;

  typedef typename Get_map<Mesh, Mesh>::type Map;

  typedef CGAL::Union_find<typename Map::Dart_const_handle> UFTree2;
  typedef typename UFTree2::handle UFTree_handle2;

  // Associate each dart of the original map, not removed, a pair of darts in the
  // reduced map.
  typedef boost::unordered_map<typename Map::Dart_const_handle,
                               std::pair<Dart_const_handle, Dart_const_handle> > TPaths;

#ifdef CGAL_PWRLE_TURN_V2
  typedef boost::unordered_map<Dart_const_handle, std::size_t> TDartIds;
#endif //CGAL_PWRLE_TURN_V2

  Minimal_quadrangulation(const Mesh& amap, bool display_time=false) :
    m_original_map(amap)
  {
    if (!m_map.is_without_boundary(1))
    {
      std::cerr<<"ERROR: the given amap has 1-boundaries; "
               <<"such a surface is not possible to process here."
               <<std::endl;
      return;
    }
    if (!m_map.is_without_boundary(2))
    {
      std::cerr<<"ERROR: the given amap has 2-boundaries; "
               <<"which are not yet considered (but this will be done later)."
               <<std::endl;
      return;
    }
 
    CGAL::Timer t, t2;
    if (display_time)
    { t.start(); t2.start(); }

    // The mapping between darts of the original map into the copied map.
    boost::unordered_map<typename Map::Dart_const_handle, Dart_handle>
      origin_to_copy;

    // The mapping between darts of the copy into darts of the original map.
    boost::unordered_map<Dart_handle, typename Map::Dart_const_handle>
      copy_to_origin;

    // We copy the original map, while keeping mappings between darts.
    // m_map.copy(m_original_map, &origin_to_copy, &copy_to_origin);

    if (display_time)
    {
      t2.stop();
      std::cout<<"[TIME] Copy map: "<<t2.time()<<" seconds"<<std::endl;

      t2.reset(); t2.start();
    }

    // We reserve the two marks (used to mark darts in m_original_map that
    // belong to T or to L)
    m_mark_T=m_original_map.get_new_mark();
    m_mark_L=m_original_map.get_new_mark();

    /* std::cout<<"Number of darts in m_map: "<<m_map.number_of_darts()
       <<"; number of darts in origin_to_copy: "<<origin_to_copy.size()
       <<"; number of darts in copy_to_origin: "<<copy_to_origin.size()
       <<std::endl; */

    // 1) We simplify m_map in a surface with only one vertex
    surface_simplification_in_one_vertex(origin_to_copy, copy_to_origin);

    if (display_time)
    {
      t2.stop();
      std::cout<<"[TIME] Simplification in one vertex: "
               <<t2.time()<<" seconds"<<std::endl;

      t2.reset(); t2.start();
    }

#ifdef CGAL_TRACE_CMAP_TOOLS
    std::cout<<"All non loop contracted: ";
    m_map.display_characteristics(std::cout) << ", valid=" 
                                             << m_map.is_valid() << std::endl;
    /* std::cout<<"Number of darts in m_map: "<<m_map.number_of_darts()
       <<"; number of darts in origin_to_copy: "<<origin_to_copy.size()
       <<"; number of darts in copy_to_origin: "<<copy_to_origin.size()
       <<std::endl; */
#endif

    // 2) Now we compute each length two path associated with each edge that does
    // not belong to the spanning tree (which are thus all the survival edges).
    compute_length_two_paths(origin_to_copy);

    if (display_time)
    {
      t2.stop();
      std::cout<<"[TIME] Computation of length two paths: "
               <<t2.time()<<" seconds"<<std::endl;

      t2.reset(); t2.start();
    }

    /* std::cout<<"Number of darts in m_map: "<<m_map.number_of_darts()
       <<"; number of darts in origin_to_copy: "<<origin_to_copy.size()
       <<"; number of darts in copy_to_origin: "<<copy_to_origin.size()
       <<std::endl; */

    /* std::cout<<"Paths are all valid 1 ? "<<(are_paths_valid()?"YES":"NO")
       <<std::endl; */

    // 3) We simplify m_map in a surface with only one face
    surface_simplification_in_one_face(origin_to_copy, copy_to_origin);

    if (display_time)
    {
      t2.stop();
      std::cout<<"[TIME] Simplification in one face: "
               <<t2.time()<<" seconds"<<std::endl;

      t2.reset(); t2.start();
    }

#ifdef CGAL_TRACE_CMAP_TOOLS
    std::cout<<"All faces merges: ";
    m_map.display_characteristics(std::cout) << ", valid=" 
                                             << m_map.is_valid() << std::endl;

    /*      std::cout<<"Number of darts in m_map: "<<m_map.number_of_darts()
            <<"; number of darts in origin_to_copy: "<<origin_to_copy.size()
            <<"; number of darts in copy_to_origin: "<<copy_to_origin.size()
            <<std::endl; */
#endif

    if (!m_map.is_empty()) // m_map is_empty if the surface is a sphere
    {
      // 4) And we quadrangulate the face, except for the torus surfaces
      if (m_map.darts().size()!=4)
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

      // Now we label all the darts of the reduced map, to allow the computation
      // of turns in constant time: only for methods V2 and V3
#if defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
      CGAL_assertion(m_map.number_of_darts()%2==0);
      Dart_handle dh1=m_map.darts().begin();
      Dart_handle dh2=m_map.template beta<2>(dh1);
      std::size_t id=0;
      for(; dh1!=NULL; dh1=(dh1==dh2?NULL:dh2)) // We have two vertices to process
      {
        Dart_handle cur_dh=dh1;
        do
        {
#ifdef CGAL_PWRLE_TURN_V2
          m_dart_ids[cur_dh]=id++;
#else //  CGAL_PWRLE_TURN_V2
          // Here we use CGAL_PWRLE_TURN_V3
          m_map.info(cur_dh)=id++;
#endif // CGAL_PWRLE_TURN_V2

          cur_dh=m_map.template beta<2, 1>(cur_dh);
        }
        while(cur_dh!=dh1);
      }

      if (display_time)
      {
        t2.stop();
        std::cout<<"[TIME] Label darts: "<<t2.time()<<" seconds"<<std::endl;
      }
#endif // defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
    }

#ifdef CGAL_TRACE_CMAP_TOOLS
    std::cout<<"After quadrangulation: ";
    m_map.display_characteristics(std::cout) << ", valid=" 
                                             << m_map.is_valid() << std::endl;

    std::cout<<"Paths are all valid ? "<<(are_paths_valid()?"YES":"NO")
             <<std::endl;
    auto marktemp=m_map.get_new_mark();
    Dart_handle dh2=NULL;
    for (auto it=m_map.darts().begin(); it!=m_map.darts().end(); ++it)
    {
      if (!m_map.is_marked(it, marktemp))
      {
        std::cout<<"Degree="<<CGAL::template degree<Map, 0>(m_map, it)<<std::endl;
        std::cout<<"Co-degree="<<CGAL::template codegree<Map, 2>(m_map, it)<<std::endl;
        dh2=it;
        do
        {
          m_map.mark(dh2, marktemp);
          std::cout<<m_map.darts().index(dh2)<<"   "
                   <<m_map.darts().index(m_map.template beta<0>(dh2))<<std::endl;
          dh2=m_map.template beta<0,2>(dh2);
        }
        while(dh2!=it);
      }
    }
    m_map.free_mark(marktemp);
    m_map.display_darts(std::cout);
#endif

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] Total time for computation of reduced map: "
               <<t.time()<<" seconds"<<std::endl;
    }

    CGAL_assertion(m_map.darts().size()==4 || are_paths_valid()); // Because torus is a special case
  }
    
  ~Minimal_quadrangulation()
  {
    m_original_map.free_mark(m_mark_T);
    m_original_map.free_mark(m_mark_L);
  }

  /// @return true iff 'path' is contractible.
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

    if (m_map.is_empty())
    { return true; } // A closed path on a sphere is always contractible.

    CGAL::Timer t;
    if (display_time)
    { t.start(); }

    bool res=false;
    if (m_map.number_of_darts()==4)
    { // Case of torus
      Path_on_surface<CMap_for_minimal_quadrangulation>
        pt=transform_original_path_into_quad_surface_for_torus(p);

      int a, b;
      count_edges_of_path_on_torus(pt, a, b);
      res=(a==0 && b==0);
    }
    else
    {
      internal::Path_on_surface_with_rle<CMap_for_minimal_quadrangulation>
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

  /// @return true iff 'path1' and 'path2' are freely homotopic.
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
    if (m_map.is_empty())
    { return true; } // Two closed paths on a sphere are always homotopic.

    CGAL::Timer t;
    if (display_time)
    { t.start(); }

    bool res=false;
    if (m_map.number_of_darts()==4)
    { // Case of torus
      Path_on_surface<CMap_for_minimal_quadrangulation>
        pt1=transform_original_path_into_quad_surface_for_torus(p1);
      Path_on_surface<CMap_for_minimal_quadrangulation>
        pt2=transform_original_path_into_quad_surface_for_torus(p2);

      int a1, a2, b1, b2;
      count_edges_of_path_on_torus(pt1, a1, b1);
      count_edges_of_path_on_torus(pt2, a2, b2);
      res=(a1==a2 && b1==b2);
    }
    else
    {
      internal::Path_on_surface_with_rle<CMap_for_minimal_quadrangulation>
        pt1=transform_original_path_into_quad_surface_with_rle(p1);
      internal::Path_on_surface_with_rle<CMap_for_minimal_quadrangulation>
        pt2=transform_original_path_into_quad_surface_with_rle(p2);
      pt1.canonize();
      pt2.canonize();
      res=(pt1==pt2); // Do here to be counted in the computation time

#ifdef CGAL_TRACE_PATH_TESTS
      std::cout<<"Length of reduced paths: "<<pt1.length()<<" and "
               <<pt2.length()<<std::endl;
#endif

      /* std::cout<<"path1="<<Path_on_surface<CMap_for_minimal_quadrangulation>(path1)<<std::endl
         <<"path2="<<Path_on_surface<CMap_for_minimal_quadrangulation>(path2)<<std::endl;

         Path_on_surface<CMap_for_minimal_quadrangulation>(path1).display_pos_and_neg_turns(); std::cout<<std::endl;
         Path_on_surface<CMap_for_minimal_quadrangulation>(path2).display_pos_and_neg_turns(); std::cout<<std::endl;

         path1.display_pos_and_neg_turns(); std::cout<<std::endl;
         path2.display_pos_and_neg_turns(); std::cout<<std::endl; */
    }

    if (display_time)
    {
      t.stop();
      std::cout<<"[TIME] are_freely_homotopic: "<<t.time()<<" seconds"<<std::endl;
    }

    return res;
  }

  /// @return true iff 'path1' and 'path2' are base point freely homotopic.
  bool are_base_point_homotopic(const Path_on_surface<Mesh>& p1,
                                const Path_on_surface<Mesh>& p2,
                                bool display_time=false) const
  {
    if (p1.is_empty() && p2.is_empty()) { return true; }
    if (p1.is_empty() || p2.is_empty()) { return false; }

    if (!m_original_map.template belong_to_same_cell<0>(p1.front(), p2.front()) ||
        !m_original_map.template belong_to_same_cell<0>(m_original_map.other_extremity(p1.back()),
                                                        m_original_map.other_extremity(p2.back())))
    {
      std::cerr<<"Error: are_base_point_homotopic requires two paths that"
               <<" share the same vertices as extremities."<<std::endl;
      return false;
    }

    if (m_map.is_empty())
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

  const CMap_for_minimal_quadrangulation& get_map() const
  { return m_map; }

protected:
  void count_edges_of_path_on_torus
  (const Path_on_surface<CMap_for_minimal_quadrangulation>& path,
   int& a, int& b) const
  {
    CGAL_assertion(m_map.number_of_darts()==4);
      
    Dart_const_handle dha=m_map.darts().begin();
    Dart_const_handle dhb=m_map.template beta<1>(dha);

    a=0; b=0;
    for (std::size_t i=0; i<path.length(); ++i)
    {
      if (path[i]==dha) { ++a; }
      else if (path[i]==m_map.template beta<2>(dha)) { --a; }
      else if (path[i]==dhb) { ++b; }
      else if (path[i]==m_map.template beta<2>(dhb)) { --b; }
    }
  }

  Path_on_surface<CMap_for_minimal_quadrangulation>
  transform_original_path_into_quad_surface_for_torus(const Path_on_surface<Mesh>& path) const
  {
    CGAL_assertion(m_map.number_of_darts()==4);
      
    Path_on_surface<CMap_for_minimal_quadrangulation> res(m_map);
    if (path.is_empty()) return res;

    Dart_const_handle cur;
    for (std::size_t i=0; i<path.length(); ++i)
    {
      if (!m_original_map.is_marked(path[i], m_mark_T))
      {
        cur=get_first_dart_of_the_path(path[i], false);
        while(cur!=get_second_dart_of_the_path(path[i], false))
        {
          res.push_back(cur, false);
          cur=m_map.template beta<1>(cur);
        }
      }
    }
    res.update_is_closed();
    CGAL_assertion(res.is_empty() || res.is_closed());
    CGAL_assertion(res.is_valid());
    return res;
  }

  internal::Path_on_surface_with_rle<CMap_for_minimal_quadrangulation>
  transform_original_path_into_quad_surface_with_rle
  (const Path_on_surface<Mesh>& path) const
  {
    internal::Path_on_surface_with_rle<CMap_for_minimal_quadrangulation>
      res(m_map
#ifdef CGAL_PWRLE_TURN_V2
          , m_dart_ids
#endif //CGAL_PWRLE_TURN_V2
          );

    if (path.is_empty()) return res;

    for (std::size_t i=0; i<path.length(); ++i)
    {
      if (!m_original_map.is_marked(path[i], m_mark_T))
      {
        res.push_back(get_first_dart_of_the_path(path[i]), false);
        res.push_back(get_second_dart_of_the_path(path[i]), false);
      }
    }
    res.update_is_closed();
    res.merge_last_flat_with_next_if_possible();
    CGAL_assertion(res.is_closed());
    CGAL_assertion(res.is_valid());
    return res;
  }

  void initialize_vertices(UFTree2& uftrees,
                           boost::unordered_map
                           <typename Map::Dart_const_handle, UFTree_handle2>&
                           vertices)
  {
    uftrees.clear();
    vertices.clear();

    typename Map::size_type treated=m_original_map.get_new_mark();
    for (typename Map::Dart_range::const_iterator
           it=m_original_map.darts().begin(), itend=m_original_map.darts().end();
         it!=itend; ++it)
    {
      if (!m_original_map.is_marked(it, treated))
      {
        UFTree_handle2 newuf=uftrees.make_set(it);
        for (typename Map::
               template Dart_of_cell_range<0>::const_iterator
               itv=m_original_map.template darts_of_cell<0>(it).begin(),
               itvend=m_original_map.template darts_of_cell<0>(it).end();
             itv!=itvend; ++itv)
        {
          vertices[itv]=newuf;
          m_original_map.mark(itv, treated);
        }
      }
    }
    m_original_map.free_mark(treated);
  }

  void initialize_faces(UFTree& uftrees,
                        boost::unordered_map<Dart_const_handle, UFTree_handle>&
                        faces)
  {
    uftrees.clear();
    faces.clear();

    typename CMap_for_minimal_quadrangulation::size_type
      treated=m_map.get_new_mark();
    for (typename CMap_for_minimal_quadrangulation::Dart_range::iterator
           it=m_map.darts().begin(), itend=m_map.darts().end(); it!=itend;
         ++it)
    {
      if (!m_map.is_marked(it, treated))
      {
        UFTree_handle newuf=uftrees.make_set(it);
        Dart_handle cur=it;
        do
        {
          faces[cur]=newuf;
          m_map.mark(cur, treated);
          cur=m_map.template beta<1>(cur);
        }
        while (cur!=it);
      }
    }
    m_map.free_mark(treated);
  }

  UFTree_handle get_uftree(const UFTree& uftrees,
                           const boost::unordered_map<Dart_const_handle,
                           UFTree_handle>& mapdhtouf,
                           Dart_const_handle dh) const
  {
    CGAL_assertion(dh!=NULL);
    CGAL_assertion(mapdhtouf.find(dh)!=mapdhtouf.end());
    return uftrees.find(mapdhtouf.find(dh)->second);
  }

  UFTree_handle2 get_uftree2(const UFTree2& uftrees,
                             const boost::unordered_map<typename Map::Dart_const_handle,
                             UFTree_handle2>& mapdhtouf,
                             typename Map::Dart_const_handle dh) const
  {
    // CGAL_assertion(dh!=NULL);
    CGAL_assertion(mapdhtouf.find(dh)!=mapdhtouf.end());
    return uftrees.find(mapdhtouf.find(dh)->second);
  }

  /// Mark the edge containing adart in the original map.
  void mark_edge(const Map& amap, typename Map::Dart_const_handle adart,
                 std::size_t amark)
  {
    amap.mark(amap.template beta<2>(adart), amark);
    amap.mark(adart, amark);
  }

  /// Erase the edge given by adart (which belongs to the map m_map) from the
  /// associative array copy_to_origin, and erase the corresponding edge
  /// (which belongs to the map m_original_map) from the array origin_to_copy
  void erase_edge_from_associative_arrays
  (Dart_handle adart,
   boost::unordered_map<typename Map::Dart_const_handle, Dart_handle>&
   origin_to_copy,
   boost::unordered_map<Dart_handle, typename Map::Dart_const_handle>&
   copy_to_origin)
  {
    origin_to_copy.erase(m_original_map.template beta<2>
                         (copy_to_origin[adart]));
    origin_to_copy.erase(copy_to_origin[adart]);

    copy_to_origin.erase(m_map.template beta<2>(adart));
    copy_to_origin.erase(adart);
  }

  /// Step 1) Transform m_map into an equivalent surface having only one
  /// vertex. All edges contracted during this step belong to the spanning
  /// tree T, and thus corresponding edges in m_original_map are marked.
  void surface_simplification_in_one_vertex
  (boost::unordered_map<typename Map::Dart_const_handle, Dart_handle>&
   origin_to_copy,
   boost::unordered_map<Dart_handle, typename Map::Dart_const_handle>&
   copy_to_origin)
  {
    UFTree2 uftrees; // uftree of vertices; one tree for each vertex,
    // contains one dart of the vertex
    boost::unordered_map<typename Map::Dart_const_handle, UFTree_handle2> vertices;
    initialize_vertices(uftrees, vertices);

    /*    m_map.set_automatic_attributes_management(false);

          for (typename CMap_for_minimal_quadrangulation::Dart_range::iterator
          it=m_map.darts().begin(), itend=m_map.darts().end();
          it!=itend; ++it)
          {
          if (m_map.is_dart_used(it) &&
          get_uftree(uftrees, vertices, it)!=
          get_uftree(uftrees, vertices, m_map.template beta<2>(it)))
          {
          mark_edge(m_original_map, copy_to_origin[it], m_mark_T);
          erase_edge_from_associative_arrays(it, origin_to_copy, copy_to_origin);

          uftrees.unify_sets(get_uftree(uftrees, vertices, it),
          get_uftree(uftrees, vertices,
          m_map.template beta<2>(it)));
          //m_map.template contract_cell<1>(it);
          Dart_handle d1=it, d2=m_map.template beta<2>(it);
          m_map.template link_beta<1>(m_map.template beta<0>(d1),
          m_map.template beta<1>(d1));
          m_map.template link_beta<1>(m_map.template beta<0>(d2),
          m_map.template beta<1>(d2));
          m_map.erase_dart(d1);
          m_map.erase_dart(d2);
          }
          }

          m_map.set_automatic_attributes_management(true); */

    /* New version that does not need to first copy the map before to simplify it
     */
    Dart_handle d1, d2;
    for (typename Map::Dart_range::const_iterator
           it=m_original_map.darts().begin(), itend=m_original_map.darts().end();
         it!=itend; ++it)
    {
      if (typename Map::Dart_const_handle(it)<m_original_map.template beta<2>(it))
      {
        if (get_uftree2(uftrees, vertices, it)!=
            get_uftree2(uftrees, vertices, m_original_map.template beta<2>(it)))
        {
          m_original_map.mark(m_original_map.template beta<2>(it), m_mark_T);
          m_original_map.mark(it, m_mark_T);

          uftrees.unify_sets(get_uftree2(uftrees, vertices, it),
                             get_uftree2(uftrees, vertices,
                                         m_original_map.template beta<2>(it)));
        }
        else
        {
          d1=m_map.create_dart();
          d2=m_map.create_dart();
          m_map.template basic_link_beta_for_involution<2>(d1, d2);
          origin_to_copy[it]=d1;
          origin_to_copy[m_original_map.template beta<2>(it)]=d2;
          copy_to_origin[d1]=it;
          copy_to_origin[d2]=m_original_map.template beta<2>(it);
        }
      }
    }

    /// Now we only need to do the basic_link_beta_1
    typename Map::Dart_const_handle dd1;
    for (typename Map::Dart_range::const_iterator
           it=m_original_map.darts().begin(), itend=m_original_map.darts().end();
         it!=itend; ++it)
    {
      if (!m_original_map.is_marked(it, m_mark_T))
      {
        dd1=m_original_map.template beta<1>(it);
        while(m_original_map.is_marked(dd1, m_mark_T))
        { dd1=m_original_map.template beta<1>(dd1); }
        m_map.basic_link_beta_1(origin_to_copy[it], origin_to_copy[dd1]);
      }
    }
  }

  /// Step 2) Compute, for each edge of m_original_map not in the spanning
  /// tree T, the pair of darts of the edge in m_copy. This pair of edges
  /// will be updated later (in surface_simplification_in_one_face() and in
  /// surface_quadrangulate() )
  void compute_length_two_paths
  (const boost::unordered_map<typename Map::Dart_const_handle, Dart_handle>&
   origin_to_copy)
  {
    paths.clear();

    for (typename Map::Dart_range::const_iterator
           it=m_original_map.darts().begin(),
           itend=m_original_map.darts().end(); it!=itend; ++it)
    {
      if (!m_original_map.is_marked(it, m_mark_T))
      {
        CGAL_assertion(!m_original_map.template is_free<2>(it));
        if (typename Map::Dart_const_handle(it)<m_original_map.template beta<2>(it))
        {
          paths[it]=std::make_pair
            (origin_to_copy.at(it),
             m_map.template beta<2>(origin_to_copy.at(it)));
          CGAL_assertion(paths[it].first!=paths[it].second);
          CGAL_assertion(paths[it].first==m_map.template beta<2>(paths[it].second));
        }
      }
    }

#ifdef CGAL_TRACE_CMAP_TOOLS
    std::cout<<"Number of darts in paths: "<<paths.size()
             <<"; number of darts in m_map: "<<m_map.number_of_darts()
             <<std::endl;
#endif
  }

  /// Step 3) Transform the 2-map into an equivalent surface having only
  /// one face. All edges removed during this step belong to the
  /// dual spanning tree L (spanning tree of the dual 2-map).
  void surface_simplification_in_one_face
  (boost::unordered_map<typename Map::Dart_const_handle, Dart_handle>&
   origin_to_copy,
   boost::unordered_map<Dart_handle, typename Map::Dart_const_handle>&
   copy_to_origin)
  {
    UFTree uftrees; // uftree of faces; one tree for each face,
    // contains one dart of the face
    boost::unordered_map<Dart_const_handle, UFTree_handle> faces;
    initialize_faces(uftrees, faces);

    m_map.set_automatic_attributes_management(false);

    typename Map::size_type toremove=m_map.get_new_mark();
    Dart_handle currentdart=NULL, oppositedart=NULL;
      
    for (typename CMap_for_minimal_quadrangulation::Dart_range::iterator
           it=m_map.darts().begin(), itend=m_map.darts().end(); it!=itend;
         ++it)
    {
      currentdart=it;
      CGAL_assertion(!m_map.template is_free<2>(currentdart)); // TODO later, support opened surfaces
      oppositedart=m_map.template beta<2>(currentdart);

      if (currentdart<oppositedart && !m_map.is_marked(currentdart, toremove))
      {
        // We remove degree two edges (we cannot have dangling edges
        // because we had previously contracted all the non loop and thus
        // we have only one vertex).
        if (get_uftree(uftrees, faces, currentdart)!=
            get_uftree(uftrees, faces, oppositedart))
        {
          // We cannot have a dangling edge
          CGAL_assertion(m_map.template beta<0>(currentdart)!=oppositedart);
          CGAL_assertion(m_map.template beta<1>(currentdart)!=oppositedart);

          uftrees.unify_sets(get_uftree(uftrees, faces, currentdart),
                             get_uftree(uftrees, faces, oppositedart));

          m_map.mark(currentdart, toremove);
          m_map.mark(oppositedart, toremove);
          mark_edge(m_original_map, copy_to_origin[currentdart], m_mark_L);
        }
      }
    }

    /* m_map.display_characteristics(std::cout) << ", valid="
       << m_map.is_valid() << std::endl;
       m_map.display_darts(std::cout)<<std::endl; */

    if (m_map.number_of_marked_darts(toremove)==m_map.number_of_darts())
    {
      // Case of sphere; all darts are removed.
      paths.clear();
      m_map.clear();
    }
    else
    {
      update_length_two_paths_before_edge_removals(toremove, copy_to_origin);

      // We remove all the edges to remove.
      for (typename CMap_for_minimal_quadrangulation::Dart_range::iterator
             it=m_map.darts().begin(), itend=m_map.darts().end(); it!=itend;
           ++it)
      {
        if (m_map.is_dart_used(it) && m_map.is_marked(it, toremove))
        {
          erase_edge_from_associative_arrays(it, origin_to_copy, copy_to_origin);
          // TODO LATER (?) OPTIMIZE AND REPLACE THE REMOVE_CELL CALL BY THE MODIFICATION BY HAND
          // OR DEVELOP A SPECIALIZED VERSION OF REMOVE_CELL
          m_map.template remove_cell<1>(it);
        }
      }
    }

    m_map.set_automatic_attributes_management(true);
    m_map.free_mark(toremove);
  }

  /// Step 4) quadrangulate the surface.
  void surface_quadrangulate()
  {
    // Here the map has only one face and one vertex.
    typename Map::size_type oldedges=m_map.get_new_mark();
    m_map.negate_mark(oldedges); // now all edges are marked
      
    // 1) We insert a vertex in the face.
    //    New edges created by the operation are not marked.
    m_map.insert_cell_0_in_cell_2(m_map.darts().begin());

    // m_map.display_darts(std::cout);

    // 2) We update the pair of darts
    // std::cout<<"************************************************"<<std::endl;
    for (typename TPaths::iterator itp=paths.begin(), itpend=paths.end();
         itp!=itpend; ++itp)
    {
      std::pair<Dart_const_handle, Dart_const_handle>& p=itp->second;
      //std::cout<<"Pair: "<<m_map.darts().index(p.first)<<", "
      //         <<m_map.darts().index(p.second)<<std::flush;
      p.first=m_map.template beta<0, 2>(p.first);
      p.second=m_map.template beta<0>(p.second);
      //std::cout<<" -> "<<m_map.darts().index(p.first)<<", "
      //         <<m_map.darts().index(p.second)<<std::endl;
    }

    // 3) We remove all the old edges.
    for (typename CMap_for_minimal_quadrangulation::Dart_range::iterator
           it=m_map.darts().begin(), itend=m_map.darts().end(); it!=itend;
         ++it)
    {
      if (m_map.is_dart_used(it) && m_map.is_marked(it, oldedges))
      { m_map.template remove_cell<1>(it); }
    }

    m_map.free_mark(oldedges);
  }

  /// Update all length two paths, before edge removal. Edges that will be
  /// removed are marked with toremove mark.
  void update_length_two_paths_before_edge_removals
  (typename Map::size_type toremove,
   const boost::unordered_map<Dart_handle,
   typename Map::Dart_const_handle>& copy_to_origin)
  {
    // std::cout<<"************************************************"<<std::endl;
    for (auto it=m_original_map.darts().begin();
         it!=m_original_map.darts().end(); ++it)
    {
      if (!m_original_map.is_marked(it, m_mark_T) &&
          !m_original_map.is_marked(it, m_mark_L) &&
          typename Map::Dart_const_handle(it)<m_original_map.template beta<2>(it))
      { // Surviving dart => belongs to the border of the face
        std::pair<Dart_const_handle, Dart_const_handle>& p=paths[it];

        Dart_handle initdart=m_map.darts().iterator_to
          (const_cast<typename CMap_for_minimal_quadrangulation::Dart &>
           (*(p.first)));
        Dart_handle initdart2=m_map.template beta<2>(initdart);
        CGAL_assertion(initdart2==p.second);
        CGAL_assertion(!m_map.is_marked(initdart, toremove));
        CGAL_assertion(!m_map.is_marked(initdart2, toremove));

        // 1) We update the dart associated with p.second
        p.second=m_map.template beta<1>(initdart);
        while (m_map.is_marked(p.second, toremove))
        { p.second=m_map.template beta<2, 1>(p.second); }

        // 2) We do the same loop, linking all the inner darts with p.second
        initdart=m_map.template beta<1>(initdart);
        while (m_map.is_marked(initdart, toremove))
        {
          CGAL_assertion(copy_to_origin.count(initdart)==1);
          typename Map::Dart_const_handle
            d1=copy_to_origin.find(initdart)->second;
          typename Map::Dart_const_handle
            d2=m_original_map.template beta<2>(d1);
          if (d1<d2) { paths[d1].first=p.second; }
          else       { paths[d2].second=p.second; }
          initdart=m_map.template beta<2, 1>(initdart);
        }

        // 3) We do the same loop but starting from initdart2
        initdart2=m_map.template beta<1>(initdart2);
        Dart_handle enddart2=initdart2;
        while (m_map.is_marked(enddart2, toremove))
        { enddart2=m_map.template beta<2, 1>(enddart2); }

        while (m_map.is_marked(initdart2, toremove))
        {
          CGAL_assertion(copy_to_origin.count(initdart2)==1);
          typename Map::Dart_const_handle
            d1=copy_to_origin.find(initdart2)->second;
          typename Map::Dart_const_handle
            d2=m_original_map.template beta<2>(d1);
          if (d1<d2) {
            CGAL_assertion(paths.count(d1)==1);
            paths[d1].first=enddart2;
          }
          else       {
            CGAL_assertion(paths.count(d2)==1);
            paths[d2].second=enddart2;
          }
          initdart2=m_map.template beta<2, 1>(initdart2);
        }
      }
    }
  }

  /// @return true iff the edge containing adart is associated with a path.
  ///         (used for debug purpose because we are suppose to be able to
  ///          test this by using directly the mark m_mark_T).
  bool is_edge_has_path(typename Map::Dart_const_handle adart) const
  {
    typename Map::Dart_const_handle
      opposite=m_original_map.template beta<2>(adart);
    if (adart<opposite)
    {
      return paths.find(adart)!=paths.end();
    }
    return paths.find(opposite)!=paths.end();
  }

  /// @return the pair of darts associated with the edge containing adart
  ///         in m_original_map.
  /// @pre the edge containing adart must not belong to T.
  std::pair<Dart_const_handle, Dart_const_handle>& get_pair_of_darts
  (typename Map::Dart_const_handle adart)
  {
    CGAL_assertion(!m_original_map.is_marked(adart, m_mark_T));
    CGAL_assertion(is_edge_has_path(adart));

    typename Map::Dart_const_handle
      opposite=m_original_map.template beta<2>(adart);
    if (adart<opposite)
    { return paths.find(adart)->second; }

    return paths.find(opposite)->second;
  }

  Dart_const_handle get_first_dart_of_the_path
  (typename Map::Dart_const_handle adart, bool withopposite=true) const
  {
    CGAL_assertion(!m_original_map.is_marked(adart, m_mark_T));
    CGAL_assertion(is_edge_has_path(adart));

    typename Map::Dart_const_handle
      opposite=m_original_map.template beta<2>(adart);
    if (adart<opposite)
    {
      const std::pair<Dart_const_handle, Dart_const_handle>&
        p=paths.find(adart)->second;
      return p.first;
    }

    const std::pair<Dart_const_handle, Dart_const_handle>&
      p=paths.find(opposite)->second;
    return (withopposite?m_map.template beta<2>(p.second):p.second);
  }

  Dart_const_handle get_second_dart_of_the_path
  (typename Map::Dart_const_handle adart, bool withopposite=true) const
  {
    CGAL_assertion(!m_original_map.is_marked(adart, m_mark_T));
    CGAL_assertion(is_edge_has_path(adart));

    typename Map::Dart_const_handle
      opposite=m_original_map.template beta<2>(adart);
    if (adart<opposite)
    {
      const std::pair<Dart_const_handle, Dart_const_handle>&
        p=paths.find(adart)->second;
      return p.second;
    }

    const std::pair<Dart_const_handle, Dart_const_handle>&
      p=paths.find(opposite)->second;
    return (withopposite?m_map.template beta<2>(p.first):p.first);
  }

  /// Test if paths are valid, i.e.:
  /// 1) all the darts of m_original_map that do not belong to T are
  ///    associated with a pair of darts;
  /// 2) all the darts of the paths belong to m_map;
  /// 3) the origin of the second dart of the pair is the extremity of the
  ///    first dart.
  /// 4) all the darts of m_map are not free (both for beta 1 and 2)
  /// 5) The two darts in a pair are different
  bool are_paths_valid() const
  {
    if (paths.empty()) { return true; }

    bool res=true;
    for (auto it=m_original_map.darts().begin(),
           itend=m_original_map.darts().end(); it!=itend; ++it)
    {
      if (!m_original_map.is_marked(it, m_mark_T))
      {
        if (!is_edge_has_path(it))
        {
          std::cout<<"ERROR: an edge that does not belong to the spanning "
                   <<"tree T has no associated path."<<std::endl;
          res=false;
        }
      }
      else
      {
        if (is_edge_has_path(it))
        {
          std::cout<<"ERROR: an edge that belongs to the spanning tree"
                   <<" T has an associated path."<<std::endl;
          res=false;
        }
      }
    }

    for (auto it=m_map.darts().begin(),
           itend=m_map.darts().end(); it!=itend; ++it)
    {
      if (m_map.is_free(it, 1))
      {
        std::cout<<"ERROR: a dart of the quandrangulated map is 1-free"
                 <<std::endl;
        res=false;
      }
      if (m_map.is_free(it, 2))
      {
        std::cout<<"ERROR: a dart of the quandrangulated map is 2-free"
                 <<std::endl;
        res=false;
      }
    }

    for (auto it=paths.begin(); it!=paths.end(); ++it)
    {
      if (!m_map.is_dart_used(it->second.first))
      {
        std::cout<<"ERROR: first dart in paths does not exist anymore in m_map."
                 <<std::endl;
        res=false;
      }
      else if (!m_map.darts().owns(it->second.first))
      {
        std::cout<<"ERROR: first dart in paths does not belong to m_map."
                 <<std::endl;
        res=false;
      }
      if (!m_map.is_dart_used(it->second.second))
      {
        std::cout<<"ERROR: second dart in paths does not exist anymore in m_map."
                 <<std::endl;
        res=false;
      }
      else if (!m_map.darts().owns(it->second.second))
      {
        std::cout<<"ERROR: second dart in paths does not belong to m_map."
                 <<std::endl;
        res=false;
      }
      if (it->second.first==it->second.second)
      {
        std::cout<<"ERROR: two darts in the same pair are equal."
                 <<std::endl;
        res=false;
      }
    }

    for (auto it=m_original_map.darts().begin(),
           itend=m_original_map.darts().end(); it!=itend; ++it)
    {
      if (!m_original_map.is_marked(it, m_mark_T))
      {
        Dart_const_handle d1=get_first_dart_of_the_path(it);
        Dart_const_handle d2=get_second_dart_of_the_path(it);
        if (d1==NULL || d2==NULL)
        {
          std::cout<<"ERROR: an edge is associated with a null dart in paths."
                   <<std::endl;
          res=false;
        }
        else
        {
          Dart_const_handle dd1=m_map.other_extremity(d1);
          CGAL_assertion(dd1!=NULL);
          if (!CGAL::belong_to_same_cell<CMap_for_minimal_quadrangulation,0>(m_map, dd1, d2))
          {
            std::cout<<"ERROR: the two darts in a path are not consecutive."
                     <<std::endl;
            res=false;
          }
        }
      }
    }

    return res;
  }

protected:
  const typename Get_map<Mesh, Mesh>::storage_type m_original_map; // The original map
  CMap_for_minimal_quadrangulation m_map; /// the transformed map
  TPaths paths; /// Pair of edges associated with each edge of m_original_map
  /// (except the edges that belong to the spanning tree T).
  std::size_t m_mark_T; /// mark each edge of m_original_map that belong to the spanning tree T
  std::size_t m_mark_L; /// mark each edge of m_original_map that belong to the dual spanning tree L

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
