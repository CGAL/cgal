// Copyright (c) 2017 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_SURFACE_MESH_CURVE_TOPOLOGY_H
#define CGAL_SURFACE_MESH_CURVE_TOPOLOGY_H 1

#include <CGAL/Union_find.h>
#include <CGAL/Random.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Path_on_surface_with_rle.h>
#include <CGAL/Path_generators.h>
#include <CGAL/Combinatorial_map_basic_operations.h>
#include <CGAL/Timer.h>
#include <boost/unordered_map.hpp>
#include <stack>

namespace CGAL {
  
  template<typename Map>
  class Surface_mesh_curve_topology
  {
  public:
    typedef typename Map::Dart_handle Dart_handle;
    typedef typename Map::Dart_const_handle Dart_const_handle;
    typedef CGAL::Union_find<Dart_handle> UFTree;
    typedef typename UFTree::handle UFTree_handle;
    
    typedef boost::unordered_map<Dart_const_handle,
                      std::pair<Dart_const_handle, Dart_const_handle> > TPaths;
    typedef boost::unordered_map<Dart_const_handle, std::size_t> TDartIds;

    Surface_mesh_curve_topology(Map& amap, bool display_time=false) :
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
      boost::unordered_map<Dart_const_handle, Dart_handle> origin_to_copy;

      // We copy the original map, while keeping a mapping between darts.
      m_map.copy(m_original_map, &origin_to_copy);

      // The mapping between darts of the copy into darts of the original map.
      boost::unordered_map<Dart_handle, Dart_const_handle> copy_to_origin;
      for (auto it=origin_to_copy.begin(); it!=origin_to_copy.end(); ++it)
      { copy_to_origin[it->second]=it->first; }

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

     /* std::cout<<"Paths are all valid 2 ? "<<(are_paths_valid()?"YES":"NO")
               <<std::endl;*/

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
        // of turns in constant time.
        CGAL_assertion(m_map.number_of_darts()%2==0);
        m_number_of_edges=m_map.number_of_darts()/2;

        Dart_handle dh1=m_map.darts().begin();
        Dart_handle dh2=m_map.template beta<2>(dh1);
        std::size_t id=0;
        for(; dh1!=NULL; dh1=(dh1==dh2?NULL:dh2)) // We have two vertices to process
        {
          Dart_handle cur_dh=dh1;
          do
          {
            m_dart_ids[cur_dh]=id++;
            cur_dh=m_map.template beta<2, 1>(cur_dh);
          }
          while(cur_dh!=dh1);
        }

        if (display_time)
        {
          t2.stop();
          std::cout<<"[TIME] Label darts: "<<t2.time()<<" seconds"<<std::endl;
        }
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

      assert(m_map.darts().size()==4 || are_paths_valid()); // Because torus is a special case
    }
    
    ~Surface_mesh_curve_topology()
    {
      m_original_map.free_mark(m_mark_T);
      m_original_map.free_mark(m_mark_L);
    }

    const Map& get_map() const
    { return m_map; }
    
    /// Canonize the path
    void canonize(Path_on_surface<Map>& path, bool display_time=false) const
    {
      if (!path.is_closed())
      {
        std::cerr<<"Error: a non closed path cannot be canonized."<<std::endl;
        return;
      }

      CGAL::Timer t;
      if (display_time)
      { t.start(); }

      Path_on_surface_with_rle<Map> path2(path);
      path2.canonize();
      path=path2;

      if (display_time)
      {
        t.stop();
        std::cout<<"[TIME] Canonize path: "<<t.time()<<" seconds"<<std::endl;
      }

      CGAL_assertion(path.is_valid());
    }

    void count_edges_of_path_on_torus(const Path_on_surface<Map>& path,
                                      std::size_t& a,
                                      std::size_t& b)
    {
      Dart_const_handle dha=m_map.darts().begin();
      if (dha>m_map.template beta<2>(dha)) { dha=m_map.template beta<2>(dha); }

      Dart_const_handle dhb=NULL;
      auto it=m_map.darts().begin();
      while(dhb==NULL)
      {
        if (it!=dha && it!=m_map.template beta<2>(dha)) { dhb=it; }
        else { ++it; }
      }

      if (dhb>m_map.template beta<2>(dhb)) { dhb=m_map.template beta<2>(dhb); }

      if (dha<dhb) { std::swap(dha, dhb); }

      a=0; b=0;

      for (int i=0; i<path.length(); ++i)
      {
        if (path[i]==dha) { ++a; }
        else if (path[i]==m_map.template beta<2>(dha)) { --a; }
        else if (path[i]==dhb) { ++b; }
        else if (path[i]==m_map.template beta<2>(dhb)) { --b; }
      }
    }

    /// @return true iff 'path' is contractible.
    bool is_contractible(const Path_on_surface<Map>& p,
                         bool display_time=false)
    {
      CGAL::Timer t;
      if (display_time)
      { t.start(); }

      Path_on_surface<Map> pt=transform_original_path_into_quad_surface(p);

      bool res=false;
      if (m_map.number_of_darts()==4)
      { // Case of torus
        std::size_t a, b;
        count_edges_of_path_on_torus(pt, a, b);
        res=(a==0 && b==0);
      }
      else
      {
        Path_on_surface_with_rle<Map> path(pt);
        path.canonize();
        res=path.is_empty();
        std::cout<<"Length of reduced paths: "<<path.length()<<std::endl;
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
    bool are_freely_homotopic(const Path_on_surface<Map>& p1,
                              const Path_on_surface<Map>& p2,
                              bool display_time=false)
    {
      CGAL::Timer t;
      if (display_time)
      { t.start(); }

      Path_on_surface<Map> pt1=transform_original_path_into_quad_surface(p1);
      Path_on_surface<Map> pt2=transform_original_path_into_quad_surface(p2);

      bool res=false;
      if (m_map.number_of_darts()==4)
      { // Case of torus
        std::size_t a1, a2, b1, b2;
        count_edges_of_path_on_torus(p1, a1, b1);
        count_edges_of_path_on_torus(p2, a2, b2);
        res=(a1==a2 && b1==b2);
      }
      else
      {
        Path_on_surface_with_rle<Map> path1(pt1);
        Path_on_surface_with_rle<Map> path2(pt2);

        path1.canonize();
        path2.canonize();
        res=(path1==path2); // Do here to be counted in the computation time

        std::cout<<"Length of reduced paths: "<<path1.length()<<" and "
                 <<path2.length()<<std::endl;
      }

      if (display_time)
      {
        t.stop();
        std::cout<<"[TIME] are_freely_homotopic: "<<t.time()<<" seconds"<<std::endl;
      }

      return res;
    }

    /// @return true iff 'path1' and 'path2' are base point freely homotopic.
    bool are_base_point_homotopic(const Path_on_surface<Map>& p1,
                                  const Path_on_surface<Map>& p2,
                                  bool display_time=false) const
    {
      if (p1.is_empty() && p2.is_empty()) { return true; }
      if (p1.is_empty() || p2.is_empty()) { return true; }

      if (!CGAL::belong_to_same_cell<Map,0>(m_map,
                                            p1.front(), p2.front()))
      {
        std::cerr<<"Error: are_base_point_homotopic requires two paths that"
                 <<" share the same first vertex."<<std::endl;
        return false;
      }

      CGAL::Timer t;
      if (display_time)
      { t.start(); }

      Path_on_surface<Map> path=p1;
      Path_on_surface<Map> path2=p2; p2.reverse();
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

  protected:
    Path_on_surface<Map> transform_original_path_into_quad_surface
    (const Path_on_surface<Map>& path) // TODO return directly a Path_on_durface_with_rle
    {
      Path_on_surface<Map> res(m_map);
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
      assert(res.is_closed());
      assert(res.is_valid());
      return res;
    }

    void initialize_vertices(UFTree& uftrees,
                             boost::unordered_map<Dart_const_handle, UFTree_handle>&
                             vertices)
    {
      uftrees.clear();
      vertices.clear();

      typename Map::size_type treated=m_map.get_new_mark();
      for (typename Map::Dart_range::iterator it=m_map.darts().begin(),
           itend=m_map.darts().end(); it!=itend; ++it)
      {
        if (!m_map.is_marked(it, treated))
        {
          UFTree_handle newuf=uftrees.make_set(it);
          for (typename Map::template Dart_of_cell_basic_range<0>::iterator
               itv=m_map.template darts_of_cell_basic<0>(it, treated).begin(),
               itvend=m_map.template darts_of_cell_basic<0>(it, treated).end();
               itv!=itvend; ++itv)
          {
            vertices[itv]=newuf;
            m_map.mark(itv, treated);
          }
        }
      }
      m_map.free_mark(treated);
    }

    void initialize_faces(UFTree& uftrees,
                          boost::unordered_map<Dart_const_handle, UFTree_handle>&
                          faces)
    {
      uftrees.clear();
      faces.clear();

      typename Map::size_type treated=m_map.get_new_mark();
      for (typename Map::Dart_range::iterator it=m_map.darts().begin(),
           itend=m_map.darts().end();
           it!=itend; ++it)
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
                             Dart_const_handle dh)
    {
      assert(dh!=NULL);
      assert(mapdhtouf.find(dh)!=mapdhtouf.end());
      return uftrees.find(mapdhtouf.find(dh)->second);
    }

    /// Mark the edge containing adart in the given map.
    void mark_edge(const Map& amap, Dart_const_handle adart, std::size_t amark)
    {
      amap.mark(amap.template beta<2>(adart), amark);
      amap.mark(adart, amark);
    }

    /// Erase the edge given by adart (which belongs to the map m_map) from the
    /// associative array copy_to_origin, and erase the corresponding edge
    /// (which belongs to the map m_original_map) from the array origin_to_copy
    void erase_edge_from_associative_arrays
    (Dart_handle adart,
     boost::unordered_map<Dart_const_handle, Dart_handle>& origin_to_copy,
     boost::unordered_map<Dart_handle, Dart_const_handle>& copy_to_origin)
    {
      origin_to_copy.erase(m_original_map.template beta<2>(copy_to_origin[adart]));
      origin_to_copy.erase(copy_to_origin[adart]);

      copy_to_origin.erase(m_map.template beta<2>(adart));
      copy_to_origin.erase(adart);
    }

    /// Step 1) Transform m_map into an equivalent surface having only one
    /// vertex. All edges contracted during this step belong to the spanning
    /// tree T, and thus corresponding edges in m_original_map are marked.
    void surface_simplification_in_one_vertex
    (boost::unordered_map<Dart_const_handle, Dart_handle>& origin_to_copy,
     boost::unordered_map<Dart_handle, Dart_const_handle>& copy_to_origin)
    {
      UFTree uftrees; // uftree of vertices; one tree for each vertex,
                      // contains one dart of the vertex
      boost::unordered_map<Dart_const_handle, UFTree_handle> vertices;
      initialize_vertices(uftrees, vertices);

      m_map.set_automatic_attributes_management(false);

      for (typename Map::Dart_range::iterator it=m_map.darts().begin(),
           itend=m_map.darts().end(); it!=itend; ++it)
      {
        if (m_map.is_dart_used(it) &&
            get_uftree(uftrees, vertices, it)!=
            get_uftree(uftrees, vertices, m_map.template beta<2>(it)))
        {
          mark_edge(m_original_map, copy_to_origin[it], m_mark_T);
          erase_edge_from_associative_arrays(it, origin_to_copy, copy_to_origin);

          uftrees.unify_sets(get_uftree(uftrees, vertices, it),
                             get_uftree(uftrees, vertices, m_map.template beta<2>(it)));
          //m_map.template contract_cell<1>(it);
          Dart_handle d1=it, d2=m_map.template beta<2>(it);
          m_map.template link_beta<1>(m_map.template beta<0>(d1), m_map.template beta<1>(d1));
          m_map.template link_beta<1>(m_map.template beta<0>(d2), m_map.template beta<1>(d2));
          m_map.erase_dart(d1);
          m_map.erase_dart(d2);
        }
      }

      m_map.set_automatic_attributes_management(true);
    }

    /// Step 2) Compute, for each edge of m_original_map not in the spanning
    /// tree T, the pair of darts of the edge in m_copy. This pair of edges
    /// will be updated later (in surface_simplification_in_one_face() and in
    /// surface_quadrangulate() )
    void compute_length_two_paths
    (const boost::unordered_map<Dart_const_handle, Dart_handle>& origin_to_copy)
    {
      paths.clear();

      for (typename Map::Dart_range::const_iterator
           it=m_original_map.darts().begin(),
           itend=m_original_map.darts().end(); it!=itend; ++it)
      {
        if (!m_original_map.is_marked(it, m_mark_T))
        {
          assert(!m_original_map.template is_free<2>(it));
          if (it<m_original_map.template beta<2>(it))
          {
            paths[it]=std::make_pair
                (origin_to_copy.at(it),
                 m_map.template beta<2>(origin_to_copy.at(it)));
            assert(paths[it].first!=paths[it].second);
            assert(paths[it].first==m_map.template beta<2>(paths[it].second));
          }
        }
      }

#ifdef CGAL_TRACE_CMAP_TOOLS
      std::cout<<"Number of darts in paths: "<<paths.size()
               <<"; number of darts in m_map: "<<m_map.number_of_darts()<<std::endl;
#endif
    }

    /// Step 3) Transform the 2-map into an equivalent surface having only
    /// one face. All edges removed during this step belong to the
    /// dual spanning tree L (spanning tree of the dual 2-map).
    void surface_simplification_in_one_face
    (boost::unordered_map<Dart_const_handle, Dart_handle>& origin_to_copy,
     boost::unordered_map<Dart_handle, Dart_const_handle>& copy_to_origin)
    {
      UFTree uftrees; // uftree of faces; one tree for each face,
                      // contains one dart of the face
      boost::unordered_map<Dart_const_handle, UFTree_handle> faces;
      initialize_faces(uftrees, faces);

      m_map.set_automatic_attributes_management(false);

      typename Map::size_type toremove=m_map.get_new_mark();
      Dart_handle currentdart=NULL, oppositedart=NULL;
      
      for (typename Map::Dart_range::iterator it=m_map.darts().begin(),
           itend=m_map.darts().end(); it!=itend; ++it)
      {
        currentdart=it;
        assert (!m_map.template is_free<2>(currentdart)); // TODO later
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
            assert(m_map.template beta<0>(currentdart)!=oppositedart);
            assert(m_map.template beta<1>(currentdart)!=oppositedart);

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
        // update_length_two_paths_before_edge_removals_v1(toremove);
        update_length_two_paths_before_edge_removals_v2(toremove, copy_to_origin);

        // We remove all the edges to remove.
        for (typename Map::Dart_range::iterator it=m_map.darts().begin(),
             itend=m_map.darts().end();
             it!=itend; ++it)
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
      
      // 1) We insert a vertex in the face (note that all points have the
      //    same geometry). New edges created by the operation are not marked.
      m_map.insert_point_in_cell_2(m_map.darts().begin(),
                                   m_map.point(m_map.darts().begin()));

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
        // WRONG ASSERTS assert(p.first!=p.second);
        // assert(p.first!=m_map.template beta<2>(p.second));
      }

      // 3) We remove all the old edges.
      for (typename Map::Dart_range::iterator it=m_map.darts().begin(),
           itend=m_map.darts().end();
           it!=itend; ++it)
      {
        if (m_map.is_dart_used(it) && m_map.is_marked(it, oldedges))
        { m_map.template remove_cell<1>(it); }
      }

      m_map.free_mark(oldedges);
    }

    /// Version1: quadratic in number of darts in the paths
    void update_length_two_paths_before_edge_removals_v1(typename Map::size_type toremove)
    {
      // std::cout<<"************************************************"<<std::endl;
      // We update the pair of darts
      for (typename TPaths::iterator itp=paths.begin(), itpend=paths.end();
           itp!=itpend; ++itp)
      {
        std::pair<Dart_const_handle, Dart_const_handle>& p=itp->second;
        //std::cout<<"Pair: "<<m_map.darts().index(p.first)<<", "
        //         <<m_map.darts().index(p.second)<<": "<<std::flush;

        //std::cout<<m_map.darts().index(p.first)<<"; "<<std::flush;
        // p.first=m_map.template beta<0>(p.first);
        Dart_const_handle initdart=p.first;

        if (!m_map.is_marked(p.first, toremove))
        {
          p.second=m_map.template beta<1>(p.first);
          initdart=p.second;
          while (m_map.is_marked(p.second, toremove))
          {
            p.second=m_map.template beta<2, 1>(p.second);
            assert(p.second!=initdart);
          }
        }
        else
        {
          while (m_map.is_marked(p.first, toremove))
          {
            p.first=m_map.template beta<2, 1>(p.first);
            //std::cout<<m_map.darts().index(p.first)<<"; "<<std::flush;
            assert(p.first!=initdart);
          }
          //std::cout<<std::endl;
          // p.second=m_map.template beta<0>(p.second);
          initdart=p.second;
          while (m_map.is_marked(p.second, toremove))
          {
            p.second=m_map.template beta<2, 1>(p.second);
            //std::cout<<m_map.darts().index(p.second)<<"; "<<std::flush;
            assert(p.second!=initdart);
          }
          //std::cout<<std::endl;
          //std::cout<<" -> "<<m_map.darts().index(p.first)<<", "
          //         <<m_map.darts().index(p.second)<<std::endl;
        }
      }
    }

    /// Version2: linear in number of darts in the paths
    void update_length_two_paths_before_edge_removals_v2(typename Map::size_type toremove,
                                                         const boost::unordered_map<Dart_handle, Dart_const_handle>& copy_to_origin)
    {
      // std::cout<<"************************************************"<<std::endl;
      for (auto it=m_original_map.darts().begin(); it!=m_original_map.darts().end(); ++it)
      {
        if (!m_original_map.is_marked(it, m_mark_T) &&
            !m_original_map.is_marked(it, m_mark_L) &&
            it<m_original_map.template beta<2>(it))
        { // Surviving dart => belongs to the border of the face
          std::pair<Dart_const_handle, Dart_const_handle>& p=paths[it];

          Dart_handle initdart=m_map.darts().iterator_to(const_cast<typename Map::Dart &>(*(p.first)));
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
            assert(copy_to_origin.count(initdart)==1);
            Dart_const_handle d1=copy_to_origin.find(initdart)->second;
            Dart_const_handle d2=m_original_map.template beta<2>(d1);
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
            assert(copy_to_origin.count(initdart2)==1);
            Dart_const_handle d1=copy_to_origin.find(initdart2)->second;
            Dart_const_handle d2=m_original_map.template beta<2>(d1);
            if (d1<d2) {
              assert(paths.count(d1)==1);
              paths[d1].first=enddart2;
            }
            else       {
              assert(paths.count(d2)==1);
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
    bool is_edge_has_path(Dart_const_handle adart) const
    {
      Dart_const_handle opposite=m_original_map.template beta<2>(adart);
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
    (Dart_const_handle adart)
    {
      assert(!m_original_map.is_marked(adart, m_mark_T));
      assert(is_edge_has_path(adart));

      Dart_const_handle opposite=m_original_map.template beta<2>(adart);
      if (adart<opposite)
      { return paths.find(adart)->second; }

      return paths.find(opposite)->second;
    }

    Dart_const_handle get_first_dart_of_the_path(Dart_const_handle adart) const
    {
      assert(!m_original_map.is_marked(adart, m_mark_T));
      assert(is_edge_has_path(adart));

      Dart_const_handle opposite=m_original_map.template beta<2>(adart);
      if (adart<opposite)
      {
        const std::pair<Dart_const_handle, Dart_const_handle>&
            p=paths.find(adart)->second;
        return p.first;
      }

      const std::pair<Dart_const_handle, Dart_const_handle>&
          p=paths.find(opposite)->second;
      return m_map.template beta<2>(p.second);
    }

    Dart_const_handle get_second_dart_of_the_path(Dart_const_handle adart) const
    {
      assert(!m_original_map.is_marked(adart, m_mark_T));
      assert(is_edge_has_path(adart));

      Dart_const_handle opposite=m_original_map.template beta<2>(adart);
      if (adart<opposite)
      {
        const std::pair<Dart_const_handle, Dart_const_handle>&
            p=paths.find(adart)->second;
        return p.second;
      }

      const std::pair<Dart_const_handle, Dart_const_handle>&
          p=paths.find(opposite)->second;
      return m_map.template beta<2>(p.first);
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
      bool res=true;
      for (auto it=m_original_map.darts().begin(),
           itend=m_original_map.darts().end(); it!=itend; ++it)
      {
        if (!m_original_map.is_marked(it, m_mark_T))
        {
          if (!is_edge_has_path(it))
          {
            std::cout<<"ERROR: an edge that does not belong to the spanning tree"
                    <<" T has no associated path."<<std::endl;
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
            assert(dd1!=NULL);
            if (m_map.vertex_attribute(dd1)!=m_map.vertex_attribute(d2))
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
    /*
    /// @return the turn between dart number i and dart number i+1 of path.
    ///         (turn is position of the second edge in the cyclic ordering of
    ///          edges starting from the first edge around the second extremity
    ///          of the first dart)
    std::size_t next_positive_turn(const Path_on_surface<Map>& path,
                                   std::size_t i) const
    {
      // OLD return path.next_positive_turn(i);
      Dart_const_handle d1=path.get_ith_dart(i);
      Dart_const_handle d2=path.get_next_dart(i);
      assert(d1!=d2);

      if (d2==m_map.template beta<2>(d1))
      { return 0; }

      std::size_t id1=m_dart_ids.at(m_map.template beta<2>(d1));
      std::size_t id2=m_dart_ids.at(d2);

      if (id1>=m_number_of_edges)
      {
        id1-=m_number_of_edges; // id of the first dart in its own vertex
        assert(id2>=m_number_of_edges);
        id2-=m_number_of_edges; // id of the second dart in its own vertex
      }

      std::size_t res=(id1<id2?id2-id1:
                               m_number_of_edges-id1+id2);
      assert(res==path.next_positive_turn(i));
      return res;
    }

    /// Same than next_positive_turn but turning in reverse orientation around vertex.
    std::size_t next_negative_turn(const Path_on_surface<Map>& path,
                                   std::size_t i) const
    {
      // OLD return path.next_negative_turn(i);
      Dart_const_handle d1=path.get_ith_dart(i);
      Dart_const_handle d2=path.get_next_dart(i);
//      Dart_const_handle d1=m_map.template beta<2>(path.get_ith_dart(i));
//      Dart_const_handle d2=m_map.template beta<2>(path.get_next_dart(i));
      assert(d1!=d2);

      if (d2==m_map.template beta<2>(d1))
      { return 0; }

      std::size_t id1=m_dart_ids.at(m_map.template beta<2>(d1));
      std::size_t id2=m_dart_ids.at(d2);

      if (id1>=m_number_of_edges)
      {
        id1-=m_number_of_edges; // id of the first dart in its own vertex
        assert(id2>=m_number_of_edges);
        id2-=m_number_of_edges; // id of the second dart in its own vertex
      }

      std::size_t res=(id1<=id2?m_number_of_edges-id2+id1:
                               id1-id2);
      assert(res==path.next_negative_turn(i));
      return res;
    }

*/

  protected:
    const Map& m_original_map; /// The original surface; not modified
    Map m_map; /// the transformed map
    TPaths paths; /// Pair of edges associated with each edge of m_original_map
                  /// (except the edges that belong to the spanning tree T).
    std::size_t m_mark_T; /// mark each edge of m_original_map that belong to the spanning tree T
    std::size_t m_mark_L; /// mark each edge of m_original_map that belong to the dual spanning tree L
    TDartIds m_dart_ids; /// Ids of each dart of the transformed map, between 0 and n-1 (n being the number of darts)
                       /// so that darts between 0...(n/2)-1 belong to the same vertex and
                       /// d1=beta<1, 2>(d0), d2=beta<1, 2>(d1)...
                       /// The same for darts between n/2...n-1 for the second vertex
                       /// Thanks to these ids, we can compute in constant time the positive and
                       /// negative turns between two consecutive darts
    std::size_t m_number_of_edges; // number of edges in the tranformed map (==number of darts / 2)
  };
  
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_CURVE_TOPOLOGY_H //
// EOF //
