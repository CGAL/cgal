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
#ifndef CGAL_COMBINATORIAL_MAP_FUNCTIONALITIES_H
#define CGAL_COMBINATORIAL_MAP_FUNCTIONALITIES_H 1

#include <stack>
#include <CGAL/Union_find.h>
#include <boost/unordered_map.hpp>
#include <CGAL/Random.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Combinatorial_map_basic_operations.h>

namespace CGAL {
  
  template<typename Map>
  class Combinatorial_map_tools
  {
  public:
    typedef typename Map::Dart_handle Dart_handle;
    typedef typename Map::Dart_const_handle Dart_const_handle;
    typedef CGAL::Union_find<Dart_handle> UFTree;
    typedef typename UFTree::handle UFTree_handle;
    
    typedef boost::unordered_map<Dart_const_handle,
                      std::pair<Dart_const_handle, Dart_const_handle> > TPaths;

    Combinatorial_map_tools(Map& amap) : m_original_map(amap)
    {
      if (!m_map.is_without_boundary(1))
      {
        std::cerr<<"ERROR: the given amap has 1-boundaries; "
                 <<"such a surface is not possible to process here."
                 <<std::endl;
      }
      if (!m_map.is_without_boundary(2))
      {
        std::cerr<<"ERROR: the given amap has 2-boundaries; "
                 <<"which are not yet considered (but this will be done later)."
                 <<std::endl;
      }
 
      // The mapping between darts of the original map into the copied map.
      boost::unordered_map<Dart_const_handle, Dart_handle> origin_to_copy;

      // We copy the original map, while keeping a mapping between darts.
      m_map.copy(m_original_map, &origin_to_copy);

      // The mapping between darts of the copy into darts of the original map.
      boost::unordered_map<Dart_handle, Dart_const_handle> copy_to_origin;
      for (auto it=origin_to_copy.begin(); it!=origin_to_copy.end(); ++it)
      { copy_to_origin[it->second]=it->first; }

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

      /* std::cout<<"Number of darts in m_map: "<<m_map.number_of_darts()
              <<"; number of darts in origin_to_copy: "<<origin_to_copy.size()
             <<"; number of darts in copy_to_origin: "<<copy_to_origin.size()
            <<std::endl; */

      /* std::cout<<"Paths are all valid 1 ? "<<(are_paths_valid()?"YES":"NO")
               <<std::endl; */

      // 3) We simplify m_map in a surface with only one face
      surface_simplification_in_one_face(origin_to_copy, copy_to_origin);

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

      // 4) And we quadrangulate the face
      surface_quadrangulate();

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

      m_map.display_darts(std::cout);

      assert(are_paths_valid());
    }
    
    ~Combinatorial_map_tools()
    {
      m_original_map.free_mark(m_mark_T);
      m_original_map.free_mark(m_mark_L);
    }

    Path_on_surface<Map> transform_original_path_into_quad_surface
    (const Path_on_surface<Map>& path)
    {
      Path_on_surface<Map> res(m_map);
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

  protected:
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

    // Mark the edge containing adart in the given map.
    void mark_edge(const Map& amap, Dart_const_handle adart, std::size_t amark)
    {
      amap.mark(amap.template beta<2>(adart), amark);
      amap.mark(adart, amark);
    }

    // Erase the edge given by adart (which belongs to the map m_map) from the
    // associative array copy_to_origin, and erase the corresponding edge
    // (which belongs to the map m_original_map) from the array origin_to_copy
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

    // Step 1) Transform m_map into an equivalent surface having only one
    // vertex. All edges contracted during this step belong to the spanning
    // tree T, and thus corresponding edges in m_original_map are marked.
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

    // Step 2) Compute, for each edge of m_original_map not in the spanning
    // tree T, the pair of darts of the edge in m_copy. This pair of edges
    // will be updated later (in surface_simplification_in_one_face() and in
    // surface_quadrangulate() )
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

    // Step 3) Transform the 2-map into an equivalent surface having only
    // one vertex. All edges removed during this step belong to the
    // dual spanning tree L (spanning tree of the dual 2-map).
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
      }
      else
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

      m_map.set_automatic_attributes_management(true);
      m_map.free_mark(toremove);
    }

    // Step 4) quadrangulate the surface.
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

  protected:
    const Map& m_original_map; // The original surface; not modified
    Map m_map; // the transformed map
    TPaths paths; // Pair of edges associated with each edge of m_original_map
                  // (except the edges that belong to the spanning tree T).
    std::size_t m_mark_T; // mark each edge of m_original_map that belong to the spanning tree T
    std::size_t m_mark_L; // mark each edge of m_original_map that belong to the dual spanning tree L
  };
  
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_FUNCTIONALITIES_H //
// EOF //
