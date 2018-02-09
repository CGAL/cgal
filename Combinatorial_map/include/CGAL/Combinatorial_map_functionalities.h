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
#include<CGAL/Random.h>

namespace CGAL {

template<typename Map>
class Path_on_surface
{
public:
  typedef typename Map::Dart_handle Dart_handle;
  typedef typename Map::Dart_const_handle Dart_const_handle;

  Path_on_surface(const Map& amap) : m_map(amap)
  {}

  unsigned int length() const
  { return m_path.size(); }

  Dart_const_handle get_ith_dart(unsigned int i) const
  {
    assert(i<m_path.size());
    return m_path[i];
  }
  
  // @return true iff the path is valid; i.e. a sequence of edges two by
  //              two adjacent.
  bool is_valid() const
  {
    for (unsigned int i=1; i<m_path.size(); ++i)
    {
      Dart_const_handle pend=m_map.other_extremity(m_path[i-1]);
      if (pend==Map::null_handle) { return false; }

      if (!CGAL::template belong_to_same_cell<Map,0>(m_map, m_path[i], pend))
      { return false; }
    }
    return true;
  }

  // @return true iff the path is empty
  bool is_empty() const
  { return m_path.empty(); }

  // @return true iff the path is closed (i.e. the second extremity of the
  //              last dart of the path is the same vertex than the one of the
  //              first dart of the path.
  bool is_closed() const
  {
    if (is_empty()) { return false; } // or true by vacuity ?
    if (!is_valid()) { return false; } // Interest ??

    Dart_const_handle pend=m_map.other_extremity(m_path.back());
    if (pend==Map::null_handle) { return false; }

    return CGAL::belong_to_same_cell<Map,0>(m_map, m_path[0], pend);
  }

  // @return true iff the path does not pass twice through a same edge
  //              or a same vertex.
  bool is_simple() const
  {
    typename Map::size_type markvertex=m_map.get_new_mark();
    typename Map::size_type markedge=m_map.get_new_mark();

    bool res=true;
    unsigned int i=0;
    for (i=0; res && i<m_path.size(); ++i)
    {
      if (m_map.is_marked(m_path[i], markvertex)) res=false;
      if (m_map.is_marked(m_path[i], markedge)) res=false;

      CGAL::mark_cell<Map, 0>(m_path[i], markvertex);
      CGAL::mark_cell<Map, 1>(m_path[i], markedge);
    }

    i=0;
    while(m_map.number_of_marked_darts(markedge)>0)
    {
      assert(i<m_path.size());
      CGAL::unmark_cell<Map, 0>(m_path[i], markvertex);
      CGAL::unmark_cell<Map, 1>(m_path[i], markedge);
      ++i;
    }

    m_map.free_mark(markvertex);
    m_map.free_mark(markedge);

    return res;
  }

  // Generate a random path with about percent % of edge
  void generate_random_path(unsigned int percent, CGAL::Random& random)
  {
    m_path.clear();
    unsigned int length=((m_map.number_of_darts()/2)*percent)/100;
    for (unsigned int i=0; i<length; ++i)
    { extend_randomly(random); }
  }
  void generate_random_path(unsigned int percent)
  {
    CGAL::Random random;
    generate_random_path(percent, random);
  }
  
  bool extend_randomly(CGAL::Random& random, bool allow_half_turn=false)
  {
    if (m_path.empty())
    {
      unsigned int index=random.get_int(0, m_map.darts().capacity());
      while (!m_map.darts().is_used(index))
      {
         ++index;
        if (index==m_map.darts().capacity()) index=0;
      }
      m_path.push_back(m_map.darts().iterator_to(m_map.darts()[index]));
      return true;
    }

    Dart_const_handle pend=m_map.opposite(m_path.back());
    if (pend==Map::null_handle)
    {
      if (!m_map.template is_free<1>(m_path.back()))
      {
        m_path.push_back(m_map.template beta<1>(m_path.back()));
        return true;
      }
      else { return false; }
    }

    typename Map::template Dart_of_cell_range<0>::const_iterator
        it=m_map.template darts_of_cell<0>(pend).begin();

     unsigned int index=random.get_int((allow_half_turn?0:1), m_map.template darts_of_cell<0>
                                       (pend).size());
     for(unsigned int i=0; i<index; ++i)
     { ++it; }

     assert(allow_half_turn || it!=pend);
     
     m_path.push_back(it);
     return true;
  }

protected:
  const Map& m_map;
  std::vector<Dart_const_handle> m_path;
};


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
        std::cerr<<"ERROR: the given amap has 1-boundaries; such a surface is not possible to process here."
                 <<std::endl;
      }
      if (!m_map.is_without_boundary(2))
      {
        std::cerr<<"ERROR: the given amap has 2-boundaries; which are not yet considered (but this will be done later)."
                 <<std::endl;
      }

 
      // A mapping between darts of the original map into the transformed map.
      boost::unordered_map<Dart_const_handle, Dart_handle> dart_mapping;

      // We copy the original map, while keeping a mapping between darts.
      m_map.copy(m_original_map, &dart_mapping);

      // We simplify m_map in a surface with only one vertex
      surface_simplification_in_one_vertex();
      std::cout<<"All non loop contracted: ";
      m_map.display_characteristics(std::cout) << ", valid=" 
                                               << m_map.is_valid() << std::endl;      
      
      // Now we compute each length two path associated with each edge that does
      // not belong to the spanning tree (which are thus all the survival edges).
      compute_length_two_paths(dart_mapping);

      // We simplify m_map in a surface with only one face
      surface_simplification_in_one_face();
      std::cout<<"All faces merges: ";
      m_map.display_characteristics(std::cout) << ", valid=" 
                                               << m_map.is_valid() << std::endl;

      // And we quadrangulate the face 
      surface_quadrangulate();
      std::cout<<"After quadrangulation: ";
      m_map.display_characteristics(std::cout) << ", valid=" 
                                               << m_map.is_valid() << std::endl;
    }
    
    void initialize_vertices(UFTree& uftrees,
                             boost::unordered_map<Dart_const_handle, UFTree_handle>& vertices)
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
                          boost::unordered_map<Dart_const_handle, UFTree_handle>& faces)
    {
      uftrees.clear();
      faces.clear();

      typename Map::size_type treated=m_map.get_new_mark();
      for (typename Map::Dart_range::iterator it=m_map.darts().begin(), itend=m_map.darts().end();
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

    // Transform the 2-map into an equivalent surface having only one vertex.
    // All edges removed during this step belong to the spanning tree L
    // (spanning tree of the dual 2-map).
    void surface_simplification_in_one_face()
    {
      UFTree uftrees; // uftree of faces; one tree for each face, contains one dart of the face
      boost::unordered_map<Dart_const_handle, UFTree_handle> faces;
      initialize_faces(uftrees, faces);

      m_map.set_automatic_attributes_management(false);

      typename Map::size_type treated=m_map.get_new_mark();
      Dart_handle currentdart=NULL, oppositedart=NULL;
      
      for (typename Map::Dart_range::iterator it=m_map.darts().begin(),
           itend=m_map.darts().end(); it!=itend;)
      {
        currentdart=it++;
        if (!m_map.is_marked(currentdart, treated))
        {
          if (m_map.template is_free<2>(currentdart))
          { m_map.mark(currentdart, treated); }
          else
          {
            oppositedart=m_map.template beta<2>(currentdart);
            m_map.mark(currentdart, treated);
            m_map.mark(oppositedart, treated);
            
            // We remove dangling edges and degree two edges.
            // The two first tests allow to keep isolated edges (case of spheres)
            if ((m_map.template beta<0>(currentdart)!=oppositedart ||
                 m_map.template beta<1>(currentdart)!=oppositedart) &&
                get_uftree(uftrees, faces, currentdart)!=
                get_uftree(uftrees, faces, oppositedart))
            {
              uftrees.unify_sets(get_uftree(uftrees, faces, currentdart),
                                 get_uftree(uftrees, faces, oppositedart));

              if (m_map.is_marked(it, treated))
              { ++it; }
              
              // 
              std::pair<Dart_const_handle, Dart_const_handle>&
                p=(paths.find((currentdart<oppositedart)?currentdart:oppositedart))->second;
              p.first=m_map.template beta<0,2>(p.first);
              p.second=m_map.template beta<0,2>(p.second);
              
              // TODO LATER (?) OPTIMIZE AND REPLACE THE REMOVE_CELL CALL BY THE MODIFICATION BY HAND
              // OR DEVELOP A SPECIALIZED VERSION OF REMOVE_CELL
              m_map.template remove_cell<1>(currentdart);
            }
          }
        }
      }
      
      m_map.set_automatic_attributes_management(true);
      m_map.free_mark(treated);
    }

    // Transform the 2-map into an equivalent surface having only one vertex.
    // All edges contracted during this step belong to the spanning tree T.
    void surface_simplification_in_one_vertex()
    {
      UFTree uftrees; // uftree of vertices; one tree for each vertex, contains one dart of the vertex
      boost::unordered_map<Dart_const_handle, UFTree_handle> vertices;
      initialize_vertices(uftrees, vertices);

      m_map.set_automatic_attributes_management(false);

      for (typename Map::Dart_range::iterator it=m_map.darts().begin(),
           itend=m_map.darts().end(); it!=itend; ++it)
      {
        if (get_uftree(uftrees, vertices, it)!=
            get_uftree(uftrees, vertices, m_map.template beta<2>(it)))
        {
          uftrees.unify_sets(get_uftree(uftrees, vertices, it),
                             get_uftree(uftrees, vertices, m_map.template beta<2>(it)));
          m_map.template contract_cell<1>(it);
        }
      }

      m_map.set_automatic_attributes_management(true);
    }
    
    void surface_quadrangulate()
    {
      // Here the map has only one face and one vertex.
      typename Map::size_type oldedges=m_map.get_new_mark();
      m_map.negate_mark(oldedges); // now all edges are marked
      
      // 1) We insert a vertex in the face (note that all points have the same geometry).
      //    New edges created by the operation are not marked.
      m_map.insert_point_in_cell_2(m_map.darts().begin(), m_map.point(m_map.darts().begin()));

      // 2) We update the pair of darts
      for (typename TPaths::iterator itp=paths.begin(), itpend=paths.end(); itp!=itpend; ++itp)
      {
        std::pair<Dart_const_handle, Dart_const_handle>& p=itp->second;
        p.first=m_map.template beta<0, 2>(p.first);
        p.second=m_map.template beta<0>(p.second);
      }

      // 3) We remove all the old edges.
      for (typename Map::Dart_range::iterator it=m_map.darts().begin(), itend=m_map.darts().end();
           it!=itend; ++it)
      {
        if (m_map.is_marked(it, oldedges))
        { m_map.template remove_cell<1>(it); }
      }

      m_map.free_mark(oldedges);
    }

    // Compute, for each edge not in the spanning tree, the pair of darts adjacent
    // to it
    // dart_mapping is a mapping between darts of m_original_map and darts in m_map
    void compute_length_two_paths(boost::unordered_map<Dart_const_handle, Dart_handle>& dart_mapping)
    {
      paths.clear();

      for (typename Map::Dart_range::const_iterator it=m_original_map.darts().begin(),
           itend=m_original_map.darts().end(); it!=itend; ++it)
      {
        if (!m_map.is_dart_used(dart_mapping[it]))
        {  // This dart was removed during surface_simplification_in_one_vertex()
          dart_mapping.erase(dart_mapping[it]);
        }
        else
        {
          if (m_map.template is_free<2>(it) ||
              it<m_map.template beta<2>(it))
          {
            paths[it]=std::make_pair(it, m_map.template beta<2>(it));
            //          paths[it]=std::make_pair(m_map.template beta<0>(it),
            //                                   m_map.template beta<2,0>(it));
          }
        }
      }
    }

  protected:
    const Map& m_original_map; // The original surface; not modified
    Map m_map; // the transformed map
    TPaths paths;
  };
  
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_FUNCTIONALITIES_H //
// EOF //
