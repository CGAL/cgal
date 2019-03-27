// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_COMBINATORIAL_MAP_2_INCREMENTAL_BUILDER_H
#define CGAL_COMBINATORIAL_MAP_2_INCREMENTAL_BUILDER_H 1

#include <vector>
#include <unordered_map>
#include <cstddef>
#include <CGAL/Path_on_surface.h>

namespace CGAL {
  struct Combinatorial_map_tag;
  struct Generalized_map_tag;
  
  template<class MAP, class Combinatorial_data_structure=
           typename MAP::Combinatorial_data_structure>
  struct Map_incremental_builder_tools
  {};
  template<class CMAP>
  struct Map_incremental_builder_tools<CMAP, Combinatorial_map_tag>
  {
    static typename CMAP::Dart_handle
    add_vertex_to_face(CMAP& cmap, typename CMAP::Dart_handle prev_dart)
    {
      typename CMAP::Dart_handle res=cmap.create_dart();
      if (prev_dart!=cmap.null_handle)
      {
        cmap.template link_beta<1>(prev_dart, res);
      }
      return res;
    }
    static void add_last_vertex_to_face(CMAP& cmap,
                                        typename CMAP::Dart_handle prev_dart,
                                        typename CMAP::Dart_handle first_dart)
    { cmap.template link_beta<1>(prev_dart, first_dart); }
    static typename CMAP::Dart_handle
    add_edge_to_face(CMAP& cmap, typename CMAP::Dart_handle prev_dart)
    { return add_vertex_to_face(cmap, prev_dart); } // For CMap no difference
  };
  template<class GMAP>
  struct Map_incremental_builder_tools<GMAP, Generalized_map_tag>
  {
    static typename GMAP::Dart_handle
    add_vertex_to_face(GMAP& gmap, typename GMAP::Dart_handle prev_dart)
    {
      typename GMAP::Dart_handle res=gmap.create_dart();
      if (prev_dart!=gmap.null_handle)
      {
        //if (gmap.template is_free<0>(prev_dart))
        { // Case of the first dart of the face
          gmap.template link_alpha<0>(prev_dart, res);
          gmap.template link_alpha<1>(res, gmap.create_dart());
          res=gmap.template alpha<1>(res);
        }
        /*else
        {
          gmap.template link_alpha<0>(gmap.template alpha<1>(prev_dart), res);
          gmap.template link_alpha<1>(res, gmap.create_dart());
        }*/
      }
      assert(gmap.is_valid());
      return res;
    }
    static void add_last_vertex_to_face(GMAP& gmap,
                                        typename GMAP::Dart_handle prev_dart,
                                        typename GMAP::Dart_handle first_dart)
    {
      // here we need to create a last dart and 0-link it
      assert(prev_dart!=gmap.null_handle);
      //if (gmap.template is_free<0>(prev_dart))
      { // Case of the first dart of the face
        gmap.template link_alpha<0>(prev_dart, gmap.create_dart());
        gmap.template link_alpha<1>(gmap.template alpha<0>(prev_dart),
                                    first_dart);
        assert(gmap.is_valid());
      }
      /*else
      {
        gmap.template link_alpha<0>(gmap.template alpha<1>(prev_dart),
                                    gmap.create_dart());
      }*/
    }
    static typename GMAP::Dart_handle
    add_edge_to_face(GMAP& gmap, typename GMAP::Dart_handle prev_dart)
    {
      typename GMAP::Dart_handle res=gmap.create_dart();
      typename GMAP::Dart_handle dh2=gmap.create_dart();
      gmap.template link_alpha<0>(res, dh2);
      if (prev_dart!=gmap.null_handle)
      {
        gmap.template link_alpha<1>(res, gmap.template alpha<0>(prev_dart));
      }
      return res;
    }
  };

  // Incremental builder
  template < class CMap_ >
  class Combinatorial_map_2_incremental_builder
  {
  public:
    typedef CMap_                            CMap;
    typedef typename CMap::Dart_handle       Dart_handle;
    typedef typename CMap::Dart_const_handle Dart_const_handle;
    typedef typename CMap::size_type         size_type;

    Combinatorial_map_2_incremental_builder(CMap & acmap) :
      cmap(acmap),
      first_dart(cmap.null_handle),
      prev_dart(cmap.null_handle),
      used_add_vertex(false),
      used_add_edge(false),
      facet_started(false),
      path_started(false),
      m_cur_path(acmap)
    {}

    /// Start a new facet.
    void begin_facet()
    {
      if (facet_started)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"you try to start a facet"
                 <<" but the previous facet is not yet ended."<<std::endl;
        return;
      }

      first_dart = cmap.null_handle;
      prev_dart  = cmap.null_handle;
      facet_started=true;
      // std::cout<<"Begin facet: "<<std::flush;
    }

    /// Add vertex 'i' in the current facet.
    void add_vertex_to_facet(size_type i)
    {
      if (!facet_started)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"you try to add a vertex to a facet"
                 <<" but the facet is not yet started."<<std::endl;
        return;
      }

      if (used_add_edge)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"you cannot mix add_edge_to_facet"
                 <<" and add_vertex_to_facet for a same builder."<<std::endl;
        return;
      }
      used_add_vertex=true;
      if (i>=vertex_to_dart_map.size())
      { vertex_to_dart_map.resize(i+1); }

      // std::cout<<i<<"  "<<std::flush;
      Dart_handle cur=Map_incremental_builder_tools<CMap>::
          add_vertex_to_face(cmap, prev_dart);

      if ( prev_dart!=cmap.null_handle )
      {
        Dart_handle opposite=find_dart_between(i, prev_vertex);
        if ( opposite!=cmap.null_handle )
        {
          CGAL_assertion( cmap.template is_free<2>(opposite) );
          cmap.template set_opposite<2>(prev_dart, opposite);
        }

        add_dart_in_vertex_to_dart_map(prev_dart, prev_vertex, i);
      }
      else
      {
        first_dart  =cur;
        first_vertex=i;
      }

      prev_dart  =cur;
      prev_vertex=i;
    }

    /// Add one edge to the current facet, given by its label (any string, using minus sign for orientation)
    void add_edge_to_facet(const std::string& s)
    {
      if (!facet_started)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"you try to add an edge to a facet"
                 <<" but the facet is not yet started."<<std::endl;
        return;
      }
      if (used_add_vertex)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"you cannot mix add_edge_to_facet"
                 <<" and add_vertex_to_facet for a same builder."<<std::endl;
        return;
      }
      used_add_edge=true;

      assert(edge_label_to_dart.count(s)==0); // Since we have an orientable surface,
                                              // we cannot use a same edge twice.

      Dart_handle cur = Map_incremental_builder_tools<CMap>::
          add_edge_to_face(cmap, prev_dart);
      Dart_handle opposite=find_dart_with_label(opposite_label(s));

      if (opposite!=NULL)
      {
        CGAL_assertion(cmap.template is_free<2>(opposite));
        cmap.template set_opposite<2>(cur, opposite);
      }

      edge_label_to_dart[s]=cur;

      if (prev_dart==cmap.null_handle)
      { first_dart=cur; }

      prev_dart=cur;
    }

    /// add the given edges to the current facet
    /// s is a sequence of labels, add all the corresponding edges into the current facet.
    void add_edges_to_facet(const std::string& s)
    {
      if (!facet_started)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"you try to add edges to a facet"
                 <<" but the facet is not yet started."<<std::endl;
        return;
      }
      std::istringstream iss(s);
      for (std::string token; std::getline(iss, token, ' '); )
      { add_edge_to_facet(token); }
    }
      
    /// add one facet, s is a sequence of labels, add all the corresponding edges into a new facet.
    void add_facet(const std::string& s)
    {
      if (facet_started)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"you try to add a new facet"
                 <<" but the previous facet is not yet ended."<<std::endl;
        return;
      }
      begin_facet();
      add_edges_to_facet(s);
      end_facet();
    }

    /// End of the facet. Return the first dart of this facet.
    Dart_handle end_facet()
    {
      if (!facet_started)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"you try to end a facet"
                 <<" but the facet is not yet started."<<std::endl;
        return NULL;
      }
      CGAL_assertion( first_dart!=cmap.null_handle && prev_dart!=cmap.null_handle );
      if (used_add_vertex)
      {
        Map_incremental_builder_tools<CMap>::
            add_last_vertex_to_face(cmap, prev_dart, first_dart);

        Dart_handle opposite=find_dart_between(first_vertex, prev_vertex);
        if ( opposite!=cmap.null_handle )
        {
          CGAL_assertion( cmap.template is_free<2>(opposite) );
          cmap.template set_opposite<2>(prev_dart, opposite);
        }

        add_dart_in_vertex_to_dart_map(prev_dart, prev_vertex, first_vertex);
      }
      else
      {
        cmap.set_next(prev_dart, first_dart);
      }
      
      facet_started=false;
      return first_dart;
    }

    /// Start a new surface
    void begin_surface()
    {
      if (facet_started) { end_facet(); }

      first_dart    = cmap.null_handle;
      prev_dart     = cmap.null_handle;

      vertex_to_dart_map.clear();
      edge_label_to_dart.clear();
      
      used_add_vertex=false;
      used_add_edge=false;
    }

    /// End of the surface. Return one dart of the created surface.
    Dart_handle end_surface()
    { return first_dart; }

    /// Start a path on the surface
    void begin_path()
    {
      if (path_started)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"you try to start a path"
                 <<" but the previous path is not yet ended."<<std::endl;
        return;
      }
      path_started=true;
      m_first_path_vertex=true;
      m_cur_path.clear();
    }

    /// Add vertex with id i at the end of the current path
    void add_vertex_to_path(size_type i)
    {
      if (!path_started)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"you try to add a vertex to a path"
                 <<" but the path is not yet started."<<std::endl;
        return;
      }
      if (used_add_edge)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"you cannot use add_vertex_to_path with a "
                 <<"combinatorial map built by adding edges."<<std::endl;
        return;
      }

      if (m_first_path_vertex)
      {
        path_prev_vertex=i;
        m_first_path_vertex=false;
      }
      else
      {
        Dart_const_handle dh=find_dart_between(path_prev_vertex, i);
        if (dh==NULL)
        {
          std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                   <<"edge ("<<path_prev_vertex<<", "<<i<<") does not exists "
                   <<"and thus cannot be added in the path."<<std::endl;
          return;
        }
        CGAL_assertion(m_cur_path.can_be_pushed(dh));
        m_cur_path.push_back(dh);
        path_prev_vertex=i;
      }
    }

    /// Add edge labeled e at the end of the current path
    void add_edge_to_path(const std::string& e)
    {
      if (!path_started)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"you try to add an edge to a path"
                 <<" but the path is not yet started."<<std::endl;
        return;
      }
      if (used_add_vertex)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"you cannot use add_edge_to_path with a "
                 <<"combinatorial map built by adding vertices."<<std::endl;
        return;
      }

      Dart_const_handle dh=find_dart_with_label(e);
      if (dh==NULL)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"edge labeled ("<<e<<") does not exists "
                 <<"and thus cannot be added in the path."<<std::endl;
        return;
      }

      if (!m_cur_path.can_be_pushed(dh))
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"edge labeled ("<<e<<") is not adjacent to the previous "
                 <<"edge of the path and thus cannot be added."<<std::endl;
        return;
      }
      m_cur_path.push_back(dh);
    }

    /// End the current path
    CGAL::Path_on_surface<CMap> end_path()
    {
      if (!path_started)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: "
                 <<"you try to end a path"
                 <<" but the path is not yet started."<<std::endl;
        return m_cur_path;
      }
      path_started=false;
      return m_cur_path;
    }

    /// A shortcut allowing to create a path directly with a sequence
    /// of vertex ids, if the map was built by adding vertices
    /// or edge labels, if the map was built by adding edges
    CGAL::Path_on_surface<CMap> create_path(const std::string& s)
    {
      begin_path();

      std::istringstream iss(s);
      for (std::string token; std::getline(iss, token, ' '); )
      { add_edge_to_path(token); }

      return end_path();
    }

  protected:
    /// @return dart between vertices i and j, NULL if this dart does not exist
    ///   can be used only if the current map is created by adding vertices
    Dart_handle find_dart_between(size_type i, size_type j)
    {
      if (i>=vertex_to_dart_map.size()) { return NULL; }

      typename std::vector<std::pair<Dart_handle, size_type> >::reverse_iterator
        it(vertex_to_dart_map[i].rbegin());
      typename std::vector<std::pair<Dart_handle, size_type> >::reverse_iterator
        itend(vertex_to_dart_map[i].rend());

      for ( ; it!=itend; ++it )
      {
        if ( it->second==j ) return (it->first);
      }
      return cmap.null_handle;
    }

    /// Add dart adart between vertices i and j.
    ///   can be used only if the current map is created by adding vertices
    void add_dart_in_vertex_to_dart_map(Dart_handle adart, size_type i, size_type j)
    {
      CGAL_assertion( adart!=cmap.null_handle );
      vertex_to_dart_map[i].push_back(std::make_pair(adart, j));
    }

    /// @return opposite label of label s
    ///    (i.e. add/remove - depending if s is positive or negative)
    std::string opposite_label(const std::string & s)
    {
      CGAL_assertion(!s.empty());
      if (s[0]=='-')
      { return s.substr(1, std::string::npos); }

      return std::string("-")+s;
    }
    
    /// @return dart with the given label, NULL if this dart does not exist
    ///   can be used only if the current map is created by adding edges
    Dart_handle find_dart_with_label(const std::string & s)
    {
      auto ite=edge_label_to_dart.find(s);
      if (ite==edge_label_to_dart.end())
      { return NULL; }

      return ite->second;
    }

  protected:
    // For vertex number i, vector of all vertices linked to i with an edge (dart of the edge and index
    // of the second vertex). TODO use a std::unordered_map instead of the second vector
    typedef std::vector<std::vector<std::pair<Dart_handle, size_type> > > Vertex_to_dart_array;
    Vertex_to_dart_array vertex_to_dart_map;

    // For edge label, its corresponding dart. For an edge a -a, stores its two darts in two different entries.
    std::unordered_map<std::string, Dart_handle> edge_label_to_dart;
    
    CMap&       cmap;
    Dart_handle first_dart;
    Dart_handle prev_dart;
    size_type   first_vertex;
    size_type   prev_vertex;
    bool        used_add_vertex;
    bool        used_add_edge;
    bool        facet_started;
    bool        path_started;
    CGAL::Path_on_surface<CMap> m_cur_path;
    size_type   path_prev_vertex;
    bool m_first_path_vertex;
  };

} //namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_2_INCREMENTAL_BUILDER_H //
// EOF //
