// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_COMBINATORIAL_MAP_INCREMENTAL_BUILDER_H
#define CGAL_COMBINATORIAL_MAP_INCREMENTAL_BUILDER_H 1

#include <vector>
#include <unordered_map>
#include <cstddef>

namespace CGAL {
  struct Combinatorial_map_tag;
  struct Generalized_map_tag;
  
  template<class CMAP, class Combinatorial_data_structure=
           typename CMAP::Combinatorial_data_structure>
  struct Add_vertex_to_face
  {
    static typename CMAP::Dart_handle run(CMAP&,
                                         typename CMAP::Dart_handle)
    {}
  };
  template<class CMAP>
  struct Add_vertex_to_face<CMAP, Combinatorial_map_tag>
  {
    static typename CMAP::Dart_handle run(CMAP& cmap,
                                         typename CMAP::Dart_handle prev_dart)
    {
      typename CMAP::Dart_handle res=cmap.create_dart();
      if (prev_dart!=cmap.null_handle)
      {
        cmap.template link_beta<1>(prev_dart, res);
      }
      return res;
    }
    static void run_for_last(CMAP&,
                             typename CMAP::Dart_handle)
    { // here nothing to do, all darts were already created.
    }
  };
  template<class CMAP>
  struct Add_vertex_to_face<CMAP, Generalized_map_tag>
  {
    static typename CMAP::Dart_handle run(CMAP& cmap,
                                         typename CMAP::Dart_handle prev_dart)
    {
      typename CMAP::Dart_handle res=cmap.create_dart();
      if (prev_dart!=cmap.null_handle)
      {
        cmap.template link_alpha<0>(prev_dart, res);
        cmap.template link_alpha<1>(res, cmap.create_dart());
        res=cmap.template alpha<1>(res);
      }
      return res;
    }
    static void run_for_last(CMAP& cmap,
                             typename CMAP::Dart_handle prev_dart)
    {
      // here we need to create a last dart and 0-link it
      assert(prev_dart!=cmap.null_handle);
      cmap.template link_alpha<0>(prev_dart, cmap.create_dart());
    }
  };

  // Incremental builder
  template < class CMap_ >
  class Combinatorial_map_incremental_builder
  {
  public:
    typedef CMap_                      CMap;
    typedef typename CMap::Dart_handle Dart_handle;
    typedef typename CMap::size_type   size_type;

    Combinatorial_map_incremental_builder(CMap & acmap) :
      cmap(acmap),
      first_dart(cmap.null_handle),
      prev_dart(cmap.null_handle),
      used_add_vertex(false),
      used_add_edge(false)
    {}

    void begin_facet()
    {
      first_dart = cmap.null_handle;
      prev_dart  = cmap.null_handle;
      // std::cout<<"Begin facet: "<<std::flush;
    }

    void add_vertex_to_facet(size_type i)
    {
      if (used_add_edge)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: you cannot mix add_edge_to_facet"
                 <<" and add_vertex_to_facet for a same builder."<<std::endl;
        return;
      }
      used_add_vertex=true;
        
      CGAL_assertion( i<new_vertices );
      // std::cout<<i<<"  "<<std::flush;
      Dart_handle cur = Add_vertex_to_face<CMap>::run(cmap, prev_dart);

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
        first_dart   = cur;
        first_vertex = i;
      }

      prev_dart   = cur;
      prev_vertex = i;
    }

    // Add one edge to the current facet, given by its label (any string, using minus sign for orientation) 
    void add_edge_to_facet(const std::string& s)
    {
      if (used_add_vertex)
      {
        std::cerr<<"Combinatorial_map_incremental_builder ERROR: you cannot mix add_edge_to_facet"
                 <<" and add_vertex_to_facet for a same builder."<<std::endl;
        return;
      }
      used_add_edge=true;

      assert(edge_label_to_dart.count(s)==0); // Since we have an orientable surface,
                                              // we cannot use a same edge twice.

      Dart_handle cur = Add_vertex_to_face<CMap>::run(cmap, prev_dart);
      auto ite=edge_label_to_dart.find(opposite_label(s));

      if (ite!=edge_label_to_dart.end())
      {
        CGAL_assertion( cmap.template is_free<2>(ite->second) );
        cmap.template set_opposite<2>(cur, ite->second);
      }

      edge_label_to_dart[s]=cur;

      if (prev_dart==cmap.null_handle)
      { first_dart=cur; }

      prev_dart=cur;
    }

    // add one facet, s is a sequence of labels, add all the corresponding edges into the current facet.
    void add_facet(const std::string& s)
    {
      begin_facet();
      
      std::istringstream iss(s);
      for (std::string token; std::getline(iss, token, ' '); )
      { add_edge_to_facet(token); }

      end_facet();
    }

    // End of the facet. Return the first dart of this facet.
    Dart_handle end_facet()
    {
      CGAL_assertion( first_dart!=cmap.null_handle && prev_dart!=cmap.null_handle );

      if (used_add_vertex)
      {
        Add_vertex_to_face<CMap>::run_for_last(cmap, prev_dart);
        cmap.set_next(prev_dart, first_dart);

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
      
      return first_dart;
    }

    void begin_surface( std::size_t v, std::size_t /*f*/, std::size_t h)
    {
      new_vertices  = 0;
      first_dart    = cmap.null_handle;
      prev_dart     = cmap.null_handle;

      vertex_to_dart_map.clear();
      vertex_to_dart_map.reserve(v);

      edge_label_to_dart.clear();
      edge_label_to_dart.reserve(h/2);
      
      used_add_vertex=false;
      used_add_edge=false;
      
      // cmap.reserve(v,h);
    }

    // End of the surface. Return one dart of the created surface.
    Dart_handle end_surface()
    { return first_dart; }

  protected:

    Dart_handle find_dart_between(size_type i, size_type j)
    {
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

    void add_dart_in_vertex_to_dart_map(Dart_handle adart, size_type i, size_type j)
    {
      CGAL_assertion( adart!=cmap.null_handle );
      vertex_to_dart_map[i].push_back(std::make_pair(adart, j));
    }

    std::string opposite_label(const std::string & s)
    {
      assert(!s.empty());
      if (s[0]=='-')
      { return s.substr(1, std::string::npos); }

      return std::string("-")+s;
    }
    
  private:
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
    size_type   new_vertices;
    bool        used_add_vertex;
    bool        used_add_edge;
  };

} //namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_INCREMENTAL_BUILDER_H //
// EOF //
