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
#ifndef CGAL_LINEAR_CELL_COMPLEX_INCREMENTAL_BUILDER_H
#define CGAL_LINEAR_CELL_COMPLEX_INCREMENTAL_BUILDER_H 1

#include <vector>
#include <cstddef>

namespace CGAL {
  template<class LCC, class Combinatorial_data_structure=
           typename LCC::Combinatorial_data_structure>
  struct Add_vertex_to_face
  {
    static typename LCC::Dart_handle run(LCC&,
                                         typename LCC::Vertex_attribute_handle,
                                         typename LCC::Dart_handle)
    {}
  };
  template<class LCC>
  struct Add_vertex_to_face<LCC, Combinatorial_map_tag>
  {
    static typename LCC::Dart_handle run(LCC& lcc,
                                         typename LCC::Vertex_attribute_handle vh,
                                         typename LCC::Dart_handle prev_dart)
    {
      typename LCC::Dart_handle res=lcc.create_dart(vh);
      if (prev_dart!=lcc.null_handle)
      {
        lcc.template link_beta<1>(prev_dart, res);
      }
      return res;
    }
    static void run_for_last(LCC&,
                             typename LCC::Vertex_attribute_handle,
                             typename LCC::Dart_handle)
    { // here nothing to do, all darts were already created.
    }
  };
  template<class LCC>
  struct Add_vertex_to_face<LCC, Generalized_map_tag>
  {
    static typename LCC::Dart_handle run(LCC& lcc,
                                         typename LCC::Vertex_attribute_handle vh,
                                         typename LCC::Dart_handle prev_dart)
    {
      typename LCC::Dart_handle res=lcc.create_dart(vh);
      if (prev_dart!=lcc.null_handle)
      {
        lcc.template link_alpha<0>(prev_dart, res);
        lcc.template link_alpha<1>(res, lcc.create_dart(vh));
        res=lcc.template alpha<1>(res);
      }
      return res;
    }
    static void run_for_last(LCC& lcc,
                             typename LCC::Vertex_attribute_handle vh,
                             typename LCC::Dart_handle prev_dart)
    {
      // here we need to create a last dart and 0-link it
      assert(prev_dart!=lcc.null_handle);
      lcc.template link_alpha<0>(prev_dart, lcc.create_dart(vh));
    }
  };

  // Incremental builder
  template < class LCC_ >
  class Linear_cell_complex_incremental_builder_3
  {
  public:
    typedef LCC_ LCC;
    typedef typename LCC::Dart_handle             Dart_handle;
    typedef typename LCC::Vertex_attribute_handle Vertex_attribute_handle;
    typedef typename LCC::Point                   Point_3;
    typedef typename LCC::size_type               size_type;

    Linear_cell_complex_incremental_builder_3(LCC & alcc) :
      lcc(alcc)
    {}

    Vertex_attribute_handle add_vertex(const Point_3& p)
    {
      Vertex_attribute_handle res = lcc.create_vertex_attribute(p);
      vertex_map.push_back(res);
      vertex_to_dart_map.push_back(std::vector<Dart_handle>());
      ++new_vertices;
      return res;
    }

    void begin_facet()
    {
      first_dart = lcc.null_handle;
      prev_dart  = lcc.null_handle;
      // std::cout<<"Begin facet: "<<std::flush;
    }

    void add_vertex_to_facet(size_type i)
    {
      CGAL_assertion( i<new_vertices );
      // std::cout<<i<<"  "<<std::flush;
      Dart_handle cur = Add_vertex_to_face<LCC>::
          run(lcc, vertex_map[i], prev_dart);

      if ( prev_dart!=lcc.null_handle )
      {
        Dart_handle opposite=
            find_dart_between(i,lcc.vertex_attribute(prev_dart));
        if ( opposite!=lcc.null_handle )
        {
          CGAL_assertion( lcc.template is_free<2>(opposite) );
          lcc.template set_opposite<2>(prev_dart, opposite);
        }

        add_dart_in_vertex_to_dart_map(prev_dart, prev_vertex);
      }
      else
      {
        first_dart   = cur;
        first_vertex = i;
      }

      prev_dart   = cur;
      prev_vertex = i;
    }

    // End of the facet. Return the first dart of this facet.
    Dart_handle end_facet()
    {
      CGAL_assertion( first_dart!=lcc.null_handle && prev_dart!=lcc.null_handle );

      Add_vertex_to_face<LCC>::run_for_last(lcc, vertex_map[first_vertex],
                                            prev_dart);

      lcc.set_next(prev_dart, first_dart);

      Dart_handle opposite =
        find_dart_between(first_vertex,lcc.vertex_attribute(prev_dart));
      if ( opposite!=lcc.null_handle )
      {
        CGAL_assertion( lcc.template is_free<2>(opposite) );
        lcc.template set_opposite<2>(prev_dart, opposite);
      }

      add_dart_in_vertex_to_dart_map(prev_dart, prev_vertex);

      return first_dart;
    }

    void begin_surface( std::size_t v, std::size_t /*f*/, std::size_t /*h*/)
    {
      new_vertices  = 0;
      first_dart    = lcc.null_handle;
      prev_dart     = lcc.null_handle;
      vertex_map.clear();
      vertex_to_dart_map.clear();
      vertex_map.reserve(v);
      vertex_to_dart_map.reserve(v);
      // lcc.reserve(v,h);
    }

    // End of the surface. Return one dart of the created surface.
    Dart_handle end_surface()
    { return first_dart; }

  protected:

    Dart_handle find_dart_between(size_type i, Vertex_attribute_handle vh)
    {
      typename std::vector<Dart_handle>::reverse_iterator
        it(vertex_to_dart_map[i].rbegin());
      typename std::vector<Dart_handle>::reverse_iterator
        itend(vertex_to_dart_map[i].rend());

      for ( ; it!=itend; ++it )
      {
        if ( lcc.vertex_attribute(lcc.next(*it))==vh ) return (*it);
      }
      return lcc.null_handle;
    }

    void add_dart_in_vertex_to_dart_map( Dart_handle adart, size_type i )
    {
      CGAL_assertion( adart!=lcc.null_handle );
      vertex_to_dart_map[i].push_back(adart);
    }

  private:
    std::vector<Vertex_attribute_handle> vertex_map;
    std::vector<std::vector<Dart_handle> > vertex_to_dart_map;

    LCC&                      lcc;
    Dart_handle               first_dart;
    Dart_handle               prev_dart;
    size_type                 first_vertex;
    size_type                 prev_vertex;
    size_type                 new_vertices;
  };

} //namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_INCREMENTAL_BUILDER_H //
// EOF //
