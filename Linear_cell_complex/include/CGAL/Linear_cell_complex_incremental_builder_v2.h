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
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_LINEAR_CELL_COMPLEX_INCREMENTAL_BUILDER_V2_H
#define CGAL_LINEAR_CELL_COMPLEX_INCREMENTAL_BUILDER_V2_H 1

#include <CGAL/Linear_cell_complex.h>
#include <vector>
#include <cstddef>

namespace CGAL {
/** This second version of incremental builder will create a LCC in similar
 * way that what is done for Polyhedron_3 and Surface_mesh:
 * first we need an LCC with 2-attributes; second darts are always created
 * by pair (two consecutives darts); third 2-attributes are created and
 * associated with faces, non existing dart having no 2-attribute associated.
 */
  template < class LCC_ >
  class Linear_cell_complex_incremental_builder_3_v2
  {
  public:
    typedef LCC_ LCC;
    typedef typename LCC::Dart_handle             Dart_handle;
    typedef typename LCC::Vertex_attribute_handle Vertex_attribute_handle;
    typedef typename LCC::Point                   Point_3;
    typedef typename LCC::size_type               size_type;

    Linear_cell_complex_incremental_builder_3_v2(LCC & alcc) :
      lcc(alcc)
    {}

    Vertex_attribute_handle add_vertex (const Point_3& p)
    {
      Vertex_attribute_handle res = lcc.create_vertex_attribute(p);
      vertex_map.push_back(res);
      vertex_to_dart_map.push_back(std::vector<Dart_handle>());
      ++new_vertices;
      return res;
    }

    void begin_facet()
    {
      // std::cout<<"Begin facet: "<<std::flush;
      first_dart = lcc.null_handle;
      prev_dart  = lcc.null_handle;
      face_attrib = lcc.template create_attribute<2>();
      begin_face = true;
    }

    void add_vertex_to_facet(size_type i)
    {
      // std::cout<<i<<"  "<<std::flush;
      CGAL_assertion( i<new_vertices );

      if ( !begin_face )
      {
        Dart_handle cur=find_dart_between(prev_vertex, i);

        if (cur==lcc.null_handle)
        {
          cur = lcc.create_dart(vertex_map[i]);
          Dart_handle opposite=lcc.create_dart(vertex_map[prev_vertex]);
          lcc.template basic_link_beta_for_involution<2>(cur, opposite);
          add_dart_in_vertex_to_dart_map( opposite, i );
        }

        lcc.template set_dart_attribute<2>(cur, face_attrib);

        if (first_dart==lcc.null_handle)
        { first_dart=cur; }
        else
        { lcc.basic_link_beta_1(prev_dart, cur); }

        prev_dart = cur;
      }
      else
      {
        first_vertex = i;
        begin_face = false;
      }

      prev_vertex = i;
    }

    // End of the facet. Return the first dart of this facet.
    Dart_handle end_facet()
    {
      // std::cout<<"  end facet."<<std::endl;
      CGAL_assertion( first_dart!=lcc.null_handle &&
                      prev_dart!=lcc.null_handle );

      add_vertex_to_facet(first_vertex);
      lcc.basic_link_beta_1(prev_dart, first_dart);
      face_attrib=NULL;
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
    {
      unsigned int nb=0;

      for (typename LCC::Dart_range::iterator it=lcc.darts().begin(),
           itend=lcc.darts().end(); it!=itend; ++it)
      {
        if (lcc.template attribute<2>(it)==NULL &&
            lcc.template is_free<1>(it))
        {
          Dart_handle other=it;
          do
          {
             other=lcc.template beta<2,0>(other);
          }
          while (lcc.template attribute<2>(lcc.template beta<2>(other))!=NULL);
          assert(lcc.template is_free<0>(lcc.template beta<2>(other)));
          lcc.basic_link_beta_1(it, lcc.template beta<2>(other));

          // For BGL halfedge graph, darts of border vertices must be border darts.
          lcc.template set_dart_of_attribute<0>(lcc.vertex_attribute(it), it);
          ++nb;
        }
      }

      return first_dart;
    }

  protected:

    // test if it exists a dart between vertex number i and vertex number j
    // If such a dart exists, returns it; otherwise returns null_handle.
    Dart_handle find_dart_between(size_type i, size_type j)
    {
      typename std::vector<Dart_handle>::reverse_iterator
        it(vertex_to_dart_map[i].rbegin());
      typename std::vector<Dart_handle>::reverse_iterator
        itend(vertex_to_dart_map[i].rend());

      Vertex_attribute_handle vh = vertex_map[j];
      for ( ; it!=itend; ++it )
      {
        if ( lcc.temp_vertex_attribute(*it)==vh )
          return (*it);
      }
      return lcc.null_handle;
    }

    // Add adart in the vertex number i (that means adart has vertex number
    // i as origin.
    void add_dart_in_vertex_to_dart_map( Dart_handle adart, size_type i )
    {
      CGAL_assertion( adart!=lcc.null_handle );
      vertex_to_dart_map[i].push_back(adart);
    }

  private:
    // The whole vertex of the lcc.
    std::vector<Vertex_attribute_handle> vertex_map;

    // A vector of vector of dart handle.
    // vertex_to_dart_map[i] contains all the darts that belong to vertex i.
    std::vector<std::vector<Dart_handle> > vertex_to_dart_map;

    LCC&        lcc;
    Dart_handle first_dart;
    Dart_handle prev_dart;
    size_type   first_vertex;
    size_type   prev_vertex;
    size_type   new_vertices;
    bool        begin_face;
    typename LCC::template Attribute_handle<2>::type face_attrib;
  };

  template < class LCC >
  void load_off_v2(LCC& alcc, std::istream& in)
  {
    File_header_OFF  m_file_header;
    File_scanner_OFF scanner( in, m_file_header.verbose());
    if ( ! in) return;
    m_file_header = scanner;  // Remember file header after return.

    Linear_cell_complex_incremental_builder_3_v2<LCC> B( alcc);
    B.begin_surface( scanner.size_of_vertices(),
                     scanner.size_of_facets(),
                     scanner.size_of_halfedges());

    typedef typename LCC::Point Point;

    // read in all vertices
    std::size_t  i;
    for ( i = 0; i < scanner.size_of_vertices(); i++) {
      Point p;
      file_scan_vertex( scanner, p);
      B.add_vertex( p);
      scanner.skip_to_next_vertex( i);
    }
    /* TODO rollback
       if ( ! in  || B.error()) {
       B.rollback();
       in.clear( std::ios::badbit);
       return;
       }
    */

    // read in all facets
    for ( i = 0; i < scanner.size_of_facets(); i++)
    {
      B.begin_facet();
      std::size_t no;
      scanner.scan_facet( no, i);
      /* TODO manage errors
         if( ! in || B.error() || no < 3) {
         if ( scanner.verbose()) {
         std::cerr << " " << std::endl;
         std::cerr << "Polyhedron_scan_OFF<Traits>::" << std::endl;
         std::cerr << "operator()(): input error: facet " << i
         << " has less than 3 vertices." << std::endl;
         }
         B.rollback();
         in.clear( std::ios::badbit);
         return;
         } */
      for ( std::size_t j = 0; j < no; j++) {
        std::size_t index;
        scanner.scan_facet_vertex_index( index, i);
        B.add_vertex_to_facet( index);
      }
      B.end_facet();
      scanner.skip_to_next_facet( i);
    }
    /* TODO manage errors
       if ( ! in  || B.error()) {
       B.rollback();
       in.clear( std::ios::badbit);
       return;
       }
       if ( B.check_unconnected_vertices()) {
       if ( ! B.remove_unconnected_vertices()) {
       if ( scanner.verbose()) {
       std::cerr << " " << std::endl;
       std::cerr << "Polyhedron_scan_OFF<Traits>::" << std::endl;
       std::cerr << "operator()(): input error: cannot "
       "succesfully remove isolated vertices."
       << std::endl;
       }
       B.rollback();
       in.clear( std::ios::badbit);
       return;
       }
       }*/
    B.end_surface();
  }

} //namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_INCREMENTAL_BUILDER_V2_H //
// EOF //
