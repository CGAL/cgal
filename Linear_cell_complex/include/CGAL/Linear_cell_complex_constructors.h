// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_LINEAR_CELL_COMPLEX_CONSTRUCTORS_H
#define CGAL_LINEAR_CELL_COMPLEX_CONSTRUCTORS_H 1

#include <CGAL/IO/File_header_OFF.h>
#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/IO/File_writer_OFF.h>
#include <CGAL/Linear_cell_complex_incremental_builder.h>

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <list>

namespace CGAL {

  /** @file Linear_cell_complex_constructors.h
   * Some construction operations for a linear cell complex from other
   * CGAL data structures.
   */

  /** Import an embedded plane graph read into a flux into a
   *  linear cell complex.
   * @param alcc the linear cell complex where the graph will be imported.
   * @param ais the istream where read the graph.
   * @return A dart created during the convertion.
   */
  template< class LCC >
  typename LCC::Dart_handle import_from_plane_graph(LCC& alcc,
                                                    std::istream& ais)
  {
    CGAL_static_assertion( LCC::dimension>=2 && LCC::ambient_dimension==2 );

    typedef typename LCC::Dart_handle Dart_handle;
    typedef typename LCC::Traits::Direction_2 Direction;
    typedef typename std::list<Dart_handle>::iterator List_iterator;
    typedef typename std::map<Direction, Dart_handle>::iterator LCC_iterator;

    // Arrays of vertices
    std::vector< typename LCC::Vertex_attribute_handle > initVertices;
    std::vector< std::list<Dart_handle> > testVertices;

    std::string txt;
    typename LCC::FT x, y;
    Dart_handle d1=alcc.null_handle;
    unsigned int v1, v2;

    unsigned int nbSommets = 0;
    unsigned int nbAretes = 0;

    ais >> nbSommets >> nbAretes;
    while (nbSommets > 0)
    {
      if (!ais.good())
      {
        std::cout << "Problem: file does not contain enough vertices."
                  << std::endl;
        return alcc.null_handle;
      }

      ais >> iformat(x) >> iformat(y);
      initVertices.push_back(alcc.create_vertex_attribute
                             (typename LCC::Point(x, y)));
      testVertices.push_back(std::list<Dart_handle>());
      --nbSommets;
    }

    while (nbAretes>0)
    {
      if (!ais.good())
      {
        std::cout << "Problem: file does not contain enough edges."
                  << std::endl;
        return alcc.null_handle;
      }

      // We read an egde (given by the number of its two vertices).
      ais >> v1 >> v2;
      --nbAretes;

      CGAL_assertion(v1 < initVertices.size());
      CGAL_assertion(v2 < initVertices.size());

      d1 = alcc.make_segment(initVertices[v1], initVertices[v2], true);

      testVertices[v1].push_back(d1);
      testVertices[v2].push_back(alcc.template opposite<2>(d1));
    }

    // LCC associating directions and darts.
    std::map<Direction, Dart_handle> tabDart;
    List_iterator it;
    LCC_iterator  it2;

    Dart_handle first = alcc.null_handle;
    Dart_handle prec = alcc.null_handle;
    typename LCC::Point sommet1, sommet2;

    for (unsigned int i=0; i<initVertices.size(); ++i)
    {
      it = testVertices[i].begin();
      if (it != testVertices[i].end()) // Si la liste n'est pas vide.
      {
        // 1. We insert all the darts and sort them depending on the direction
        tabDart.clear();

        sommet1 = alcc.point(*it);
        sommet2 = alcc.point(alcc.other_extremity(*it));

        tabDart.insert(std::pair<Direction, Dart_handle>
                       (typename LCC::Traits::Construct_direction_2()
                        (typename LCC::Traits::Construct_vector()
                         (sommet1,sommet2)), *it));

        ++it;
        while (it!=testVertices[i].end())
        {
          sommet2 = alcc.point(alcc.other_extremity(*it));
          tabDart.insert(std::pair<Direction, Dart_handle>
                         (typename LCC::Traits::Construct_direction_2()
                          (typename LCC::Traits::Construct_vector()
                           (sommet1,sommet2)), *it));
          ++it;
        }

        // 2. We run through the array of darts and 1 links darts.
        it2 = tabDart.begin();
        first = it2->second;
        prec = first;
        ++it2;

        while (it2!=tabDart.end())
        {
          alcc.set_next(alcc.template opposite<2>(it2->second), prec);
          prec = it2->second;
          ++it2;
        }
        alcc.set_next(alcc.template opposite<2>(first), prec);
      }
    }

    // We return a dart from the imported object.
    return first;
  }

  template < class LCC >
  typename LCC::Dart_handle
  import_from_plane_graph(LCC& alcc, const char* filename)
  {
    std::ifstream input(filename);
    if (!input.is_open()) return alcc.null_handle;
    return import_from_plane_graph(alcc, input);
  }

  template < class LCC >
  bool load_off(LCC& alcc, std::istream& in)
  {
    File_header_OFF  m_file_header;
    File_scanner_OFF scanner( in, m_file_header.verbose());
    if (!in) return false;
    m_file_header = scanner;  // Remember file header after return.

    Linear_cell_complex_incremental_builder_3<LCC> B(alcc);
    B.begin_surface(scanner.size_of_vertices(),
                    scanner.size_of_facets(),
                    scanner.size_of_halfedges());

    typedef typename LCC::Point Point;

    // read in all vertices
    std::size_t  i;
    for (i = 0; i < scanner.size_of_vertices(); i++)
    {
      Point p;
      file_scan_vertex(scanner, p);
      B.add_vertex(p);
      scanner.skip_to_next_vertex(i);
    }
    /* TODO rollback
       if ( ! in  || B.error()) {
       B.rollback();
       in.clear( std::ios::badbit);
       return;
       }
    */

    // read in all facets
    for (i=0; i<scanner.size_of_facets(); i++)
    {
      B.begin_facet();
      std::size_t no;
      scanner.scan_facet(no, i);
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
      for (std::size_t j=0; j<no; j++)
      {
        std::size_t index;
        scanner.scan_facet_vertex_index(index, i);
        B.add_vertex_to_facet(index);
      }
      B.end_facet();
      scanner.skip_to_next_facet(i);
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
       "successfully remove isolated vertices."
       << std::endl;
       }
       B.rollback();
       in.clear( std::ios::badbit);
       return;
       }
       }*/
    B.end_surface();

    return true;
  }

  template < class LCC >
  bool load_off(LCC& alcc, const char* filename)
  {
    std::ifstream input(filename);
    if (!input.is_open())
    { return false; }

    return load_off(alcc, input);
  }

  /** Export the alcc in off file format. If dimension>2, export all faces but only once.
   *  @pre all faces are closed (i.e. form by closed cycles of edges)
   */
  template < class LCC >
  bool write_off(LCC& alcc, std::ostream& out)
  {
    if (!alcc.are_all_faces_closed())
    {
      std::cerr<<"Impossible to write in off a map having open faces."<<std::endl;
      return false;
    }

    File_header_OFF header(false);
    header.set_binary(is_binary(out));
    header.set_no_comments(!is_pretty(out));
    File_writer_OFF writer( header);
    writer.header().set_polyhedral_surface(true);
    writer.header().set_halfedges(alcc.number_of_darts());

    // Print header.
    writer.write_header(out,
                        alcc.number_of_vertex_attributes(),
                        alcc.number_of_halfedges(),
                        alcc.template one_dart_per_cell<2>().size() );

    typedef typename LCC::Vertex_attribute_range::iterator VCI;
    VCI vit, vend = alcc.vertex_attributes().end();
    for (vit=alcc.vertex_attributes().begin(); vit!=vend; ++vit)
    {
      writer.write_vertex(::CGAL::to_double(vit->point().x()),
                          ::CGAL::to_double(vit->point().y()),
                          ::CGAL::to_double(vit->point().z()));
    }

    typedef Inverse_index< VCI > Index;
    Index index( alcc.vertex_attributes().begin(),
                 alcc.vertex_attributes().end());
    writer.write_facet_header();

    typename LCC::size_type m = alcc.get_new_mark();

    for (typename LCC::Dart_range::iterator itall = alcc.darts().begin(),
           itallend = alcc.darts().end(); itall!=itallend; ++itall)
    {
      if (!alcc.is_marked(itall, m))
      {
        std::size_t n = 0;
        typename LCC::Dart_handle cur=itall;
        do
        {
          ++n;
          assert(alcc.is_next_exist(cur));
          cur=alcc.next(cur);
        }
        while(cur!=itall);

        CGAL_assertion( n>=3 );
        writer.write_facet_begin(n);

        // Second we write the indices of vertices.
        do
        {
          // TODO case with index
          writer.write_facet_vertex_index(index[VCI(alcc.vertex_attribute(cur))]);
          alcc.mark(cur, m);
          alcc.mark(alcc.other_orientation(cur), m); // for GMap only, for CMap
          assert(alcc.is_next_exist(cur));           // marks the same dart twice
          cur=alcc.next(cur);
        }
        while(cur!=itall);

        writer.write_facet_end();
      }
    }
    writer.write_footer();
    alcc.free_mark(m);
    return true;
  }

  template < class LCC >
  bool write_off(LCC& alcc, const char* filename)
  {
    std::ofstream output(filename);
    if (!output.is_open())
    { return false; }
    
    return write_off(alcc, output);
  }

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_CONSTRUCTORS_H //
// EOF //
