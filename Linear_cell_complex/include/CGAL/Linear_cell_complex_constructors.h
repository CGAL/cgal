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

#include <CGAL/IO/OFF.h>
#include <CGAL/Linear_cell_complex_incremental_builder_3.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/assertions.h>

#include <algorithm>
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


  /**
   * Imports a plane-embedded graph from a list of points and edges represented as pairs of vertex indices
   */
  template< class LCC >
  typename LCC::Dart_descriptor import_from_plane_graph(LCC& alcc,
                                                   const std::vector<typename LCC::Point>& vertices,
                                                   const std::vector<size_t>& edge_indices)
  {
    typedef typename LCC::Traits::Construct_direction_2 Construct_direction_2;
    typedef typename LCC::Traits::Construct_vector Construct_vector;
    typedef typename LCC::Dart_descriptor Dart_descriptor;
    typedef typename LCC::Traits::Direction_2 Direction;
    typedef typename std::map<Direction, Dart_descriptor>::iterator LCC_iterator;
    typedef typename std::list<Dart_descriptor>::iterator List_iterator;
    typedef typename LCC::Point Point;

    CGAL_static_assertion( LCC::dimension>=2 && LCC::ambient_dimension==2 );
    CGAL_assertion(edge_indices.size() % 2 == 0);

    std::vector< typename LCC::Vertex_attribute_descriptor > initVertices;
    initVertices.reserve(vertices.size());
    std::transform(vertices.begin(),
                   vertices.end(),
                   std::back_inserter(initVertices),
                   [&](const Point& point) {
      return alcc.create_vertex_attribute(point);
    });

    std::vector< std::list<Dart_descriptor> > testVertices{vertices.size(), std::list<Dart_descriptor>()};

    Dart_descriptor d1 = alcc.null_descriptor;
    for (std::size_t i = 0; (i + 1) < edge_indices.size(); i += 2) {
      const auto& v1 = edge_indices[i];
      const auto& v2 = edge_indices[i + 1];
      CGAL_assertion(v1 < initVertices.size());
      CGAL_assertion(v2 < initVertices.size());

      d1 = alcc.make_segment(initVertices[v1], initVertices[v2], true);

      testVertices[v1].push_back(d1);
      testVertices[v2].push_back(alcc.template opposite<2>(d1));
    }

    // LCC associating directions and darts.
    std::map<Direction, Dart_descriptor> tabDart;
    List_iterator it;
    LCC_iterator  it2;

    Dart_descriptor first = alcc.null_descriptor;
    Dart_descriptor prec = alcc.null_descriptor;

    for (unsigned int i=0; i<initVertices.size(); ++i)
    {
      it = testVertices[i].begin();
      if (it != testVertices[i].end()) // Si la liste n'est pas vide.
      {
        // 1. We insert all the darts and sort them depending on the direction
        tabDart.clear();

        Point vertex1 = alcc.point(*it);
        Point vertex2 = alcc.point(alcc.other_extremity(*it));

        tabDart.insert(std::pair<Direction, Dart_descriptor>
                       (Construct_direction_2()
                        (Construct_vector()
                        (vertex1,vertex2)), *it));

        ++it;
        while (it!=testVertices[i].end())
        {
          vertex2 = alcc.point(alcc.other_extremity(*it));
          tabDart.insert(std::pair<Direction, Dart_descriptor>
                         (Construct_direction_2()
                          (Construct_vector()
                           (vertex1,vertex2)), *it));
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

  /**
   * Imports a plane-embedded graph from a file into a LinearCellComplex.
   *
   * @param alcc the linear cell complex where the graph will be imported.
   * @param ais the istream where read the graph.
   * @return A dart created during the conversion.
   */
  template< class LCC >
  typename LCC::Dart_descriptor import_from_plane_graph(LCC& alcc,
                                                    std::istream& ais)
  {
    using FT = typename LCC::FT;
    using Point = typename LCC::Point;

    std::vector<Point> vertices;
    unsigned int numVertices = 0;
    unsigned int numEdges = 0;
    ais >> numVertices >> numEdges;
    while (numVertices > 0)
    {
      if (!ais.good())
      {
        std::cout << "Problem: file does not contain enough vertices."
                  << std::endl;
        return alcc.null_descriptor;
      }

      FT x, y;
      ais >> IO::iformat(x) >> IO::iformat(y);
      vertices.push_back(Point{x, y});
      --numVertices;
    }

   std::vector<size_t> edge_indices;
    while (numEdges>0)
    {
      if (!ais.good())
      {
        std::cout << "Problem: file does not contain enough edges."
                  << std::endl;
        return alcc.null_descriptor;
      }

      // We read an edge (given by the number of its two vertices).
      unsigned int v1, v2;
      ais >> v1 >> v2;
      --numEdges;

      CGAL_assertion(v1 < vertices.size());
      CGAL_assertion(v2 < vertices.size());

      edge_indices.push_back(v1);
      edge_indices.push_back(v2);
    }

    return import_from_plane_graph(alcc, vertices, edge_indices);
  }

  template < class LCC >
  typename LCC::Dart_descriptor
  import_from_plane_graph(LCC& alcc, const char* filename)
  {
    std::ifstream input(filename);
    if (!input.is_open()) return alcc.null_descriptor;
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
    B.begin_surface();

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
      std::size_t no=0;
      scanner.scan_facet(no, i);
      /* TODO manage errors
         if( ! in || B.error() || no < 3) {
         if ( scanner.verbose()) {
         std::cerr << " " << std::endl;
         std::cerr << "Polyhedron_scan_OFF<Traits>::" << std::endl;
         std::cerr << "operator()(): input error: facet " << i
         << " has fewer than 3 vertices." << std::endl;
         }
         B.rollback();
         in.clear( std::ios::badbit);
         return;
         } */
      for (std::size_t j=0; j<no; j++)
      {
        std::size_t index=0;
        scanner.scan_facet_vertex_index(index, j+1, i);
        if(! in){
          return false;
        }
        B.add_vertex_to_facet(static_cast<typename LCC::size_type>(index));
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
    header.set_binary(IO::is_binary(out));
    header.set_no_comments(!IO::is_pretty(out));
    File_writer_OFF writer( header);
    writer.header().set_polyhedral_surface(true);
    writer.header().set_halfedges(alcc.number_of_darts());

    // Print header.
    writer.write_header(out,
                        alcc.number_of_vertex_attributes(),
                        alcc.number_of_halfedges(),
                        alcc.template one_dart_per_cell<2>().size());

    typedef typename LCC::Vertex_attribute_range::const_iterator VCI;
    VCI vit, vend = alcc.vertex_attributes().end();

    // TODO FOR index we do not need the Unique_hash_map.
    size_t i=0;
    CGAL::Unique_hash_map< typename LCC::Vertex_attribute_const_descriptor,
        size_t, typename LCC::Hash_function > index;
    for (vit=alcc.vertex_attributes().begin(); vit!=vend; ++vit)
    {
      writer.write_vertex(::CGAL::to_double(vit->point().x()),
                          ::CGAL::to_double(vit->point().y()),
                          ::CGAL::to_double(vit->point().z()));
      index[vit]=i++; // TODO for index
    }

    writer.write_facet_header();
    typename LCC::size_type m = alcc.get_new_mark();

    for (typename LCC::Dart_range::iterator itall = alcc.darts().begin(),
           itallend = alcc.darts().end(); itall!=itallend; ++itall)
    {
      if (!alcc.is_marked(itall, m))
      {
        std::size_t n = 0;
        typename LCC::Dart_descriptor cur=itall;
        do
        {
          ++n;
          CGAL_assertion(alcc.is_next_exist(cur));
          cur=alcc.next(cur);
        }
        while(cur!=itall);

        CGAL_assertion( n>=3 );
        writer.write_facet_begin(n);

        // Second we write the indices of vertices.
        do
        {
          writer.write_facet_vertex_index(index[alcc.vertex_attribute(cur)]); // TODO for index
          alcc.mark(cur, m);
          alcc.mark(alcc.other_orientation(cur), m); // for GMap only, for CMap
          CGAL_assertion(alcc.is_next_exist(cur));           // marks the same dart twice
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
