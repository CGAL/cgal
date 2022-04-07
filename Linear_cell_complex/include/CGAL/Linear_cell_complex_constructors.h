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
#include <CGAL/Linear_cell_complex_incremental_builder.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/assertions.h>

#include <algorithm>
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

    // TODO FOR index do we need the Unique_hash_map ?
    // size_t i=0;
    // CGAL::Unique_hash_map< typename LCC::Vertex_attribute_const_handle,
    //     size_t, typename LCC::Hash_function > index;
    for (vit=alcc.vertex_attributes().begin(); vit!=vend; ++vit)
    {
      writer.write_vertex(::CGAL::to_double(vit->point().x()),
                          ::CGAL::to_double(vit->point().y()),
                          ::CGAL::to_double(vit->point().z()));
      // TODO for index ?? index[i++]=vit;
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
          CGAL_assertion(alcc.is_next_exist(cur));
          cur=alcc.next(cur);
        }
        while(cur!=itall);

        CGAL_assertion( n>=3 );
        writer.write_facet_begin(n);

        // Second we write the indices of vertices.
        do
        {
          // TODO case with index
         // TODO For index ? writer.write_facet_vertex_index(index[alcc.vertex_attribute(itf)]);// ?
          writer.write_facet_vertex_index(index[VCI(alcc.vertex_attribute(cur))]);
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
