// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_POLYHEDRON_GEOMVIEW_OSTREAM_H
#define CGAL_IO_POLYHEDRON_GEOMVIEW_OSTREAM_H

#include <CGAL/license/Polyhedron.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/assertions.h>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/boost/graph/IO/Generic_facegraph_printer.h>

namespace CGAL {

class Polyhedron_writer_geomview
{
  Geomview_stream* out;

public:
  void write_header(Geomview_stream& os,
                    std::size_t vertices,
                    std::size_t,
                    std::size_t facets,
                    bool = false,
                    bool = false,
                    bool = false)
  {
    out = &os;

    // Print header.
    out->set_ascii_mode();
    *out << "(geometry " << out->get_new_id("polyhedron")
         << " {appearance {}{ ";
    out->set_binary_mode();
    *out << "OFF BINARY\n"  << int(vertices) << int(facets) << 0 ;
  }

  void write_footer()
  {
    *out << "}})";
    out->set_ascii_mode();
  }

  void write_vertex( const double& x, const double& y, const double& z) { *out << x << y << z; }
  void write_vertex_normal(const double, const double, const double) { }
  void write_vertex_color(const double, const double, const double) { }
  void write_vertex_texture(const double, const double) { }

  void write_facet_header() {}
  void write_facet_begin(std::size_t no) { *out << int(no); }
  void write_facet_vertex_index( std::size_t index) { *out << int(index); }
  void write_face_color(const double, const double, const double) { }
  void write_facet_end()
  {
    double r = out->fcr(),
           g = out->fcg(),
           b = out->fcb();
    *out << 4 << r << g << b << 1.0;
  }
};

template <class Traits,
          class Items,
          template < class T, class I, class A>
          class HDS, class Alloc>
Geomview_stream& operator<<(Geomview_stream &gv,
                            const Polyhedron_3<Traits, Items, HDS, Alloc>& P)
{
  IO::internal::Generic_facegraph_printer<Geomview_stream,
                                          Polyhedron_3<Traits, Items, HDS, Alloc>,
                                          Polyhedron_writer_geomview> printer(gv);
  printer(P);

  return gv;
}

} //namespace CGAL

#endif // CGAL_IO_POLYHEDRON_GEOMVIEW_OSTREAM_H //
// EOF //
