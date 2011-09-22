// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_FILE_WRITER_GOCAD_H
#define CGAL_IO_FILE_WRITER_GOCAD_H 1

#include <CGAL/IO/binary_file_io.h>
#include <CGAL/IO/File_header_gocad.h>
#include <iostream>
#include <cstddef>

namespace CGAL {

class File_writer_gocad {
    int m_vertex_index;
    std::ostream*           m_out;
    File_header_gocad         m_header;
public:
    File_writer_gocad(std::string fname, std::string color) : m_vertex_index(0), m_header(fname, color) {}
    File_writer_gocad( const File_header_gocad& h) : m_vertex_index(0), m_header( h) {}

    std::ostream&           out()          { return *m_out;   }
    File_header_gocad&        header()       { return m_header; }
    const File_header_gocad&  header() const { return m_header; }

    void write_header( std::ostream& o,
                       std::size_t   vertices,
                       std::size_t   halfedges,
                       std::size_t   facets,
                       bool          normals = false)
  {
    m_out = &o;
    out() << header();
  }


    void write_footer() 
  {}

    void write_vertex( const double& x, const double& y, const double& z) {
      out() << "VRTX " << m_vertex_index++ << ' ' << x << ' ' << y << ' ' << z << '\n';
    }

    void write_facet_header() {
                out() << '\n';            
    }

    void write_facet_begin( std::size_t no) {
      CGAL_assertion(no == 3);
      out() << "TRGL " << ' ';
    }
    void write_facet_vertex_index( std::size_t index) {
            out() << ' ' << index;
    }
    void write_facet_end() {
            out() << '\n';
    }
};

} //namespace CGAL
#endif // CGAL_IO_FILE_WRITER_GOCAD_H //
// EOF //
