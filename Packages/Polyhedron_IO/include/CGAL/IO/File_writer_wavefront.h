// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_FILE_WRITER_WAVEFRONT_H
#define CGAL_IO_FILE_WRITER_WAVEFRONT_H 1
#ifndef CGAL_IO_BINARY_FILE_IO_H
#include <CGAL/IO/binary_file_io.h>
#endif // CGAL_IO_BINARY_FILE_IO_H
#ifndef CGAL_PROTECT_IOSTREAM
#include <iostream>
#define CGAL_PROTECT_IOSTREAM
#endif
#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif

CGAL_BEGIN_NAMESPACE

class File_writer_wavefront {
    std::ostream*  m_out;
    std::size_t    m_facets;
public:
    std::ostream& out() const { return *m_out; }
    void write_header( std::ostream& out,
                       std::size_t vertices,
                       std::size_t halfedges,
                       std::size_t facets);
    void write_footer() const {
        out() << "\n# End of Wavefront obj format #" << std::endl;
    }
    void write_vertex( const double& x, const double& y, const double& z) {
        out() << "v " << x << ' ' << y << ' ' << z << '\n';
    }
    void write_facet_header() {
        out() << "\n# " << m_facets << " facets\n";
        out() << "# ------------------------------------------\n\n";
    }
    void write_facet_begin( std::size_t)            { out() << "f "; }
    void write_facet_vertex_index( std::size_t idx) { out() << ' ' << idx+1; }
    void write_facet_end()                          { out() << '\n'; }
};

CGAL_END_NAMESPACE
#endif // CGAL_IO_FILE_WRITER_WAVEFRONT_H //
// EOF //
