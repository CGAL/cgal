// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_FILE_WRITER_VRML_2_H
#define CGAL_IO_FILE_WRITER_VRML_2_H 1

#include <CGAL/basic.h>
#include <iostream>
#include <cstddef>

namespace CGAL {

class CGAL_EXPORT File_writer_VRML_2 {
    std::ostream*      m_out;
    std::size_t        m_facets;
public:
    File_writer_VRML_2() {}
    std::ostream& out() const { return *m_out; }
    void write_header( std::ostream& o,
                       std::size_t vertices,
                       std::size_t halfedges,
                       std::size_t facets);
    void write_footer() const;
    void write_vertex( const double& x, const double& y, const double& z) {
        out() << "                                "
              << x << ' ' << y << ' ' << z << ',' << '\n';
    }
    void write_facet_header() const;
    void write_facet_begin( std::size_t) {
        out() << "                            ";
    }
    void write_facet_vertex_index( std::size_t idx) { out() << idx << ',';}
    void write_facet_end() { out() << "-1,\n"; }
};

} //namespace CGAL

#ifdef CGAL_HEADER_ONLY
#include <CGAL/IO/File_writer_VRML_2_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_IO_FILE_WRITER_VRML_2_H //
// EOF //
