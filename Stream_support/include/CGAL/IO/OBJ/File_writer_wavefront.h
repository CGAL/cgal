// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_OBJ_FILE_WRITER_WAVEFRONT_H
#define CGAL_IO_OBJ_FILE_WRITER_WAVEFRONT_H

#include <CGAL/IO/io.h>

#include <fstream>
#include <iostream>

namespace CGAL {

class File_writer_wavefront
{
  std::ostream* m_os;
  std::size_t m_facets;

public:
  std::ostream& out() const { return *m_os; }

  void write_header(std::ostream& o,
                    std::size_t vertices,
                    std::size_t halfedges,
                    std::size_t facets,
                    const bool /*colors*/ = false,
                    const bool /*normals*/ = false,
                    const bool /*textures*/ = false)
  {
    m_os = &o;
    m_facets = facets;

    // Print header.
    out() << "# file written from a CGAL tool in Wavefront obj format\n";
    out() << "# " << vertices << " vertices\n";
    out() << "# " << halfedges << " halfedges\n";
    out() << "# " << facets << " facets\n\n";

    out() << "\n# " << vertices << " vertices\n";
    out() << "# ------------------------------------------\n\n";
  }

  void write_footer() const { out() << "\n# End of Wavefront obj format #" << std::endl; }

  void write_vertex(const double x, const double y, const double z) {
    out() << "v " << IO::oformat(x) << ' ' << IO::oformat(y) << ' ' << IO::oformat(z) << '\n';
  }

  void write_vertex_normal(const double x, const double y, const double z) {
    out() << "vn " << IO::oformat(x) << ' ' << IO::oformat(y) << ' ' << IO::oformat(z) << '\n';
  }

  void write_vertex_color(const double, const double, const double) { }
  void write_vertex_texture(const double, const double) { }

  void write_facet_header()
  {
    out() << "\n# " << m_facets << " facets\n";
    out() << "# ------------------------------------------\n\n";
  }
  void write_facet_begin(std::size_t) { out() << "f "; }
  void write_facet_vertex_index(std::size_t idx) { out() << ' ' << idx+1; }
  void write_face_color(const double, const double, const double) { }
  void write_facet_end() { out() << '\n'; }
};

} // namespace CGAL

#endif // CGAL_IO_OBJ_FILE_WRITER_WAVEFRONT_H
