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

#ifndef CGAL_IO_FILE_WRITER_INVENTOR_H
#define CGAL_IO_FILE_WRITER_INVENTOR_H

#include <CGAL/IO/OI/Inventor_ostream.h>

#include <iostream>

namespace CGAL {

class File_writer_inventor
{
  Inventor_ostream_base m_os;
  std::size_t m_facets;

public:
  File_writer_inventor() {}
  std::ostream& out() const { return m_os.os(); }

  void write_header(Inventor_ostream_base& o,
                    std::size_t vertices,
                    std::size_t halfedges,
                    std::size_t facets,
                    const bool /*colors*/ = false,
                    const bool /*normals*/ = false,
                    const bool /*textures*/ = false)
  {
    m_os = o;
    m_facets = facets;

    out() << "# " << vertices << " vertices\n";
    out() << "# " << halfedges << " halfedges\n";
    out() << "# " << facets << " facets\n\n";
    out() << "Separator {\n"
             "    Coordinate3 {\n"
             "        point   [" << std::endl;
  }

  void write_footer() const
  {
    out() << "        ] #coordIndex\n"
             "    } #IndexedFaceSet\n"
             "} #Separator" << std::endl;
  }

  void write_vertex( const double x, const double y, const double z)
  {
    out() << "            " << IO::oformat(x) << ' ' << IO::oformat(y) << ' ' << IO::oformat(z) << ',' <<'\n';
  }

  void write_vertex_normal(const double, const double, const double) { }
  void write_vertex_color(const double, const double, const double) { }
  void write_vertex_texture(const double, const double) { }

  void write_facet_header() const
  {
    out() << "        ] #point\n"
             "    } #Coordinate3\n"
             "    # " << m_facets << " facets\n"
             "    IndexedFaceSet {\n"
             "        coordIndex [\n";
  }

  void write_facet_begin( std::size_t) { out() << "            "; }
  void write_facet_vertex_index( std::size_t idx) { out() << idx << ','; }
  void write_face_color(const double, const double, const double) { }
  void write_facet_end() { out() << "-1,\n"; }
};

} // namespace CGAL

#endif // CGAL_IO_FILE_WRITER_INVENTOR_H
