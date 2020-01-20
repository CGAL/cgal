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
#define CGAL_IO_FILE_WRITER_INVENTOR_H 1

#include <CGAL/basic.h>

#include <iostream>
#include <cstddef>

namespace CGAL {

class CGAL_EXPORT File_writer_inventor
{
  std::ostream*      m_out;
  std::size_t        m_facets;

public:
  File_writer_inventor() {}
  std::ostream& out() const { return *m_out; }

  void write_header(std::ostream& o,
                    std::size_t vertices,
                    std::size_t halfedges,
                    std::size_t facets)
  {
    m_out = &o;
    m_facets = facets;

    out() << "# " << vertices  << " vertices\n";
    out() << "# " << halfedges << " halfedges\n";
    out() << "# " << facets    << " facets\n\n";
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

  void write_vertex( const double& x, const double& y, const double& z) {
    out() << "            " << x << ' ' << y << ' ' << z << ',' <<'\n';
  }

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
  void write_facet_end() { out() << "-1,\n"; }
};

} //namespace CGAL

#endif // CGAL_IO_FILE_WRITER_INVENTOR_H
