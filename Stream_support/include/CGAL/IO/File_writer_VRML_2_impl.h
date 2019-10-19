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

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/IO/File_writer_VRML_2.h>

namespace CGAL {

CGAL_INLINE_FUNCTION
void
File_writer_VRML_2::
write_header( std::ostream& o,
              std::size_t   vertices,
              std::size_t   halfedges,
              std::size_t   facets) {
    m_out    = &o;
    m_facets = facets;

    out() << "        #-- Begin of Polyhedron_3\n";
    out() << "        # " << vertices  << " vertices\n";
    out() << "        # " << halfedges << " halfedges\n";
    out() << "        # " << facets    << " facets\n";
    out() << "        Group {\n"
             "            children [\n"
             "                Shape {\n"
             "                    appearance Appearance { material "
                                               "USE Material }\n"
             "                    geometry IndexedFaceSet {\n"
             "                        convex FALSE\n"
             "                        solid  FALSE\n"
             "                        coord  Coordinate {\n"
             "                            point [" << std::endl;
}

void
File_writer_VRML_2::
write_facet_header() const {
    out() << "                            ] #point\n"
             "                        } #coord Coordinate\n"
             "                        coordIndex  [" << std::endl;
}

void
File_writer_VRML_2::
write_footer() const {
    out() << "                        ] #coordIndex\n"
             "                    } #geometry\n"
             "                } #Shape\n"
             "            ] #children\n"
             "        } #Group" << std::endl;
}

} //namespace CGAL
// EOF //
