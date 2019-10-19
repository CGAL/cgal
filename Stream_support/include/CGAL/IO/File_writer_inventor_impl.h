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

#include <CGAL/IO/File_writer_inventor.h>

namespace CGAL {

CGAL_INLINE_FUNCTION
void
File_writer_inventor::
write_header( std::ostream& o,
              std::size_t   vertices,
              std::size_t   halfedges,
              std::size_t   facets){
    m_out    = &o;
    m_facets = facets;
    out() << "# " << vertices  << " vertices\n";
    out() << "# " << halfedges << " halfedges\n";
    out() << "# " << facets    << " facets\n\n";
    out() << "Separator {\n"
             "    Coordinate3 {\n"
             "        point   [" << std::endl;
}

CGAL_INLINE_FUNCTION
void
File_writer_inventor::
write_facet_header() const {
    out() << "        ] #point\n"
             "    } #Coordinate3\n"
             "    # " << m_facets << " facets\n"
             "    IndexedFaceSet {\n"
             "        coordIndex [\n";
}

CGAL_INLINE_FUNCTION
void
File_writer_inventor::
write_footer() const {
    out() << "        ] #coordIndex\n"
             "    } #IndexedFaceSet\n"
             "} #Separator" << std::endl;
}

} //namespace CGAL
// EOF //
