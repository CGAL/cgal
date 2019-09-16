// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
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
