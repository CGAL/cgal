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
#include <CGAL/IO/File_writer_wavefront.h>

namespace CGAL {

CGAL_INLINE_FUNCTION
void
File_writer_wavefront::
write_header( std::ostream& o,
              std::size_t   vertices,
              std::size_t   halfedges,
              std::size_t   facets){
    m_out    = &o;
    m_facets = facets;
    // Print header.
    out() << "# file written from a CGAL tool in Wavefront obj format\n";
    out() << "# " << vertices  << " vertices\n";
    out() << "# " << halfedges << " halfedges\n";
    out() << "# " << facets    << " facets\n\n";

    out() << "\n# " << vertices << " vertices\n";
    out() << "# ------------------------------------------\n\n";
}

} //namespace CGAL
// EOF //
