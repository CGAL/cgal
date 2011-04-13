// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : File_writer_inventor.C
// chapter       : $CGAL_Chapter: Support Library ... $
// package       : $CGAL_Package: Polyhedron_IO 2.11 (04 Feb 2000) $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Writer for polyhedral surfaces OpenInventor format (.iv)
// ============================================================================

#ifndef CGAL_IO_FILE_WRITER_INVENTOR_H
#include <CGAL/IO/File_writer_inventor.h>
#endif // CGAL_IO_FILE_WRITER_INVENTOR_H

CGAL_BEGIN_NAMESPACE

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

void
File_writer_inventor::
write_facet_header() const {
    out() << "        ] #point\n"
             "    } #Coordinate3\n"
             "    # " << m_facets << " facets\n"
             "    IndexedFaceSet {\n"
             "        coordIndex [\n";
}

void
File_writer_inventor::
write_footer() const {
    out() << "        ] #coordIndex\n"
             "    } #IndexedFaceSet\n"
             "} #Separator" << std::endl;
}

CGAL_END_NAMESPACE
// EOF //
