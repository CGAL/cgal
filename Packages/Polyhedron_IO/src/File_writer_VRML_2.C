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
// file          : File_writer_VRML_2.C
// chapter       : $CGAL_Chapter: Support Library ... $
// package       : $CGAL_Package: Polyhedron_IO 2.11 (04 Feb 2000) $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Writer for polyhedral surfaces VRML_2.0 format (.wrl)
// ============================================================================

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#ifndef CGAL_PROTECT_IOSTREAM
#include <iostream>
#define CGAL_PROTECT_IOSTREAM
#endif
#ifndef CGAL_IO_FILE_WRITER_VRML_2_H
#include <CGAL/IO/File_writer_VRML_2.h>
#endif // CGAL_IO_FILE_WRITER_VRML_2_H

CGAL_BEGIN_NAMESPACE

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

CGAL_END_NAMESPACE
// EOF //
