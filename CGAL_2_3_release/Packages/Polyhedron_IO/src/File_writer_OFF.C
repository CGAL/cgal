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
// file          : File_writer_OFF.C
// chapter       : $CGAL_Chapter: Support Library ... $
// package       : $CGAL_Package: Polyhedron_IO 2.11 (04 Feb 2000) $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Writer for polyhedral surfaces in object file format (OFF)
// ============================================================================

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#ifndef CGAL_PROTECT_IOSTREAM
#include <iostream>
#define CGAL_PROTECT_IOSTREAM
#endif
#ifndef CGAL_IO_FILE_WRITER_OFF_H
#include <CGAL/IO/File_writer_OFF.h>
#endif // CGAL_IO_FILE_WRITER_OFF_H

CGAL_BEGIN_NAMESPACE

void
File_writer_OFF::
write_header( std::ostream& o,
              std::size_t   vertices,
              std::size_t   halfedges,
              std::size_t   facets,
              int           normals) {
    m_out = &o;
    m_header.set_vertices(  vertices);
    // Don't. This halfdges aren't trusted:
    // m_header.set_halfedges( halfedges);
    (void)halfedges;
    m_header.set_facets(  facets);
    m_header.set_normals( normals);
    // Print header.
    out() << m_header;
}
CGAL_END_NAMESPACE
// EOF //
