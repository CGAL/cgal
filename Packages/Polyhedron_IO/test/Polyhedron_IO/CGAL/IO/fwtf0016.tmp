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
// file          : Writer_OFF.h
// package       : $CGAL_Package: $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// STL compliant interface to write OFF.
// ============================================================================

#ifndef CGAL_IO_WRITER_OFF_H
#define CGAL_IO_WRITER_OFF_H 1
#ifndef CGAL_IO_FILE_WRITER_OFF_H
#include <CGAL/IO/File_writer_OFF.h>
#endif // CGAL_IO_FILE_WRITER_OFF_H
#ifndef CGAL_IO_GENERIC_WRITER_H
#include <CGAL/IO/Generic_writer.h>
#endif // CGAL_IO_GENERIC_WRITER_H

template < class Pt >
class CGAL_Writer_OFF : public CGAL_Generic_writer<CGAL_File_writer_OFF,Pt> {
public:
    CGAL_Writer_OFF( ostream& out, size_t v, size_t h, size_t f,
                    bool binary = false, bool nocomments = false,
                    bool skel = false)
    : CGAL_Generic_writer<CGAL_File_writer_OFF,Pt>(
        CGAL_File_writer_OFF( binary, nocomments, skel), out, v, h, f
    ) {}
};
#endif // CGAL_IO_WRITER_OFF_H //
// EOF //
