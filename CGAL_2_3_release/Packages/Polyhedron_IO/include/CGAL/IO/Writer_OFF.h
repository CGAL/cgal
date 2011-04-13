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
// chapter       : $CGAL_Chapter: Support Library ... $
// package       : $CGAL_Package: Polyhedron_IO 2.11 (04 Feb 2000) $
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

CGAL_BEGIN_NAMESPACE

template < class Pt >
class Writer_OFF : public Generic_writer<File_writer_OFF,Pt> {
public:
    Writer_OFF( std::ostream& out,
                std::size_t   v,
                std::size_t   h,
                std::size_t   f,
                bool          verbose = false)
    : Generic_writer<File_writer_OFF,Pt>(
        File_writer_OFF( verbose), out, v, h, f
    ) {}
    Writer_OFF( std::ostream& out,
                std::size_t   v,
                std::size_t   h,
                std::size_t   f,
                const File_header_OFF& header)
    : Generic_writer<File_writer_OFF,Pt>(
        File_writer_OFF( header), out, v, h, f
    ) {}
};

CGAL_END_NAMESPACE
#endif // CGAL_IO_WRITER_OFF_H //
// EOF //
