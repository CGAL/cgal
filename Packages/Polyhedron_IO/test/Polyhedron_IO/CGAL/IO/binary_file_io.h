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
// file          : binary_file_io.h
// package       : $CGAL_Package: $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Binary read and write on streams for Integer32 and float
// ============================================================================

#ifndef CGAL_IO_BINARY_FILE_IO_H
#define CGAL_IO_BINARY_FILE_IO_H 1
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

// Forward declarations.
class ostream;

void inline
CGAL__Binary_write_integer32(ostream& out, CGAL_Integer32 i) {
    out.write( (char*)(&i), 4);
}
void inline
CGAL__Binary_write_float32(ostream& out, float f) {
    out.write( (char*)(&f), 4);
}

void inline
CGAL__Binary_read_integer32(istream& in, CGAL_Integer32& i) {
    in.read( (char*)(&i), 4);
}
void inline
CGAL__Binary_read_float32(istream& in, float& f) {
    in.read( (char*)(&f), 4);
}

void inline
CGAL__swap_to_big_endian( CGAL_UInteger32& u) {
    if( CGAL_LITTLE_ENDIAN) {
        u = ((u >> 24) | (u << 24) |
             ((u >> 8) & 0xff00) | ((u << 8) & 0xff0000));
    }
}

void inline
CGAL__swap_to_big_endian( CGAL_Integer32& i) {
    CGAL_UInteger32& u = (CGAL_UInteger32&)i;
    CGAL__swap_to_big_endian( u);
}

void inline
CGAL__swap_to_big_endian( float& f) {
    CGAL_UInteger32& u = (CGAL_UInteger32&)f;
    CGAL__swap_to_big_endian( u);
}

void inline
CGAL__Binary_write_big_endian_integer32(ostream& out, CGAL_Integer32 i) {
    CGAL__swap_to_big_endian( i);
    out.write( (char*)(&i), 4);
}
void inline
CGAL__Binary_write_big_endian_float32(ostream& out, float f) {
    CGAL__swap_to_big_endian( f);
    out.write( (char*)(&f), 4);
}

void inline
CGAL__Binary_read_big_endian_integer32(istream& in, CGAL_Integer32& i) {
    in.read( (char*)(&i), 4);
    CGAL__swap_to_big_endian( i);
}
void inline
CGAL__Binary_read_big_endian_float32(istream& in, float& f) {
    in.read( (char*)(&f), 4);
    CGAL__swap_to_big_endian( f);
}
#endif // CGAL_IO_BINARY_FILE_IO_H //
// EOF //
