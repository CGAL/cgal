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
// file          : Polyhedron_geomview_ostream.h
// package       : $CGAL_Package: Polyhedron_IO 1.1 (01 Mar 1998) $
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Output stream operator for Polyhedrons into Geomview_stream.
// ============================================================================

#ifndef CGAL_IO_POLYHEDRON_GEOMVIEW_OSTREAM_H
#define CGAL_IO_POLYHEDRON_GEOMVIEW_OSTREAM_H 1

#ifndef CGAL_IO_GEOMVIEW_STREAM_H
#include <CGAL/IO/Geomview_stream.h>
#endif // CGAL_IO_GEOMVIEW_STREAM_H
#ifndef CGAL_IO_GENERIC_PRINT_POLYHEDRON_H
#include <CGAL/IO/generic_print_polyhedron.h>
#endif // CGAL_IO_GENERIC_PRINT_POLYHEDRON_H
#ifndef CGAL_POLYHEDRON_3_H
#include <CGAL/Polyhedron_3.h>
#endif // CGAL_POLYHEDRON_3

class CGAL_Polyhedron_writer_geomview {
    CGAL_Geomview_stream*  out;
    size_t   _facets;
public:
    CGAL_Polyhedron_writer_geomview( CGAL_Geomview_stream& geo) : out(&geo) {}
    void header( ostream&, size_t vertices, size_t, size_t facets) {
	// ignore ostream. Output goes to Geomview_stream.
	_facets = facets;
	// Print header.
	out->set_ascii_mode();
	*out << "(geometry polyhedron  {appearance {}{ ";
	out->set_binary_mode();
	*out << "OFF BINARY\n"  << int(vertices) << int(facets) << 0 ;
    }
    void footer() {
        *out << "}})";
	out->set_ascii_mode();
    }
    void write_vertex( const double& x, const double& y, const double& z) {
	*out << x << y << z;
    }
    void write_facet_header() {}
    void write_facet_begin( size_t no) {
	*out << int(no);
    }
    void write_facet_vertex_index( size_t index) {
	*out << int(index);
    }
    void write_facet_end() {
	double r = out->fcr(),
	       g = out->fcg(),
	       b = out->fcb();
	*out << 4 << r << g << b << 1.0;
    }
};


template < class Traits, class HDS >
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream &gv,
           const CGAL_Polyhedron_3<Traits,HDS> &P)
{
    CGAL_Polyhedron_writer_geomview  writer(gv);
    CGAL_generic_print_polyhedron( cerr, P, writer);
    return gv;
}


#endif // CGAL_IO_POLYHEDRON_GEOMVIEW_OSTREAM_H //
// EOF //

