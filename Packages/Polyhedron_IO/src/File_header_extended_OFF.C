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
// file          : src/File_header_extended_OFF.C
// package       : Polyhedron_IO 2.11 (04 Feb 2000)
// chapter       : Support Library
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>
// maintainer    :
// coordinator   : INRIA, Sophia Antipolis
//
// Extended file header information of an object file format (OFF) file
// ============================================================================

#include <CGAL/IO/File_header_extended_OFF.h>
#include <CGAL/basic.h>

#include <cstdlib>
#include <cctype>
#include <cstring>
#include <iostream>
#include <sstream>
#include <algorithm>


CGAL_BEGIN_NAMESPACE

bool File_header_extended_OFF::
is_POL()  const {
    return is_OFF() && polyhedral_surface();
}

bool File_header_extended_OFF::
is_CBP()  const {
    return is_POL() && triangulated() && non_empty_facets() &&
        normalized_to_sphere() && radius() <= 1.0;
}

bool File_header_extended_OFF::
is_TRN()  const { return is_CBP() && terrain(); }

int  File_header_extended_OFF::
is_CBPn() const {
    if ( is_POL() && triangulated() && non_empty_facets() &&
         normalized_to_sphere() && rounded() &&
         (radius() <= ( 1l << rounded_bits())))
        return rounded_bits();
    else
        return 0;
}

int  File_header_extended_OFF::
is_TRNn() const { return ( terrain() ? is_CBPn() : 0); }


// The proper file suffix with respect to file format.
std::string File_header_extended_OFF::
suffix() const {
    if ( is_TRNn()) {
        std::ostringstream out;
        out << "trn" << m_rounded_bits << '\0';
        return out.str();
    }
    if ( is_TRN())
        return std::string("trn");
    if ( is_CBPn()) {
        std::ostringstream out;
        out << "cbp" << m_rounded_bits << '\0';
        return out.str();
    }
    if ( is_CBP())
        return std::string("cbp");
    if ( is_POL())
        return std::string("pol");
    return std::string("off");
}

// The proper format name.
std::string File_header_extended_OFF::
format_name() const {
    if ( is_TRNn()) {
        std::ostringstream out;
        out << "TRN" << m_rounded_bits << '\0';
        return out.str();
    }
    if ( is_TRN())
        return std::string("TRN");
    if ( is_CBPn()) {
        std::ostringstream out;
        out << "CBP" << m_rounded_bits << '\0';
        return out.str();
    }
    if ( is_CBP())
        return std::string("CBP");
    if ( is_POL())
        return std::string("POL");
    return std::string("OFF");
}

File_header_extended_OFF& File_header_extended_OFF::
operator+=( const File_header_extended_OFF& header) {
    m_verbose              = m_verbose || header.m_verbose;
    m_polyhedral_surface   = m_polyhedral_surface &&
                             header.m_polyhedral_surface;
    m_halfedges           += header.m_halfedges;
    m_triangulated         = m_triangulated && header.m_triangulated;
    m_non_empty_facets     = m_non_empty_facets &&
                             header.m_non_empty_facets;
    m_terrain              = m_terrain && header.m_terrain;
    m_normalized_to_sphere = m_normalized_to_sphere &&
                             header.m_normalized_to_sphere;
    m_radius               = std::max(m_radius, header.m_radius);
    m_rounded              = m_rounded && header.m_rounded;
    m_rounded_bits         = std::max( m_rounded_bits,
                                       header.m_rounded_bits);
    m_off_header           = m_off_header && header.m_off_header;
    return *this;
}

#define OUT(item) out << "# " #item " " << h.item() << '\n'
#define OUTBOOL(item) out << "# " #item " " << (h.item() ? '1':'0') << '\n'

// Write extended header incl. CGAL/ENDCBP keywords.
std::ostream& operator<<( std::ostream& out,
                          const File_header_extended_OFF& h) {
    out << "#CBP\n";
    OUTBOOL( polyhedral_surface);
    OUT(     halfedges);
    OUTBOOL( triangulated);
    OUTBOOL( non_empty_facets);
    OUTBOOL( terrain);
    OUTBOOL( normalized_to_sphere);
    OUT(     radius);
    OUTBOOL( rounded);
    OUT(     rounded_bits);
    out << "# ENDCBP\n" << std::endl;
    return out;
}
#undef OUT
#undef OUTBOOL

#define IN(item,type)                         \
    else if ( CGAL_CLIB_STD::strcmp( keyword, #item) == 0) { \
        type t;                               \
        in >> t;                              \
        h.set_##item( t);                     \
    }

#define INBOOL(item)                          \
    else if ( CGAL_CLIB_STD::strcmp( keyword, #item) == 0) { \
        in >> c;                              \
        h.set_##item( c == '1');              \
    }

// Scan extended header. The CBP keyword must be read already.
std::istream& operator>>( std::istream& in, File_header_extended_OFF& h) {
    const int max_keyword = 42;
    char c;
    char keyword[max_keyword] = "";
    in >> keyword;
    while ( in && CGAL_CLIB_STD::strcmp( keyword, "ENDCBP") != 0) {
        if ( CGAL_CLIB_STD::strcmp( keyword, "#") == 0)
            ;
        INBOOL( polyhedral_surface)
        IN(     halfedges, int)
        INBOOL( triangulated)
        INBOOL( non_empty_facets)
        INBOOL( terrain)
        INBOOL( normalized_to_sphere)
        IN(     radius, double)
        INBOOL( rounded)
        IN(     rounded_bits, int)
        else if ( h.verbose()) {
            std::cerr << "warning: File_header_extended_OFF: unknown key '"
                      << keyword << "'." << std::endl;
        }
        in >> keyword;
    }
    in >> skip_until_EOL >> skip_comment_OFF;
    return in;
}
#undef IN
#undef INBOOL

CGAL_END_NAMESPACE
// EOF //
