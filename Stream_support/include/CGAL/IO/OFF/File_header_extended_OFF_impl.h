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

#include <CGAL/IO/OFF/File_header_extended_OFF.h>
#include <CGAL/basic.h>

#include <cstdlib>
#include <cctype>
#include <cstring>
#include <iostream>
#include <sstream>
#include <algorithm>


namespace CGAL {

CGAL_INLINE_FUNCTION
bool File_header_extended_OFF::
is_POL()  const {
    return is_OFF() && polyhedral_surface();
}

CGAL_INLINE_FUNCTION
bool File_header_extended_OFF::
is_CBP()  const {
    return is_POL() && triangulated() && non_empty_facets() &&
        normalized_to_sphere() && radius() <= 1.0;
}

CGAL_INLINE_FUNCTION
bool File_header_extended_OFF::
is_TRN()  const { return is_CBP() && terrain(); }

CGAL_INLINE_FUNCTION
int  File_header_extended_OFF::
is_CBPn() const {
    if ( is_POL() && triangulated() && non_empty_facets() &&
         normalized_to_sphere() && rounded() &&
         (radius() <= double( 1l << rounded_bits())))
        return rounded_bits();
    else
        return 0;
}

CGAL_INLINE_FUNCTION
int  File_header_extended_OFF::
is_TRNn() const { return ( terrain() ? is_CBPn() : 0); }


// The proper file suffix with respect to file format.
CGAL_INLINE_FUNCTION
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
CGAL_INLINE_FUNCTION
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

CGAL_INLINE_FUNCTION
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
    m_radius               = (std::max)(m_radius, header.m_radius);
    m_rounded              = m_rounded && header.m_rounded;
    m_rounded_bits         = (std::max)( m_rounded_bits,
                                       header.m_rounded_bits);
    m_off_header           = m_off_header && header.m_off_header;
    return *this;
}

} //namespace CGAL
// EOF //
