// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#include <CGAL/basic.h>
#include <cstdlib>
#include <iostream>
#include <CGAL/IO/binary_file_io.h>
#include <CGAL/IO/File_scanner_OFF.h>

namespace CGAL {

void
File_scanner_OFF::
skip_to_next_vertex( std::size_t current_vertex) {
    CGAL_assertion( current_vertex < size_of_vertices());
    if ( binary()) {
        float f;
        if ( has_normals() && ! normals_read) {
            I_Binary_read_big_endian_float32( m_in, f);
            I_Binary_read_big_endian_float32( m_in, f);
            I_Binary_read_big_endian_float32( m_in, f);
            if ( is_homogeneous())
                I_Binary_read_big_endian_float32( m_in, f);
        }
        if ( has_colors()) {
            // It is not well stated in the Geomview manual
            // how color is coded following a vertex. It is
            // parsed similar to the optional color for facets.
	    boost::int32_t k;
            I_Binary_read_big_endian_integer32( m_in, k);
            if (k<0 || k>4) {
                m_in.clear( std::ios::badbit);
                if ( verbose()) {
                    std::cerr << " " << std::endl;
                    std::cerr << "File_scanner_OFF::" << std::endl;
                    std::cerr << "skip_to_next_vertex(): input error: bad "
                                 " number of color indices at vertex "
                              << current_vertex << "." << std::endl;
                }
                set_off_header( false);
                return;
            }
            while (k--) {
                float dummy;
                I_Binary_read_big_endian_float32( m_in, dummy);
            }
        }
    } else {
        if ( has_normals() && ! normals_read) {
            double dummy;
            if ( is_homogeneous()) {
                m_in >> dummy >> dummy >> dummy >> dummy;
            } else {
                m_in >> dummy >> dummy >> dummy;
            }
        }
        if ( has_colors()) { // skip color entries (1 to 4)
            m_in >> skip_until_EOL;
        }
    }
    if( ! m_in) {
        if ( verbose()) {
            std::cerr << " " << std::endl;
            std::cerr << "File_scanner_OFF::" << std::endl;
            std::cerr << "skip_to_next_vertex(): input error: cannot read "
                         "OFF file beyond vertex " << current_vertex << "."
                      << std::endl;
        }
        set_off_header( false);
        return;
    }
    normals_read = false;
}

void
File_scanner_OFF::
skip_to_next_facet( std::size_t current_facet) {
    // Take care of trailing informations like color triples.
    if ( binary()) {
        boost::int32_t k;
        I_Binary_read_big_endian_integer32( m_in, k);
        if (k<0 || k>4) {
            m_in.clear( std::ios::badbit);
            if ( verbose()) {
                std::cerr << " " << std::endl;
                std::cerr << "File_scanner_OFF::" << std::endl;
                std::cerr << "skip_to_next_facet(): input error: bad "
                             "number of color indices at vertex "
                          << current_facet << "." << std::endl;
            }
            set_off_header( false);
            return;
        }
        while (k--) {
            float dummy;
            I_Binary_read_big_endian_float32( m_in, dummy);
        }
    } else {
        m_in >> skip_until_EOL;
    }
}

} //namespace CGAL
// EOF //
