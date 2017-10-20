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
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/basic.h>
#include <cstdlib>
#include <cctype>
#include <cstring>
#include <iostream>
#include <CGAL/IO/binary_file_io.h>
#include <CGAL/IO/File_header_OFF.h>
#include <algorithm>
#include <boost/cstdint.hpp>

namespace CGAL {

CGAL_INLINE_FUNCTION
File_header_OFF::File_header_OFF( bool verbose)
:   File_header_extended_OFF( verbose),
    n_vertices(0),
    n_facets(0),
    m_skel(false),
    m_binary(false),
    m_no_comments(false),
    m_offset(0),
    m_colors(false),
    m_normals(false),
    m_tag4(false),
    m_tagDim(false),
    m_dim(3)
{}
CGAL_INLINE_FUNCTION
File_header_OFF::File_header_OFF(
                     bool binary, bool noc, bool skel, bool verbose)
:   File_header_extended_OFF( verbose),
    n_vertices(0),
    n_facets(0),
    m_skel(skel),
    m_binary(binary),
    m_no_comments(noc),
    m_offset(0),
    m_colors(false),
    m_normals(false),
    m_tag4(false),
    m_tagDim(false),
    m_dim(3)
{}
//CGAL_INLINE_FUNCTION
//File_header_OFF::File_header_OFF( int v, int h, int f,
//                                          bool verbose)
//:   File_header_extended_OFF( verbose),
//    n_vertices(v),
//    n_facets(f),
//    m_skel(false),
//    m_binary(false),
//    m_no_comments(false),
//    m_offset(0),
//    m_colors(false),
//    m_normals(false),
//    m_tag4(false),
//    m_tagDim(false),
//    m_dim(3)
//{
//    set_halfedges(h);
//}
CGAL_INLINE_FUNCTION
File_header_OFF::File_header_OFF( std::size_t v, std::size_t h, std::size_t f,
                     bool binary, bool noc, bool skel, bool verbose)
:   File_header_extended_OFF( verbose),
    n_vertices(v),
    n_facets(f),
    m_skel(skel),
    m_binary(binary),
    m_no_comments(noc),
    m_offset(0),
    m_colors(false),
    m_normals(false),
    m_tag4(false),
    m_tagDim(false),
    m_dim(3)
{
    set_halfedges(h);
}
CGAL_INLINE_FUNCTION
File_header_OFF::File_header_OFF(
                     const File_header_extended_OFF& ext_header)
:   File_header_extended_OFF( ext_header),
    n_vertices(0),
    n_facets(0),
    m_skel(false),
    m_binary(false),
    m_no_comments(false),
    m_offset(0),
    m_colors(false),
    m_normals(false),
    m_tag4(false),
    m_tagDim(false),
    m_dim(3)
{}
CGAL_INLINE_FUNCTION
File_header_OFF::File_header_OFF(
                     const File_header_extended_OFF& ext_header,
                     bool binary, bool noc, bool skel)
:   File_header_extended_OFF( ext_header),
    n_vertices(0),
    n_facets(0),
    m_skel(skel),
    m_binary(binary),
    m_no_comments(noc),
    m_offset(0),
    m_colors(false),
    m_normals(false),
    m_tag4(false),
    m_tagDim(false),
    m_dim(3)
{}
CGAL_INLINE_FUNCTION
File_header_OFF::File_header_OFF(
                                 std::size_t v, std::size_t h, std::size_t f,
                     const File_header_extended_OFF& ext_header)
:   File_header_extended_OFF( ext_header),
    n_vertices(v),
    n_facets(f),
    m_skel(false),
    m_binary(false),
    m_no_comments(false),
    m_offset(0),
    m_colors(false),
    m_normals(false),
    m_tag4(false),
    m_tagDim(false),
    m_dim(3)
{
    set_halfedges(h);
}
CGAL_INLINE_FUNCTION
File_header_OFF::File_header_OFF(
                                 std::size_t v, std::size_t h, std::size_t f,
                     const File_header_extended_OFF& ext_header,
                     bool binary, bool noc, bool skel)
:   File_header_extended_OFF( ext_header),
    n_vertices(v),
    n_facets(f),
    m_skel(skel),
    m_binary(binary),
    m_no_comments(noc),
    m_offset(0),
    m_colors(false),
    m_normals(false),
    m_tag4(false),
    m_tagDim(false),
    m_dim(3)
{
    set_halfedges(h);
}

CGAL_INLINE_FUNCTION
File_header_OFF& File_header_OFF::
operator+=( const File_header_OFF& header) {
    (File_header_extended_OFF&)(*this) = header;
    n_vertices += header.n_vertices;
    n_facets   += header.n_facets;
    return *this;
}

// Write header.
CGAL_INLINE_FUNCTION
std::ostream& operator<<( std::ostream& out, const File_header_OFF& h) {
    if ( h.comments()) {
        out << "# Output of a CGAL tool\n";
        out << static_cast<const File_header_extended_OFF&>( h);
    }
    if ( h.has_normals())
        out << 'N';
    if ( h.skel())
        out << "SKEL";
    else
        out << "OFF";
    if ( h.binary()) {
        out << " BINARY\n";
        I_Binary_write_big_endian_integer32( out, static_cast<int>(h.size_of_vertices()));
        I_Binary_write_big_endian_integer32( out, static_cast<int>(h.size_of_facets()));
        if ( h.off())
            I_Binary_write_big_endian_integer32( out, 0);
    } else {
        out << '\n';
        out << h.size_of_vertices() << ' '<< h.size_of_facets();
        if ( h.off())
            out << " 0";
        if ( h.comments()) {
            out << "\n\n# " << h.size_of_vertices() << " vertices\n";
            out << "# ------------------------------------------\n";
        }
        out << std::endl;
    }
    return out;
}

// Scan header. Marks streams badbit if not in SKEL format nor in OFF.
CGAL_INLINE_FUNCTION
std::istream& operator>>( std::istream& in, File_header_OFF& h) {
    // read in the first character and scan for comments, `OFF', or `NOFF',
    // or `SKEL', or `4SKEL'.
    h.set_off_header( false);
    char c;
    while ( (in >> c) && c == '#') {
        if ( in.get(c) && c == 'C' &&
             in.get(c) && c == 'B' &&
             in.get(c) && c == 'P') {
            in >> static_cast<File_header_extended_OFF&>( h);
        } else if ( c != '\n')
            in >> skip_until_EOL;
    }
    if ( ! in)
        return in;
    h.set_skel( false);
    h.set_binary( false);
    h.set_index_offset( 1);
    h.set_colors( false);
    h.set_normals( false);
    h.set_homogeneous( false);
    h.set_dimensional( false);
    h.set_dimension( 3);

    const int max_keyword = 42;
    char keyword[max_keyword] = "";
    int i = 0;
    keyword[i++] = c;
    while( i < max_keyword - 1 && in.get(c) && std::isalnum(c))
        keyword[i++] = c;
    keyword[i] = '\0';
    if ( i < 2 || (std::isdigit(keyword[0]) && keyword[0] != '4')
               || std::isdigit(keyword[1])) {
        h.set_vertices( std::atoi( keyword));
    } else {
        h.set_index_offset( 0);
        int j = 0;
        if ( j < i && keyword[j] == 'C') {
            h.set_colors( true);
            j++;
        }
        if ( j < i && keyword[j] == 'N') {
            h.set_normals( true);
            j++;
        }
        if ( j < i && keyword[j] == '4') {
            h.set_homogeneous( true);
            j++;
        }
        if ( j < i && keyword[j] == 'n') {
            h.set_dimensional( true);
            j++;
        }
        if ( i-j != 3 || keyword[j]   != 'O'
                      || keyword[j+1] != 'F'
                      || keyword[j+2] != 'F') {
            if ( i-j != 4 || keyword[j]   != 'S'
                          || keyword[j+1] != 'K'
                          || keyword[j+2] != 'E'
                          || keyword[j+3] != 'L') {
                in.clear( std::ios::badbit);
                if ( h.verbose()) {
                    std::cerr << " " << std::endl;
                    std::cerr << "error: File_header_OFF: "
                                  "wrong format: neither OFF nor SKEL."
                              << std::endl;
                }
                return in;
            } else {
                h.set_skel( true);
            }
        }
        in >> skip_comment_OFF >> c;
        if ( std::isdigit(c)) {
            in.putback(c);
            int n;
            in >> n;
            h.set_vertices(n);
        } else {
            i = 0;
            keyword[i++] = c;
            while( i < max_keyword - 1 && in.get(c) &&
                   std::isalnum(c))
                keyword[i++] = c;
            keyword[i] = '\0';
            if ( std::strcmp( keyword, "BINARY") == 0) {
                h.set_binary( true);
                if ( c != '\n')
                    in >> skip_until_EOL;
            } else {
                in.clear( std::ios::badbit);
                if ( h.verbose()) {
                    std::cerr << " " << std::endl;
                    std::cerr << "error: File_header_OFF(): "
                                 "wrong format: neither OFF nor SKEL."
                              << std::endl;
                }
                return in;
            }
        }
    }
    // Read remaining size value(s).
    int n_h;
    if ( h.binary()) {
        boost::int32_t a, b, c;
        I_Binary_read_big_endian_integer32( in, a);
        if ( h.n_dimensional()) {
            h.set_dimension( a);
            I_Binary_read_big_endian_integer32( in, a);
        }
        I_Binary_read_big_endian_integer32( in, b);
        if ( h.off())
            I_Binary_read_big_endian_integer32( in, c);
        else
            c = 0;
        h.set_vertices( a);
        if (b<0){
          in.clear( std::ios::badbit );
          if ( h.verbose()) {
              std::cerr << " " << std::endl;
              std::cerr << "error: File_header_OFF(): File contains < 0 facets."
                        << std::endl;
          }
          return in;
        }
        h.set_facets( b);
        n_h = c;
    } else {
        int n;
        if ( h.n_dimensional()) {
            h.set_dimension( static_cast<int>(h.size_of_vertices()));
            in >> n;
            h.set_vertices(n);
        }
        in >> n;
        if (n < 0){
          in.clear( std::ios::badbit );
          if ( h.verbose()) {
              std::cerr << " " << std::endl;
              std::cerr << "error: File_header_OFF(): File contains < 0 facets."
                        << std::endl;
          }
          return in;
        }
        h.set_facets(n);
        if ( h.off())
            in >> n_h;
        else
            n_h = 0;
    }
    if ( n_h == 0)
        h.set_index_offset( 0);
    if ( ! in || h.size_of_vertices() <= 0 ) {
        in.clear( std::ios::badbit);
        if ( h.verbose()) {
            std::cerr << " " << std::endl;
            std::cerr << "error: File_header_OFF(): File contains <= 0 vertices."
                      << std::endl;
        }
        return in;
    }
    if ( h.size_of_halfedges() == 0) {
        // be careful, because border edges count twice
        h.set_halfedges( 2 * n_h);
        // check against the Eulerian equation for connected planar graphs.
        // We do not know the number of holes we must represent as dummy
        // facets and we do not know the genus of the surface.
        // So we add 12 and a factor of 5 percent.
        if (    h.size_of_halfedges() == 0
             || double(h.size_of_halfedges()) > double(h.size_of_vertices()
                + h.size_of_facets() - 2 + 12) * 2.1
             || h.size_of_halfedges() < (h.size_of_vertices()
                + h.size_of_facets() - 2) * 2
        )
            h.set_halfedges( int( double(h.size_of_vertices() +
                                         h.size_of_facets() - 2 + 12) * 2.1
                                 )
                            );
    }
    h.set_off_header( h.off());
    return in;
}

} //namespace CGAL
// EOF //
