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


#ifndef CGAL_IO_FILE_HEADER_OFF_H
#define CGAL_IO_FILE_HEADER_OFF_H 1

#include <CGAL/IO/File_header_extended_OFF.h>
#include <iostream>

namespace CGAL {

// Info structure for OFF file headers
// ===================================
class CGAL_EXPORT File_header_OFF : public File_header_extended_OFF {
private:
    // Publicly accessible file informations.
    std::size_t  n_vertices;
  std::size_t n_facets;
    bool m_skel;        // SKEL format instead of OFF.
    bool m_binary;      // OFF in binary format.
    bool m_no_comments; // no comments in output.
  std::size_t  m_offset;      // index offset for vertices, usually 0.

    // Publicly accessible but not that well supported file informations.
    bool m_colors;      // COFF detected.
    bool m_normals;     // NOFF format stores also normals at vertices.

    // More privately used file informations to scan the file.
    bool m_tag4;        // 4OFF detected.
    bool m_tagDim;      // nOFF detected (will not be supported).
    int  m_dim;         // dimension for nOFF (will not be supported).
public:
    typedef  File_header_OFF           Self;
    typedef  File_header_extended_OFF  Base;

    explicit File_header_OFF( bool verbose = false);
    File_header_OFF( bool binary, bool noc, bool skel,
                     bool verbose = false);
    //File_header_OFF( int v, int h, int f, bool verbose = false);
    File_header_OFF( std::size_t v, std::size_t h, std::size_t f,
                     bool binary, bool noc, bool skel,
                     bool verbose = false);
    File_header_OFF( const File_header_extended_OFF& ext_header);
    File_header_OFF( const File_header_extended_OFF& ext_header,
                     bool binary, bool noc, bool skel);
    File_header_OFF( std::size_t v, std::size_t h, std::size_t f,
                     const File_header_extended_OFF& ext_header);
    File_header_OFF( std::size_t v, std::size_t h, std::size_t f,
                     const File_header_extended_OFF& ext_header,
                     bool binary, bool noc, bool skel);

    Self& operator= ( const Base& base) { (Base&)(*this) = base;
                                          return *this;
                                        }
    std::size_t  size_of_vertices()   const { return n_vertices; }
    std::size_t  size_of_facets()     const { return n_facets; }

    bool skel()               const { return m_skel; }   // SKEL format.
    bool off()                const { return ! m_skel; } // OFF format.
    bool binary()             const { return m_binary; } // binary format.
    bool ascii()              const { return ! m_binary; } // ASCII format.
    bool no_comments()        const { return m_no_comments; }
    bool comments()           const { return ! m_no_comments; }

    std::size_t  index_offset()       const { return m_offset; }
    bool has_colors()         const { return m_colors; } // COFF detected.
    bool has_normals()        const { return m_normals;} // NOFF format.
    bool is_homogeneous()     const { return m_tag4; }   // 4OFF detected.
                           // nOFF detected. (will not be supported).
    bool n_dimensional()      const { return m_tagDim; }
                           // dimension for nOFF (will not be supported).
    int  dimension()          const { return m_dim; }

    void set_vertices( std::size_t n)       { n_vertices = n; }
    void set_facets( std::size_t n)         { n_facets   = n; }

    void set_skel( bool b)          { m_skel        = b; }
    void set_binary( bool b)        { m_binary      = b; }
    void set_no_comments( bool b)   { m_no_comments = b; }
    void set_index_offset( std::size_t i)   { m_offset      = i; }

    void set_colors( bool b)        { m_colors      = b; }
    void set_normals( bool b)       { m_normals     = b;}
    void set_homogeneous( bool b)   { m_tag4        = b; }
    void set_dimensional( bool b)   { m_tagDim      = b; }
    void set_dimension( int i)      { m_dim         = i; }
    Self& operator+=( const Self& header);
};

// Write header.
CGAL_EXPORT std::ostream& operator<<( std::ostream& out, const File_header_OFF& h);

// Scan header. Marks streams badbit if not in SKEL format nor in OFF.
CGAL_EXPORT std::istream& operator>>( std::istream& in, File_header_OFF& h);

} //namespace CGAL
#endif // CGAL_IO_FILE_HEADER_OFF_H //
// EOF //
