// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_PRINT_OFF_H
#define CGAL_IO_PRINT_OFF_H 1

#include <CGAL/license/Polyhedron.h>


#include <CGAL/basic.h>
#include <CGAL/IO/File_writer_OFF.h>
#include <CGAL/IO/generic_print_polyhedron.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>

namespace CGAL {

template <class Polyhedron>
void print_polyhedron_with_header_OFF( std::ostream& out, 
                                       const Polyhedron& P,
                                       const File_header_OFF& header) {
    File_writer_OFF  writer( header);
    writer.header().set_polyhedral_surface( true);
    writer.header().set_halfedges( P.size_of_halfedges());
    generic_print_polyhedron( out, P, writer);
}


template <class Polyhedron>
void print_polyhedron_OFF( std::ostream& out, 
                           const Polyhedron& P,
                           bool verbose = false) {
    File_header_OFF header( verbose);
    header.set_binary( is_binary( out));
    header.set_no_comments( ! is_pretty( out));
    print_polyhedron_with_header_OFF( out, P, header);
}


// Deprecated global functions, replaced with functions above

template < class Traits,
           class Items,
           template < class T, class I, class A>
           class HDS, class Alloc>
void
print_OFF( std::ostream& out,
           const Polyhedron_3<Traits,Items,HDS,Alloc>& P,
           bool verbose = false) {
    File_writer_OFF  writer( verbose);
    writer.header().set_binary( is_binary( out));
    writer.header().set_no_comments( ! is_pretty( out));
    writer.header().set_polyhedral_surface( true);
    writer.header().set_halfedges( P.size_of_halfedges());
    generic_print_polyhedron( out, P, writer);
}

template < class Traits,
           class Items,
           template < class T, class I, class A>
           class HDS, class Alloc>
void
print_OFF( std::ostream& out,
           const Polyhedron_3<Traits,Items,HDS,Alloc>& P,
           const File_header_OFF& header) {
    File_writer_OFF  writer( header);
    writer.header().set_polyhedral_surface( true);
    writer.header().set_halfedges( P.size_of_halfedges());
    generic_print_polyhedron( out, P, writer);
}

} //namespace CGAL
#endif // CGAL_IO_PRINT_OFF_H //
// EOF //
