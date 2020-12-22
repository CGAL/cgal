// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_PRINT_WAVEFRONT_H
#define CGAL_IO_PRINT_WAVEFRONT_H 1

#include <CGAL/license/Polyhedron.h>


#include <CGAL/IO/File_writer_wavefront.h>
#include <CGAL/IO/generic_print_polyhedron.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>

namespace CGAL {

template <class Polyhedron>
void print_polyhedron_wavefront( std::ostream& out, const Polyhedron& P) {
    File_writer_wavefront  writer;
    generic_print_polyhedron( out, P, writer);
}

// Deprecated global functions, replaced with functions above

template < class Traits,
           class Items,
           template < class T, class I, class A>
           class HDS, class Alloc>
void
print_wavefront( std::ostream& out,
                 const Polyhedron_3<Traits,Items,HDS,Alloc>& P) {
    File_writer_wavefront  writer;
    generic_print_polyhedron( out, P, writer);
}

} //namespace CGAL
#endif // CGAL_IO_PRINT_WAVEFRONT_H //
// EOF //
