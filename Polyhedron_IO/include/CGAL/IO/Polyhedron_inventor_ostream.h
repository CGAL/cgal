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

#ifndef CGAL_IO_POLYHEDRON_INVENTOR_OSTREAM_H
#define CGAL_IO_POLYHEDRON_INVENTOR_OSTREAM_H 1

#include <CGAL/license/Polyhedron.h>


#include <CGAL/IO/Inventor_ostream.h>
#include <CGAL/IO/File_writer_inventor.h>
#include <CGAL/IO/generic_print_polyhedron.h>
#include <CGAL/Polyhedron_3.h>

namespace CGAL {

template < class Traits,
           class Items,
           template < class T, class I, class A>
           class HDS, class Alloc>
Inventor_ostream_base&
operator<<( Inventor_ostream_base& out,
            const Polyhedron_3<Traits,Items,HDS,Alloc>& P) {
    File_writer_inventor  writer;
    generic_print_polyhedron( out.os(), P, writer);
    return out;
}

} //namespace CGAL
#endif // CGAL_IO_POLYHEDRON_INVENTOR_OSTREAM_H //
// EOF //
