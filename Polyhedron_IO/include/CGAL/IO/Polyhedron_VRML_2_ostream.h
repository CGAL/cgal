// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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

#ifndef CGAL_IO_POLYHEDRON_VRML_2_OSTREAM_H
#define CGAL_IO_POLYHEDRON_VRML_2_OSTREAM_H 1

#include <CGAL/basic.h>
#include <CGAL/IO/VRML_2_ostream.h>
#include <CGAL/IO/File_writer_VRML_2.h>
#include <CGAL/IO/generic_print_polyhedron.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>

namespace CGAL {

template < class Traits,
           class Items,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
           class HDS, class Alloc>
VRML_2_ostream&
operator<<( VRML_2_ostream& out,
            const Polyhedron_3<Traits,Items,HDS,Alloc>& P) {
    File_writer_VRML_2  writer;
    generic_print_polyhedron( out.os(), P, writer);
    return out;
}

} //namespace CGAL
#endif // CGAL_IO_POLYHEDRON_VRML_2_OSTREAM_H //
// EOF //
