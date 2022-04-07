// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_POLYHEDRON_VRML_2_OSTREAM_H
#define CGAL_IO_POLYHEDRON_VRML_2_OSTREAM_H

#include <CGAL/license/Polyhedron.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/IO/VRML.h>
#include <CGAL/boost/graph/IO/Generic_facegraph_printer.h>

namespace CGAL {

template <class Traits,
          class Items,
          template <class T, class I, class A>
          class HDS, class Alloc>
VRML_2_ostream& operator<<(VRML_2_ostream& out,
                           const Polyhedron_3<Traits, Items, HDS, Alloc>& P)
{
  IO::internal::Generic_facegraph_printer<VRML_2_ostream,
                                          Polyhedron_3<Traits, Items, HDS, Alloc>,
                                          File_writer_VRML_2> printer(out);
  printer(P);

  return out;
}

} // namespace CGAL

#endif // CGAL_IO_POLYHEDRON_VRML_2_OSTREAM_H
