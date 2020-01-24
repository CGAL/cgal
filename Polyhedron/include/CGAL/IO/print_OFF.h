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

#ifndef CGAL_IO_PRINT_OFF_H
#define CGAL_IO_PRINT_OFF_H

#include <CGAL/license/Polyhedron.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/IO/OFF.h>
#include <CGAL/IO/generic_print_polyhedron.h>

#include <fstream>

namespace CGAL {

template <class Polyhedron, class Vpm>
void print_polyhedron_with_header_OFF(std::ostream& out,
                                      const Polyhedron& P,
                                      const File_header_OFF& header,
                                      const Vpm& vpm)
{
  File_writer_OFF writer(header);
  writer.header().set_polyhedral_surface(true);
  writer.header().set_halfedges(P.size_of_halfedges());
  generic_print_polyhedron(out, P, writer, vpm);
}

template <class Polyhedron>
void print_polyhedron_with_header_OFF(std::ostream& out,
                                      const Polyhedron& P,
                                      const File_header_OFF& header)
{
  print_polyhedron_with_header_OFF(out, P, header, get(CGAL::vertex_point, P));
}

template <class Polyhedron>
void print_polyhedron_OFF(std::ostream& out,
                          const Polyhedron& P,
                          bool verbose = false)
{
  File_header_OFF header(verbose);
  header.set_binary(is_binary(out));
  header.set_no_comments(!is_pretty(out));
  print_polyhedron_with_header_OFF(out, P, header);
}

} //namespace CGAL

#endif // CGAL_IO_PRINT_OFF_H
