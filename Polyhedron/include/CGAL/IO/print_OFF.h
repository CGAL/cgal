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

#ifndef CGAL_POLYHEDRON_IO_PRINT_OFF_H
#define CGAL_POLYHEDRON_IO_PRINT_OFF_H

#include <CGAL/license/Polyhedron.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/IO/OFF.h>
#include <CGAL/boost/graph/IO/Generic_facegraph_printer.h>

#include <fstream>

namespace CGAL {

template <class Polyhedron, class VPM>
bool print_polyhedron_with_header_OFF(std::ostream& out,
                                      const Polyhedron& P,
                                      const File_header_OFF& header,
                                      const VPM& vpm)
{
  File_writer_OFF writer(header);
  writer.header().set_polyhedral_surface(true);
  writer.header().set_halfedges(P.size_of_halfedges());

  IO::internal::Generic_facegraph_printer<std::ostream,
                                          Polyhedron,
                                          File_writer_OFF> printer(out, writer);

  return printer(P, parameters::vertex_point_map(vpm));
}

template <class Polyhedron>
bool print_polyhedron_with_header_OFF(std::ostream& out,
                                      const Polyhedron& P,
                                      const File_header_OFF& header)
{
  return print_polyhedron_with_header_OFF(out, P, header, get(CGAL::vertex_point, P));
}

template <class Polyhedron>
bool print_polyhedron_OFF(std::ostream& out,
                          const Polyhedron& P,
                          bool verbose = false)
{
  File_header_OFF header(verbose);
  header.set_binary(IO::is_binary(out));
  header.set_no_comments(!IO::is_pretty(out));

  return print_polyhedron_with_header_OFF(out, P, header);
}

} // namespace CGAL

#endif // CGAL_POLYHEDRON_IO_PRINT_OFF_H
