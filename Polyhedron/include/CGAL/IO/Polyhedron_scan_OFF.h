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

#ifndef CGAL_IO_POLYHEDRON_SCAN_OFF_H
#define CGAL_IO_POLYHEDRON_SCAN_OFF_H

#include <CGAL/license/Polyhedron.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/IO/OFF.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <iostream>
#include <cstddef>

namespace CGAL {

template < class HDS>
class Polyhedron_scan_OFF
  : public Modifier_base<HDS>
{
protected:
  std::istream& m_in;
  File_header_OFF m_file_header;

public:
  typedef HDS Halfedge_data_structure;

  // DEFINITION
  //
  // Polyhedron_scan_OFF<Traits> is a polyhedral surface builder.
  // It scans a polyhedron given in OFF from a stream and appends it
  // incrementally using the incremental builder.

  Polyhedron_scan_OFF( std::istream& in, bool verbose = false)
    : m_in(in), m_file_header(verbose)
  { }

  // Activation
  void operator()( HDS& hds);

  const File_header_OFF& header() const { return m_file_header; }
};

template < class HDS >
void Polyhedron_scan_OFF<HDS>:: operator()(HDS& target)
{
  File_scanner_OFF scanner( m_in, m_file_header.verbose());
  if(! m_in )
  {
    if(scanner.verbose())
    {
      std::cerr << " " << std::endl;
      std::cerr << "Polyhedron_scan_OFF<HDS>::" << std::endl;
      std::cerr << "operator(): input error: file format is not in "
                   "OFF." << std::endl;
    }
    m_in.clear( std::ios::badbit);
    return;
  }
  m_file_header = scanner; // Remember file header after return.

  Polyhedron_incremental_builder_3<HDS> B( target, scanner.verbose());
  B.begin_surface( scanner.size_of_vertices(),
                   scanner.size_of_facets(),
                   scanner.size_of_halfedges());

  typedef typename HDS::Traits     Traits;
  typedef typename Traits::Point_3 Point;

  // read in all vertices
  std::size_t  i;
  for(i=0; i<scanner.size_of_vertices(); ++i)
  {
    Point p;
    file_scan_vertex( scanner, p);
    B.add_vertex( p);
    if(scanner.has_vcolors())
    {
      IO::Color c;
      file_scan_color(scanner, c);
    }
      scanner.skip_to_next_vertex(i);
  }

  if(!m_in || B.error())
  {
    B.rollback();
    m_in.clear( std::ios::badbit);
    return;
  }

  // read in all facets
  for(i=0; i<scanner.size_of_facets(); ++i)
  {
    B.begin_facet();
    std::size_t no = 0;
    scanner.scan_facet( no, i);
    if(! m_in || B.error() || no < 3)
    {
      if(scanner.verbose())
      {
        std::cerr << " " << std::endl;
        std::cerr << "Polyhedron_scan_OFF<Traits>::" << std::endl;
        std::cerr << "operator()(): input error: facet " << i
                  << " has fewer than 3 vertices." << std::endl;
      }

      B.rollback();
      m_in.clear( std::ios::badbit);
      return;
    }

    for(std::size_t j=0; j<no; ++j)
    {
      std::size_t index = 0;
      scanner.scan_facet_vertex_index( index, j+1, i); //current_entry = j + 1 for the size entry
      if(!m_in)
        return;
      B.add_vertex_to_facet( index);
    }
    B.end_facet();
    scanner.skip_to_next_facet( i);
  }

  if(! m_in  || B.error())
  {
    B.rollback();
    m_in.clear( std::ios::badbit);
    return;
  }

  if(B.check_unconnected_vertices())
  {
    if(! B.remove_unconnected_vertices())
    {
      if(scanner.verbose()) {
        std::cerr << " " << std::endl;
        std::cerr << "Polyhedron_scan_OFF<Traits>::" << std::endl;
        std::cerr << "operator()(): input error: cannot "
                     "successfully remove isolated vertices."
                  << std::endl;
      }

      B.rollback();
      m_in.clear( std::ios::badbit);
      return;
    }
  }

  B.end_surface();
}

} // namespace CGAL

#endif // CGAL_IO_POLYHEDRON_SCAN_OFF_H
