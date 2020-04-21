// Copyright (c) 2009   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_PERIODIC_3_TRIANGULATION_3_IO_H
#define CGAL_PERIODIC_3_TRIANGULATION_3_IO_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/array.h>

#include <iostream>
#include <utility>

namespace CGAL {

template <class Stream, class Triangulation>
Stream &write_triangulation_to_off(Stream &out, Triangulation &t) {
  typedef typename Triangulation::Point Point;
  typedef typename Triangulation::Domain Domain;
  typedef typename Triangulation::Triangulation_data_structure::size_type size_type;

  size_type number_of_cells = t.tds().number_of_cells();

  out << "OFF "
      << "\n" << 4*number_of_cells
      << " "  << 4*number_of_cells
      << " "  << 0
      << std::endl;

  if (t.is_1_cover()) {
    for (typename Triangulation::Cell_iterator it = t.cells_begin();
         it != t.cells_end(); it++) {
      for (int i=0; i<4; i++) {
        Point p = t.point(t.periodic_point(it,i));
        out << p.x() << " "
            << p.y() << " "
            << p.z() << std::endl;
      }
    }
  } else {
    for (typename Triangulation::Cell_iterator it = t.cells_begin();
         it != t.cells_end(); it++) {
      for (int i=0; i<4; i++) {
        typename Triangulation::Vertex_handle vh;
        typename Triangulation::Offset off;
        t.get_vertex(it, i, vh, off);
        Point p = t.point(t.periodic_point(it, i));

        out << p.x() << " "
            << p.y() << " "
            << p.z() << std::endl;
      }
    }
  }

  for (size_type i=0; i<number_of_cells; i++) {
    out << "3 " << i*4   << " " << i*4+1 << " " << i*4+2 << std::endl;
    out << "3 " << i*4   << " " << i*4+1 << " " << i*4+3 << std::endl;
    out << "3 " << i*4   << " " << i*4+2 << " " << i*4+3 << std::endl;
    out << "3 " << i*4+1 << " " << i*4+2 << " " << i*4+3 << std::endl;
  }

  return out;
}

template<class Stream, class Triangulation, class Cell_iterator>
Stream &write_cells_to_off(Stream &out, Triangulation &t, int number_of_cells,
                           Cell_iterator cit, Cell_iterator cells_end) {
  typedef typename Triangulation::Point Point;
  out << "OFF "
      << "\n" << 4*number_of_cells
      << " "  << 4*number_of_cells
      << " "  << 0
      << std::endl;

  while (cit != cells_end) {
    for (int i=0; i<4; i++) {
      Point p = t.get_point(*cit,i);
      out << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }
    ++cit;
  }

  for (int i=0; i<number_of_cells; i++) {
    out << "3 " << i*4 << " " << i*4+1 << " " << i*4+2 << std::endl;
    out << "3 " << i*4 << " " << i*4+1 << " " << i*4+3 << std::endl;
    out << "3 " << i*4 << " " << i*4+2 << " " << i*4+3 << std::endl;
    out << "3 " << i*4+1 << " " << i*4+2 << " " << i*4+3 << std::endl;
  }

  return out;
}

template <class Stream, class Triangulation>
Stream& draw_dual_to_off(Stream &os, Triangulation &t) {
  os << "OFF " << "\n"
     << 2*t.number_of_facets() << " "
     << t.number_of_facets() << " 0" << std::endl;

  t.draw_dual(os);

  for(unsigned int i=0; i < t.number_of_facets(); i++) {
    os << "3 " << i*2 << " " << i*2+1 << " " << i*2 << std::endl;
  }
  return os;
}

} // namespace CGAL

#endif //CGAL_PERIODIC_3_TRIANGULATION_3_IO_H
