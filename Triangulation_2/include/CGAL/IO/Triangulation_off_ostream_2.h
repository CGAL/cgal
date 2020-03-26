// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Clement Jamin


#ifndef CGAL_TRIANGULATION_OFF_OSTREAM_2_H
#define CGAL_TRIANGULATION_OFF_OSTREAM_2_H

#include <CGAL/Triangulation_2.h>
#include <sstream>
#include <iostream>

namespace CGAL {

template < class GT, class TDS >
std::ostream &
export_triangulation_2_to_off(std::ostream & os,
                              const Triangulation_2<GT,TDS> & tr)
{
  typedef Triangulation_2<GT,TDS>                       Tr;
  typedef typename Tr::Vertex_handle                    Vertex_handle;
  typedef typename Tr::Vertex_iterator                  Vertex_iterator;
  typedef typename Tr::Finite_vertices_iterator         Finite_vertex_iterator;
  typedef typename Tr::Finite_faces_iterator            Finite_faces_iterator;

  size_t n = tr.number_of_vertices();

  std::stringstream output;

  // write the vertices
  std::map<Vertex_handle, int> index_of_vertex;
  int i = 0;
  for(Finite_vertex_iterator it = tr.finite_vertices_begin();
      it != tr.finite_vertices_end(); ++it, ++i)
  {
    output << it->point().x() << " " << it->point().y() << " 0" << std::endl;
    index_of_vertex[it.base()] = i;
  }
  CGAL_assertion( i == n );

  size_t number_of_triangles = 0;

  for (Finite_faces_iterator fit = tr.finite_faces_begin() ;
        fit != tr.finite_faces_end() ; ++fit)
  {
    output << "3 "
           << index_of_vertex[fit->vertex(0)] << " "
           << index_of_vertex[fit->vertex(1)] << " "
           << index_of_vertex[fit->vertex(2)]
           << std::endl;
    ++number_of_triangles;
  }

  os << "OFF \n"
     << n << " "
     << number_of_triangles << " 0\n"
     << output.str();

  return os;
}

} //namespace CGAL

#endif // CGAL_TRIANGULATION_OFF_OSTREAM_2_H