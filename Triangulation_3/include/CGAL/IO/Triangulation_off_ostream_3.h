// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:  $
// $Id:  $
//
// Author(s)     : Clement Jamin


#ifndef CGAL_TRIANGULATION_OFF_OSTREAM_3_H
#define CGAL_TRIANGULATION_OFF_OSTREAM_3_H

#include <CGAL/license/Triangulation.h>

#include <CGAL/Triangulation_3.h>
#include <sstream>
#include <iostream>

namespace CGAL {

template < class GT, class TDS >
std::ostream &
export_triangulation_3_to_off(std::ostream & os, 
                              const Triangulation_3<GT,TDS> & tr,
                              bool export_surface_only = false)
{
  typedef Triangulation_3<GT,TDS>                       Tr;
  typedef typename Tr::Vertex_handle                    Vertex_handle;
  typedef typename Tr::Finite_vertices_iterator         Finite_vertex_iterator;
  typedef typename Tr::All_cells_iterator               Cells_iterator;
  typedef typename Tr::Finite_cells_iterator            Finite_cells_iterator;

  size_t n = tr.number_of_vertices();

  std::stringstream output;
  
  // write the vertices
  std::map<Vertex_handle, int> index_of_vertex;
  std::size_t i = 0;
  for(Finite_vertex_iterator it = tr.finite_vertices_begin(); 
      it != tr.finite_vertices_end(); ++it, ++i)
  {
    output << it->point().x() << " " 
           << it->point().y() << " " 
           << it->point().z() << std::endl;
    index_of_vertex[it.base()] = i;
  }
  CGAL_assertion( i == n );
  
  size_t number_of_triangles = 0;

  if (export_surface_only)
  {
    for (Cells_iterator cit = tr.cells_begin() ;
          cit != tr.cells_end() ; ++cit)
    {
      if (tr.is_infinite(cit))
      {
        output << "3 ";
        for (int i = 0 ; i < 4 ; ++i)
        {
          if (!tr.is_infinite(cit->vertex(i)))
            output << index_of_vertex[cit->vertex(i)] << " ";
        }
        output << std::endl;
        ++number_of_triangles;
      }
    }
  }
  else
  {
    for (Finite_cells_iterator cit = tr.finite_cells_begin() ;
          cit != tr.finite_cells_end() ; ++cit)
    {
      output << "3 "
             << index_of_vertex[cit->vertex(0)] << " "
             << index_of_vertex[cit->vertex(1)] << " "
             << index_of_vertex[cit->vertex(2)]
             << std::endl;
      output << "3 "
             << index_of_vertex[cit->vertex(0)] << " "
             << index_of_vertex[cit->vertex(2)] << " "
             << index_of_vertex[cit->vertex(3)]
             << std::endl;
      output << "3 "
             << index_of_vertex[cit->vertex(1)] << " "
             << index_of_vertex[cit->vertex(2)] << " "
             << index_of_vertex[cit->vertex(3)]
             << std::endl;
      output << "3 "
             << index_of_vertex[cit->vertex(0)] << " "
             << index_of_vertex[cit->vertex(1)] << " "
             << index_of_vertex[cit->vertex(3)]
             << std::endl;
      number_of_triangles += 4;
    }
  }

  os << "OFF \n"
     << n << " "
     << number_of_triangles << " 0\n"
     << output.str();

  return os;
}

} //namespace CGAL

#endif // CGAL_TRIANGULATION_OFF_OSTREAM_3_H
