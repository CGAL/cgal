// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_FILE_MEDIT_H
#define CGAL_FILE_MEDIT_H

#include <iostream>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

namespace CGAL {

template < class Tr>
void
output_pslg_to_medit (std::ostream& os, const Tr & T) {
  typedef typename Tr::Finite_cells_iterator Finite_cells_iterator;
  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Point Point;

  // Header.

  os << "MeshVersionFormatted 1" << std::endl
     << "Dimension 3" << std::endl;

  // Vertices
  
  os << "Vertices" << std::endl
     << T.number_of_vertices() << std::endl;

  os << std::setprecision(20);
 
  std::map<Vertex_handle, int> V;
  int inum = 1;
  for( Finite_vertices_iterator
      vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit) 
    {
      V[vit] = inum++;
      Point p = static_cast<Point>(vit->point());
      os << p.x() << " " << p.y() << " " << p.z() << " " 
	 << 0 << std::endl; //ref
    }

  // Facets
  os << "Triangles" << std::endl
     << number_of_facets_on_surface(T) << std::endl;
  for( Finite_facets_iterator fit = T.finite_facets_begin(); 
       fit != T.finite_facets_end(); ++fit)
    if ((*fit).first->is_facet_on_surface((*fit).second)==true)
      {
	for (int i=0; i<4; i++)
          if (i != (*fit).second)
	    os << V[(*fit).first->vertex(i)] << " ";
	
	os << "0" << std::endl; //ref
      }

  // Tetrahedra
  os << "Tetrahedra" << std::endl
     << number_of_cells_in_domain(T) << std::endl;
  for( Finite_cells_iterator cit = T.finite_cells_begin(); 
       cit != T.finite_cells_end(); ++cit)
    if( cit->is_in_domain() )
      {
	for (int i=0; i<4; i++)
	  os << V[cit->vertex(i)] << " ";
	os << "1" << std::endl; //ref
      }
  
  // End
  os << "End" << std::endl;
}

template < class Tr>
int number_of_cells_in_domain(const Tr& T) {
  int result=0;
  for (typename Tr::Finite_cells_iterator cit = T.finite_cells_begin(); 
       cit != T.finite_cells_end(); ++cit)
    if (cit->is_in_domain ())
      ++result;
  return result;
}

} // end namespace CGAL

#endif // CGAL_FILE_MEDIT_H
