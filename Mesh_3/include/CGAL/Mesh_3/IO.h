// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_MESH_3_IO_H
#define CGAL_MESH_3_IO_H

#include <CGAL/IO/File_medit.h> // for Debug
#include <iostream>
#include <map>

namespace CGAL { namespace Mesh_3 {

/** Ouput a C2T3 corresponding to a mesh to a std::ostream.
    This function assumes that the dimension of the triangulation is 3.
*/
template <class C2T3>
void output_mesh(std::ostream& os, const C2T3& c2t3)
{
  typedef typename C2T3::Triangulation_3 Tr;
  typedef typename Tr::Cell_iterator Cell_iterator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Point Point;

  const Tr& tr = c2t3.triangulation();
  os << tr;


//   int number_of_vertices = tr.number_of_vertices();

//   os << number_of_vertices;
//   if (is_ascii(os))
//     os << std::endl;

//   std::map<Vertex_handle, int> indices_of_vertices;

//   indices_of_vertices[tr.infinite_vertex()] = 0;

//   int i = 0;
//   for (Finite_vertices_iterator vit=tr.finite_vertices_begin();
//        vit!=tr.finite_vertices_end();
//        ++vit)
//   {
//     indices_of_vertices[vit] = ++i;
//     os << *vit;
//     if (is_ascii(os))
// 	os << std::endl;
//   }

//   tr.tds().print_cells(os, indices_of_vertices);

//   for(Cell_iterator cit=tr.cells_begin(); cit != tr.cells_end(); ++cit)
//   { 
//     os << *cit;
//     if (is_ascii(os))
//       os << std::endl;
//   }
} // end output_mesh

template <class C2T3>
bool
input_mesh(std::istream& is, 
           C2T3 & c2t3,
           bool debug = false, 
           std::ostream* debug_str = &std::cout)
{
  typedef typename C2T3::Triangulation_3 Tr;
  typedef typename Tr::Triangulation_data_structure Tds;
  typedef typename Tr::Vertex_handle Vertex_handle;  
  typedef typename Tr::Cell_handle Cell_handle;  

  details::Debug debug_stream(debug,
                              debug_str,
                              "CGAL::Mesh_3::input_mesh()"
                              " input error: ");

  Tr& tr = c2t3.triangulation();
  is >> tr;
  return (is.good());

//   Tds& tds = tr.tds();

//   tr.clear();

//   tds.set_dimension(3);

//   int number_of_vertices;
//   is >> number_of_vertices;

//   if( !is)
//     return debug_stream << "Cannot read number of vertices!\n";

//   std::map<int, Vertex_handle> vertex_handles;

//   vertex_handles[0] = tr.infinite_vertex();

//   for(int i = 1; i <= number_of_vertices; ++i)
//   {
//     vertex_handles[i] = tds.create_vertex();
//     is >> *vertex_handles[i];
//     if( !is )
//       return debug_stream << "Cannot read vertex:" << i << "!\n";
//   }

//   std::map<int, Cell_handle > C;

//   int m;
//   tds.read_cells(is, vertex_handles, m, C);

//   if( !is )
//     return debug_stream << "Cannot read cells!\n";

//   std::cerr << "Reading infos\n";
//   for (int j=0 ; j < m; j++)
//   {
//     is >> *(C[j]);
//     if( !is )
//       return debug_stream << "Cannot read informations for cell:"
//                           << j << "!\n";
//   }

//    CGAL_triangulation_assertion( tr.is_valid(true, 3) );

//   return is;
} // end input_mesh

} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_IO_H
