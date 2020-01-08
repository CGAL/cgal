// Copyright (c) 1997-2010  INRIA Sophia-Antipolis (France).
// Copyright (c) 2011       GeometryFactory Sarl (France)
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
// Author(s)     : Laurent Rineau, Maxime Gimeno
//

// Adapted from operator>>(std::istream&, Triangulation_3&) from
// <CGAL/Triangulation_3.h>

#ifndef CGAL_TRIANGULATION_FILE_INPUT_3_H
#define CGAL_TRIANGULATION_FILE_INPUT_3_H

#include <CGAL/basic.h>

namespace CGAL {

template <typename Tr_src, 
          typename Tr_tgt,
          typename ConvertVertex,
          typename ConvertCell>
std::istream& file_input(std::istream& is, Tr_tgt &tr,
                         ConvertVertex convert_vertex = ConvertVertex(),
                         ConvertCell convert_cell = ConvertCell())
  // reads
  // the dimension
  // the number of finite vertices
  // the non combinatorial information on vertices (point, etc)
  // the number of cells
  // the cells by the indices of their vertices in the preceding list
  // of vertices, plus the non combinatorial information on each cell
  // the neighbors of each cell by their index in the preceding list of cells
  // when dimension < 3 : the same with faces of maximal dimension
{
  typedef Tr_tgt Triangulation;
  typedef typename Triangulation::Vertex_handle  Vertex_handle;
  typedef typename Triangulation::Cell_handle    Cell_handle;

  typedef typename Tr_src::Vertex Vertex1;
  typedef typename Tr_src::Cell Cell1;

  tr.clear();
  tr.tds().cells().clear();

  std::size_t n;
  int d;
  if(is_ascii(is))
     is >> d >> n;
  else {
    read(is, d);
    read(is, n);
  }
  if(!is) return is;
  tr.tds().set_dimension(d);

  std::vector< Vertex_handle > V(n+1);
  V[0] = tr.infinite_vertex();
  // the infinite vertex is numbered 0

  for (std::size_t i=1; i <= n; i++) {
    //V[i] = tr.tds().create_vertex();   
    Vertex1 v;
    if(!(is >> v)) return is;
    Vertex_handle vh=tr.tds().create_vertex( convert_vertex(v) );
    V[i] = vh;
    convert_vertex(v, *V[i]);
  }

  std::vector< Cell_handle > C;

  std::size_t m;
  tr.tds().read_cells(is, V, m, C);

  for (std::size_t j=0 ; j < m; j++) {
    Cell1 c;
    if(!(is >> c)) return is;
    Cell_handle ch=tr.tds().create_cell(convert_cell(c));
    C[j] = ch;
    convert_cell(c, *ch);
  }

  CGAL_triangulation_assertion( tr.is_valid(false) );
  return is;
}

} // end namespace CGAL


#endif // TRIANGULATION_FILE_INPUT_H
