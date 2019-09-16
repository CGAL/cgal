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
// Author(s)     : Laurent Rineau
//

// Adapted from operator>>(std::istream&, Triangulation_3&) from
// <CGAL/Triangulation_3.h>

#ifndef CGAL_TRIANGULATION_FILE_INPUT_3_H
#define CGAL_TRIANGULATION_FILE_INPUT_3_H

#include <CGAL/basic.h>

namespace CGAL {

template <typename Tr1, 
          typename Tr2,
          typename Update_vertex,
          typename Update_cell>
std::istream& file_input(std::istream& is, Tr2 &tr,
                         Update_vertex update_vertex = Update_vertex(),
                         Update_cell update_cell = Update_cell())
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
  typedef Tr2 Triangulation;
  typedef typename Triangulation::Vertex_handle  Vertex_handle;
  typedef typename Triangulation::Cell_handle    Cell_handle;

  typedef typename Tr1::Vertex Vertex1;
  typedef typename Tr1::Cell Cell1;

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
    V[i] = tr.tds().create_vertex();
    Vertex1 v;
    if(!(is >> v)) return is;
    if(!update_vertex(v, *V[i])) {
      is.setstate(std::ios_base::failbit);
      return is;
    }
  }

  std::vector< Cell_handle > C;

  std::size_t m;
  tr.tds().read_cells(is, V, m, C);

  for (std::size_t j=0 ; j < m; j++) {
    Cell1 c;
    if(!(is >> c)) return is;
    if(!update_cell(c, *(C[j]))) {
      is.setstate(std::ios_base::failbit);
      return is;
    }
  }

  CGAL_triangulation_assertion( tr.is_valid(false) );
  return is;
}

} // end namespace CGAL

#endif // CGAL_TRIANGULATION_FILE_INPUT_3_H
