// Copyright (c) 1997-2010  INRIA Sophia-Antipolis (France).
// Copyright (c) 2011, 2020 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau, Maxime Gimeno

// Adapted from operator>>(std::istream&, Triangulation_3&) from
// <CGAL/Triangulation_3.h>


#ifndef CGAL_TR_OR_TDS_FILE_INPUT_H
#define CGAL_TR_OR_TDS_FILE_INPUT_H

#include <CGAL/license/TDS_3.h>

#include <CGAL/license/Triangulation_3.h>


#include <CGAL/basic.h>

namespace CGAL {
namespace internal{

template <typename Tr_src,
          typename Tr_tgt,
          typename ConvertVertex,
          typename ConvertCell>
std::istream& file_input(std::istream& is, Tr_tgt &tr, bool is_tds,
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

  // If this is used for a TDS, the vertices are processed from 0 to n.
  // Else, we make V[0] the infinite vertex and work from 1 to n+1.
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

  std::size_t V_size = is_tds ? n : n+1;
  std::vector< Vertex_handle > V(V_size);

  // the infinite vertex is numbered 0
  if(!is_tds)
    V[0] = tr.infinite_vertex();

  for (std::size_t i=is_tds ? 0 : 1; i < V_size; ++i) {
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
    convert_cell(c, *C[j]);
  }

  CGAL_triangulation_assertion( tr.is_valid(false) );
  return is;
}

} //end internal
} // end namespace CGAL


#endif // CGAL_TR_OR_TDS_FILE_INPUT_H
