// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Tran Kai Frank DA

#ifndef CGAL_ALPHA_SHAPE_CELL_BASE_3_H
#define CGAL_ALPHA_SHAPE_CELL_BASE_3_H

#include <vector>
#include <CGAL/utility.h>
#include <CGAL/Triangulation_cell_base_3.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Cb = Triangulation_cell_base_3<Gt> >
class Alpha_shape_cell_base_3
  : public Cb
{
public:
  typedef typename Cb::Vertex_handle   Vertex_handle;
  typedef typename Cb::Cell_handle     Cell_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other   Cb2;
    typedef Alpha_shape_cell_base_3<Gt, Cb2>                Other;
  };

  typedef typename Gt::FT                            Coord_type;
  typedef Triple<Coord_type, Coord_type, Coord_type> Interval3;

private:

  Interval3 vec_facet[4];
  Interval3 vec_edge[4][4];
  Coord_type A;

public:
  
  Alpha_shape_cell_base_3() 
    : Cb() {}
  
  Alpha_shape_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                          Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3) {}
  
  Alpha_shape_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                          Vertex_handle v2, Vertex_handle v3,
                          Cell_handle n0, Cell_handle n1,
                          Cell_handle n2, Cell_handle n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3) {}


  const Coord_type & get_alpha() const
    {
      return A;
    }
  
  void set_alpha(const Coord_type & AA)
    {
      A = AA;
    }

  const Interval3 & get_facet_ranges(int i) const
    {
      return vec_facet[i];
    }

  void set_facet_ranges(int i, const Interval3& Inter)
    {
      vec_facet[i]=Inter;
    }
  
  const Interval3 & get_edge_ranges(int i, int j) const
    {
      return vec_edge[i][j];
    }

  void set_edge_ranges(int i, int j, const Interval3& Inter)
    {
      vec_edge[i][j]=Inter;
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_ALPHA_SHAPE_CELL_BASE_3_H
