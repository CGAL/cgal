// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_FIXED_ALPHA_SHAPE_CELL_BASE_3_H
#define CGAL_FIXED_ALPHA_SHAPE_CELL_BASE_3_H

#include <CGAL/license/Alpha_shapes_3.h>


#include <vector>
#include <CGAL/Compact_container.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/internal/Classification_type.h>

namespace CGAL {

template < class Gt, class Cb = Delaunay_triangulation_cell_base_3<Gt> >
class Fixed_alpha_shape_cell_base_3
  : public Cb
{
public:
  typedef typename Cb::Vertex_handle   Vertex_handle;
  typedef typename Cb::Cell_handle     Cell_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other   Cb2;
    typedef Fixed_alpha_shape_cell_base_3<Gt, Cb2>                Other;
  };


private:
  typedef internal::Classification_type Classification_type;

  Classification_type facet_status[4];
  Classification_type status_;

public:

  Fixed_alpha_shape_cell_base_3()
    : Cb() {}

  Fixed_alpha_shape_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                          Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3) {}

  Fixed_alpha_shape_cell_base_3(Vertex_handle v0, Vertex_handle v1,
                          Vertex_handle v2, Vertex_handle v3,
                          Cell_handle n0, Cell_handle n1,
                          Cell_handle n2, Cell_handle n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3) {}


  Classification_type get_facet_classification_type(int i) const {return facet_status[i];}
  void set_facet_classification_type(int i, Classification_type status) { facet_status[i]=status;  }
  Classification_type get_classification_type() { return status_;}
  void set_classification_type(Classification_type status) {status_=status;}
};

} //namespace CGAL

#endif // CGAL_FIXED_ALPHA_SHAPE_CELL_BASE_3_H
