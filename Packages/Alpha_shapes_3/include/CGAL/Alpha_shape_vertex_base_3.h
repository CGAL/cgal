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

#ifndef CGAL_ALPHA_SHAPE_VERTEX_BASE_3_H
#define CGAL_ALPHA_SHAPE_VERTEX_BASE_3_H

#include <utility>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>

CGAL_BEGIN_NAMESPACE

template <class Gt, class Vb = Triangulation_vertex_base_3<Gt> >
class Alpha_shape_vertex_base_3
  : public Vb
{
public:

  typedef typename Vb::Cell_handle    Cell_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other   Vb2;
    typedef Alpha_shape_vertex_base_3<Gt, Vb2>              Other;
  };

  typedef typename Gt::Point_3 Point;
  typedef typename Gt::FT      NT;
  typedef Alpha_status<NT>     Alpha_status;

private:
  Alpha_status* _as;


public:

  Alpha_shape_vertex_base_3()
    : Vb() {}
  
  Alpha_shape_vertex_base_3(const Point& p)
    : Vb(p) {}
  
  Alpha_shape_vertex_base_3(const Point& p, Cell_handle c)
    : Vb(p, c) {}

  Alpha_status*  get_alpha_status() const { return _as;}
  void set_alpha_status(Alpha_status* as) {_as= as;}

};

CGAL_END_NAMESPACE

#endif // CGAL_ALPHA_SHAPE_VERTEX_BASE_3_H
