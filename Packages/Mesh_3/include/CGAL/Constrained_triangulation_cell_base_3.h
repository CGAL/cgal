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

#include <CGAL/Triangulation_cell_base_3.h>

CGAL_BEGIN_NAMESPACE

template <class Gt, class Cb = Triangulation_cell_base_3<Gt> >
class Constrained_triangulation_cell_base_3 : public Cb
{
public:
  typedef typename Cb::Vertex_handle  Vertex_handle;
  typedef typename Cb::Cell_handle    Cell_handle;

  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other  Cb2;
    typedef Constrained_triangulation_cell_base_3<Gt, Cb2> Other;
  };

  Constrained_triangulation_cell_base_3()
    : Cb() {}

  Constrained_triangulation_cell_base_3(Vertex_handle v0, Vertex_handle v1,
					Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3) {}

  Constrained_triangulation_cell_base_3(Vertex_handle v0, Vertex_handle v1,
					Vertex_handle v2, Vertex_handle v3,
					Cell_handle   n0, Cell_handle   n1,
					Cell_handle   n2, Cell_handle   n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3) {}

protected:
  bool C[4];

public:
  bool is_constrained(int i) const
    { 
      CGAL_triangulation_precondition( i == 0 || i == 1 ||
				       i == 2 || i ==3 );
      return(C[i]);
    }

  void set_constrained(int i, bool b)
    {
      CGAL_triangulation_precondition( i == 0 || i == 1 ||
				       i == 2 || i ==3 );
      C[i] = b;
    }
};

CGAL_END_NAMESPACE
