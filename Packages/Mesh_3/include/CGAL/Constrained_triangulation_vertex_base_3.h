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

#include <CGAL/Triangulation_vertex_base_3.h>
#include <set>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Vb = CGAL::Triangulation_vertex_base_3<Gt> >
class Constrained_triangulation_vertex_base_3
  : public Vb
{
public:
  typedef typename Vb::Vertex_handle  Vertex_handle;
  typedef typename Vb::Cell_handle    Cell_handle;
  typedef typename Vb::Point          Point;

  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
    typedef Constrained_triangulation_vertex_base_3<Gt, Vb2> Other;
  };

  Constrained_triangulation_vertex_base_3()
    : Vb() {}

  Constrained_triangulation_vertex_base_3(const Point & p)
    : Vb(p) {}

  Constrained_triangulation_vertex_base_3(const Point & p, Cell_handle c)
    : Vb(p, c) {}

protected:
  int surface_number;

public:
  bool is_adjacent_by_constraint(const Vertex_handle v)
    {
      return form_constraint_with_vertices.find(v) !=
	form_constraint_with_vertices.end();
    }

  void set_is_adjacent_by_constraint(const Vertex_handle v,
				     const bool b)
    { 
      if (b) form_constraint_with_vertices.insert(v);
      else form_constraint_with_vertices.erase(v);
    }
};
CGAL_END_NAMESPACE
