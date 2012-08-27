// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
//
// Author(s)     : Rebufat Francois

#include <cassert>

template <class Vertex>
void
_test_vertex_tds_3(const Vertex &)
{
  typedef typename Vertex::Triangulation_data_structure  Tds;
  typedef typename Tds::Cell_handle                      Cell_handle;
  typedef typename Tds::Vertex_handle                    Vertex_handle;

  Tds tds;

  Cell_handle c1 = tds.create_cell();
  Vertex_handle v1 = tds.create_vertex();
  v1->set_cell(c1);
  assert(v1->cell() == c1);
  c1->set_vertex(0, v1);
  assert(tds.is_valid(v1));

  Cell_handle c2 = tds.create_cell();
  v1->set_cell(c2);
  c2->set_vertex(0, v1);
  assert(tds.is_valid(v1));

  // Unicity of the default constructed handle.
  Vertex_handle v = Vertex_handle();
  v = Vertex_handle();
  assert(v == Vertex_handle());
  assert(v1 != Vertex_handle());

  Cell_handle c = Cell_handle();
  c = Cell_handle();
  assert(c == Cell_handle());
  assert(c1 != Cell_handle());

  // We want the following comparisons to work for use in std::set<>...
  bool b1 = v<v1;  (void) b1;
  bool b2 = c<c1;  (void) b2;
}
