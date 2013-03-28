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
// Author(s)     : Nico Kruithof

#include <cassert>
#include <iostream>
#include <CGAL/use.h>

template <class Triangulation>
void
_test_cls_triangulation_simplex_3(const Triangulation &)
{
  typedef Triangulation                        Cls;

  typedef typename Cls::Simplex                Simplex;

  typedef typename Cls::Point                  Point;

  typedef typename Cls::Vertex                 Vertex;
  typedef typename Cls::Cell                   Cell;
  typedef typename Cls::Facet                  Facet;
  typedef typename Cls::Edge                   Edge;

  typedef typename Cls::size_type              size_type;
  typedef typename Cls::difference_type        difference_type;

  typedef typename Cls::Vertex_handle          Vertex_handle;
  typedef typename Cls::Cell_handle            Cell_handle;

  typedef typename Cls::Cell_circulator        Cell_circulator;
  typedef typename Cls::Facet_circulator       Facet_circulator;

  typedef typename Cls::Cell_iterator          Cell_iterator;
  typedef typename Cls::Facet_iterator         Facet_iterator;
  typedef typename Cls::Edge_iterator          Edge_iterator;
  typedef typename Cls::Vertex_iterator        Vertex_iterator;

  typedef typename Cls::Finite_vertices_iterator    Finite_vertices_iterator;
  typedef typename Cls::Finite_edges_iterator       Finite_edges_iterator;
  typedef typename Cls::Finite_facets_iterator      Finite_facets_iterator;
  typedef typename Cls::Finite_cells_iterator       Finite_cells_iterator;

  CGAL_USE_TYPE(Vertex);
  CGAL_USE_TYPE(Cell);
  CGAL_USE_TYPE(size_type);
  CGAL_USE_TYPE(difference_type);
  CGAL_USE_TYPE(Cell_circulator);
  CGAL_USE_TYPE(Facet_circulator);
  CGAL_USE_TYPE(Cell_iterator);
  CGAL_USE_TYPE(Facet_iterator);
  CGAL_USE_TYPE(Edge_iterator);
  CGAL_USE_TYPE(Vertex_iterator);
  //########################################################################
  Cls t;

  // Initialise to a 3D triangulation:
  t.insert(Point(0,0,0));
  t.insert(Point(1,0,0));
  t.insert(Point(0,1,0));
  t.insert(Point(0,0,1));

  {  // Check vertices:
    Finite_vertices_iterator vit = t.finite_vertices_begin();
    Vertex_handle vh = vit;
    
    Simplex s1 = vh;

    Simplex s2(vit);
    Simplex s3(vh);

    Simplex s4(s1);
    
    assert(s1.dimension() == 0);
    assert(s1 == s2);
    assert(s1 == s3);
    assert(s1 == s4);
    Vertex_handle vh2 = s1;
    assert(vh == vh2);
  }

  {  // Check edges
    Finite_edges_iterator eit = t.finite_edges_begin();
    Edge e = *eit;
    
    Simplex s1 = *eit;
    Simplex s2 = e;

    Simplex s3(*eit);
    Simplex s4(e);

    Simplex s5(s1);
    
    assert(s1.dimension() == 1);
    assert(s1 == s2);
    assert(s1 == s3);
    assert(s1 == s4);
    assert(s1 == s5);
    Edge e2 = s1;
    assert(e == e2);
  }

  {  // Check facets
    Finite_facets_iterator fit = t.finite_facets_begin();
    Facet f = *fit;
    
    Simplex s1 = *fit;
    Simplex s2 = f;

    Simplex s3(*fit);
    Simplex s4(f);

    Simplex s5(s1);
    
    assert(s1.dimension() == 2);
    assert(s1 == s2);
    assert(s1 == s3);
    assert(s1 == s4);
    assert(s1 == s5);
    Facet f2 = s1;
    assert(f == f2);
  }

  {  // Check cells
    Finite_cells_iterator cit = t.finite_cells_begin();
    Cell_handle ch = cit;
    
    Simplex s1 = Simplex(cit);
    Simplex s2 = ch;

    Simplex s3(cit);
    Simplex s4(ch);

    Simplex s5(s1);
    
    assert(s1.dimension() == 3);
    assert(s1 == s2);
    assert(s1 == s3);
    assert(s1 == s4);
    assert(s1 == s5);
    Cell_handle ch2 = s1;
    assert(ch == ch2);
  }
}
