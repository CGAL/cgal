// Copyright (c) 2005  Stanford University (USA).
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_REGULAR_TRIANGULATION_3_WATCHER_BASE_H
#define CGAL_KDS_REGULAR_TRIANGULATION_3_WATCHER_BASE_H
#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/Delaunay_triangulation_visitor_base_3.h>

CGAL_KDS_BEGIN_NAMESPACE

struct Regular_triangulation_visitor_base_3: 
  public Delaunay_triangulation_visitor_base_3 {
  //typedef Tr Triangulation;
  Regular_triangulation_visitor_base_3(){}

  template <class Key, class Cell>
  void pre_move(Key, Cell){}

  template <class Key, class Cell>
  void post_move(Key, Cell){}
  
  template <class Key, class Cell>
  void pre_push(Key, Cell){}

  template <class Vertex_handle>
  void post_push(Vertex_handle){}

  template <class Vertex_handle>
  void pre_pop(Vertex_handle){}

  template <class Key, class Cell>
  void post_pop(Key, Cell){}
};

CGAL_KDS_END_NAMESPACE

#endif
