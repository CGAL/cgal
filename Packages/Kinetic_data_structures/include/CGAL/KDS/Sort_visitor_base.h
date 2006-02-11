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

#ifndef CGAL_KDS_DELAUNAY_SORT_WATCHER_BASE_H
#define CGAL_KDS_DELAUNAY_SORT_WATCHER_BASE_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_NAMESPACE

struct Sort_visitor_base
{
  Sort_visitor_base(){}
  template <class Vertex_handle>
  void remove_vertex(Vertex_handle) {
  }

  template <class Vertex_handle>
  void create_vertex(Vertex_handle) {
  }

  template <class Vertex_handle>
  void modify_vertex(Vertex_handle) {
  }

  template <class Vertex_handle>
  void before_swap(Vertex_handle, Vertex_handle) {

  }
  template <class Vertex_handle>
  void after_swap(Vertex_handle, Vertex_handle) {

  }
};

CGAL_KDS_END_NAMESPACE
#endif
