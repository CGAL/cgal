// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$ $Date$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_DELAUNAY_SORT_WATCHER_BASE_H
#define CGAL_KINETIC_DELAUNAY_SORT_WATCHER_BASE_H
#include <CGAL/Kinetic/basic.h>

CGAL_KINETIC_BEGIN_NAMESPACE

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

CGAL_KINETIC_END_NAMESPACE
#endif
