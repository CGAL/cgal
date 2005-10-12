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

#ifndef CGAL_KDS_DELAUNAY_TRIANGULATION_3_WATCHER_BASE_H
#define CGAL_KDS_DELAUNAY_TRIANGULATION_3_WATCHER_BASE_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_NAMESPACE

struct Delaunay_triangulation_visitor_base_3 {
  //typedef Tr Triangulation;
  Delaunay_triangulation_visitor_base_3(){}

  template <class Vertex_handle>
  void delete_vertex(Vertex_handle) {
  }

  template <class Vertex_handle>
  void new_vertex(Vertex_handle) {
  }

  template <class Vertex_handle>
  void change_vertex(Vertex_handle) {
  }

  template <class It>
  void new_cells(It, It) {
  }

  template <class It>
  void delete_cells(It, It) {
  }

  template <class Edge>
  void pre_edge_flip(Edge){

  }
  template <class Edge>
  void post_facet_flip(Edge){

  }

  template <class Facet>
  void pre_facet_flip(Facet){

  }
  
  template <class Facet>
  void post_edge_flip(Facet){
  }
};

CGAL_KDS_END_NAMESPACE

#endif
