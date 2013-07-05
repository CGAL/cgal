// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_DELAUNAY_TRIANGULATION_2_WATCHER_BASE_H
#define CGAL_KINETIC_DELAUNAY_TRIANGULATION_2_WATCHER_BASE_H
#include <CGAL/Kinetic/basic.h>

namespace CGAL { namespace Kinetic {

struct Delaunay_triangulation_visitor_base_2
{
    Delaunay_triangulation_visitor_base_2(){}
  
  template <class Point_key>
  void pre_insert_vertex(Point_key ) {
    }

    template <class Vertex_handle>
    void post_insert_vertex(Vertex_handle ) {
    }

    template <class Vertex_handle>
    void pre_remove_vertex(Vertex_handle ) {
    }


    template <class Point_key>
    void post_remove_vertex(Point_key ) {
    }

    template <class Vertex_handle>
    void change_vertex(Vertex_handle ) {
    }

    template <class Face_handle> 
    void create_face(Face_handle ) {
    }

    template <class Face_handle>
    void destroy_face(Face_handle ) {
    }


    template <class E>
    void pre_flip(E) {

    }
    template <class E>
    void post_flip(E) {

    }

  template <class P>
  bool filter_certificate(P , P , P ) const {
    return true;
  }
  template <class P>
  bool filter_certificate(P , P , P , P ) const {
    return true;
  }
};

} } //namespace CGAL::Kinetic
#endif
