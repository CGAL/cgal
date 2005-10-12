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

#ifndef CGAL_KDS_DELAUNAY_TRIANGULATION_2_WATCHER_BASE_H
#define CGAL_KDS_DELAUNAY_TRIANGULATION_2_WATCHER_BASE_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_NAMESPACE

struct Delaunay_triangulation_visitor_base_2 {
  Delaunay_triangulation_visitor_base_2(){}
  template <class Ok>
  void delete_vertex(Ok) {
  }

  template <class Ok>
  void new_vertex(Ok) {
  }

  template <class Ok>
  void change_vertex(Ok) {
  }

  template <class It>
  void new_faces(It, It) {
  }

  template <class It>
  void delete_faces(It, It) {
  }

  template <class E>
  void pre_flip(E){

  }
  template <class E>
  void post_flip(E){

  }
};

CGAL_KDS_END_NAMESPACE

#endif
