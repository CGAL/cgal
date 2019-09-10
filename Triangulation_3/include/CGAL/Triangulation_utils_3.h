// Copyright (c) 1999  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>

#ifndef CGAL_TRIANGULATION_UTILS_3_H
#define CGAL_TRIANGULATION_UTILS_3_H

#include <CGAL/license/Triangulation_3.h>


#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>

namespace CGAL {

// We use the following template class in order to avoid having a static data
 // member of a non-template class which would require src/Triangulation_3.C .
template < class T = void >
struct Triangulation_utils_base_3
{
  static const char tab_next_around_edge[4][4];
  static const int tab_vertex_triple_index[4][3];

  static const int index_increment_map[4];
  static const int index_jump_map[4];
  static const int index_decrement_map[4];
  
  // copied from Triangulation_utils_2.h to avoid package dependency
  static const int ccw_map[3];
  static const int cw_map[3];
};

template < class T >
const char Triangulation_utils_base_3<T>::tab_next_around_edge[4][4] = {
      {5, 2, 3, 1},
      {3, 5, 0, 2},
      {1, 3, 5, 0},
      {2, 0, 1, 5} };

template < class T >
const int Triangulation_utils_base_3<T>::tab_vertex_triple_index[4][3] = {
 {1, 3, 2},
 {0, 2, 3},
 {0, 3, 1},
 {0, 1, 2}
};

template < class T >
const int Triangulation_utils_base_3<T>::index_increment_map[4] = { 1, 2, 3, 0 };

template < class T >
const int Triangulation_utils_base_3<T>::index_jump_map[4] = { 2, 3, 0, 1 };

template < class T >
const int Triangulation_utils_base_3<T>::index_decrement_map[4] = { 3, 0, 1, 2 };

template < class T >
const int Triangulation_utils_base_3<T>::ccw_map[3] = {1, 2, 0};

template < class T >
const int Triangulation_utils_base_3<T>::cw_map[3] = {2, 0, 1};

// We derive from Triangulation_cw_ccw_2 because we still use cw() and ccw()
// in the 2D part of the code.  Ideally, this should go away when we re-use
// T2D entirely.

struct Triangulation_utils_3
  : public Triangulation_utils_base_3<>
{
  static int ccw(const int i) 
    {
      CGAL_triangulation_precondition( i >= 0 && i < 3);
      return ccw_map[i];
    }

  static int cw(const int i)
    {
      CGAL_triangulation_precondition( i >= 0 && i < 3);
      return cw_map[i];
    }

  static int next_around_edge(const int i, const int j)
  {
    // index of the next cell when turning around the
    // oriented edge vertex(i) vertex(j) in 3d
    CGAL_triangulation_precondition( ( i >= 0 && i < 4 ) &&
		                     ( j >= 0 && j < 4 ) &&
		                     ( i != j ) );
    return tab_next_around_edge[i][j];
  }


  static int vertex_triple_index(const int i, const int j)
  {
    // indexes of the  jth vertex  of the facet of a cell
    // opposite to vertx i
      CGAL_triangulation_precondition( ( i >= 0 && i < 4 ) &&
		                     ( j >= 0 && j < 3 ) );
    return tab_vertex_triple_index[i][j];
  }

  // Get the index of the next vertex or facet.
  static int increment_index( int li ) {
      CGAL_triangulation_precondition( li >= 0 && li < 4 );
      return index_increment_map[ li ];
  }

  // Get the index of the vertex or facet two places further.
  static int jump_index( int li ) {
      CGAL_triangulation_precondition( li >= 0 && li < 4 );
      return index_jump_map[ li ];
  }

  // Get the index of the previous vertex or facet.
  static int decrement_index( int li ) {
      CGAL_triangulation_precondition( li >= 0 && li < 4 );
      return index_decrement_map[ li ];
  }
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_UTILS_3_H
