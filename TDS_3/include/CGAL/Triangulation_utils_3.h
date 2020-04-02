// Copyright (c) 1999  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>

#ifndef CGAL_TRIANGULATION_UTILS_3_H
#define CGAL_TRIANGULATION_UTILS_3_H

#include <CGAL/license/TDS_3.h>


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

};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_UTILS_3_H
