// Copyright (c) 1999  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>

#ifndef CGAL_TRIANGULATION_UTILS_3_H
#define CGAL_TRIANGULATION_UTILS_3_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Triangulation_utils_2.h>

CGAL_BEGIN_NAMESPACE

// We use the following template class in order to avoid having a static data
// member of a non-template class which would require src/Triangulation_3.C .
template < class T = void >
struct Triangulation_utils_base_3
{
  static const char tab_next_around_edge[4][4];

  static unsigned int random_value, count, val;
};

template < class T >
const char Triangulation_utils_base_3<T>::tab_next_around_edge[4][4] = {
      {5, 2, 3, 1},
      {3, 5, 0, 2},
      {1, 3, 5, 0},
      {2, 0, 1, 5} };

template < class T >
unsigned int Triangulation_utils_base_3<T>::random_value = 0;

template < class T >
unsigned int Triangulation_utils_base_3<T>::count = 0;

template < class T >
unsigned int Triangulation_utils_base_3<T>::val;


// We derive from Triangulation_cw_ccw_2 because we still use cw() and ccw()
// in the 2D part of the code.  Ideally, this should go away when we re-use
// T2D entirely.

struct Triangulation_utils_3
  : public Triangulation_cw_ccw_2,
    public Triangulation_utils_base_3<>
{
  static int next_around_edge(const int i, const int j)
  {
    // index of the next cell when turning around the
    // oriented edge vertex(i) vertex(j) in 3d
    CGAL_triangulation_precondition( ( i >= 0 && i < 4 ) &&
		                     ( j >= 0 && j < 4 ) &&
		                     ( i != j ) );
    return tab_next_around_edge[i][j];
  }

  // rand_4() outputs pseudo random unsigned ints < 4.
  // We compute random 16 bit values, that we slice/shift to make it faster.
  static unsigned int rand_4()
  {
      if (count==0)
      {
          count = 16;
          random_value = (421 * random_value + 2073) % 32749;
          val = random_value;
      }
      count--;
      unsigned int ret = val & 3;
      val = val >> 1;
      return ret;
  }

  static unsigned int rand_3()
  {
      unsigned int i = rand_4();
      return i==3 ? 0 : i;
  }
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_UTILS_3_H
