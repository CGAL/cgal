// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Triangulation_utils_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_UTILS_3_H
#define CGAL_TRIANGULATION_UTILS_3_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

struct Triangulation_utils_3
{
  static const char tab_next_around_edge[4][4];

  static int ccw(int i)
    {
      CGAL_triangulation_precondition( i >= 0 && i < 3 );
      return (i==2) ? 0 : i+1;
    }
  
  static int cw(int i)
    {
      CGAL_triangulation_precondition( i >= 0 && i < 3 );
      return (i==0) ? 2 : i-1;
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

  static unsigned int random_value, count, val;

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
