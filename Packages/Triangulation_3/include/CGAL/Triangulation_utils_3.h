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
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_UTILS_3_H
#define CGAL_TRIANGULATION_UTILS_3_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

struct Triangulation_utils_3
{
  static const char tab_next_around_edge[4][4];

  int ccw(const int i) const
    {
      CGAL_triangulation_precondition( 3 > (unsigned) i );
      return (i==2) ? 0 : i+1;
    }
  
  int cw(const int i) const
    {
      CGAL_triangulation_precondition( 3 > (unsigned) i );
      return (i==0) ? 2 : i-1;
    }

  int next_around_edge(const int i, const int j) const
  {
    // index of the next cell when turning around the
    // oriented edge vertex(i) vertex(j) in 3d
    CGAL_triangulation_precondition( ( 4 > (unsigned) i) &&
		                     ( 4 > (unsigned) j) &&
		                     ( i != j ) );
    return tab_next_around_edge[i][j];
  }

  unsigned int rand_4() const
  {
      static unsigned int random_value = 0;
      static unsigned int count = 0;
      static unsigned int val;
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

  unsigned int rand_3() const
  {
      unsigned int i = rand_4();
      return i==3 ? 0 : i;
  }
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_UTILS_3_H
