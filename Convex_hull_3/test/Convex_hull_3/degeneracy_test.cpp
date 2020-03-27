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
// file          : test/Convex_hull_3/degeneracy_test.C
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Convex_hull_2/ch_assertions.h>
#include <CGAL/convex_hull_3.h>
#include <vector>
#include <fstream>

typedef double                                     NumberType;
typedef CGAL::Simple_cartesian<NumberType>         Kernel;
typedef CGAL::Polyhedron_3<Kernel>                 Polyhedron;
typedef Kernel::Point_3                            Point_3;

int main()
{
  /* read points from file */
  std::ifstream F("Point_3_list.txt");
  std::vector<Point_3> V;
  std::copy( std::istream_iterator<Point_3>(F),
             std::istream_iterator<Point_3>(),
             std::back_inserter(V));

  Polyhedron P; /* define polyhedron to hold convex hull */

  /* compute convex hull */
  CGAL::convex_hull_3( V.begin(), V.end(), P);
  return 0;
}
