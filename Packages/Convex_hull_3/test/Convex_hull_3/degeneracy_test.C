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
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  
// ============================================================================
 

#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/ch_assertions.h>
#include <CGAL/convex_hull_3.h>
#include <vector>
#include <fstream>

typedef double                                                      NumberType;

/* representation class */
typedef CGAL::Cartesian<NumberType>                                 RepCls;

/* define polyhedron type */
typedef CGAL::Halfedge_data_structure_polyhedron_default_3<RepCls>  HDS;
typedef CGAL::Polyhedron_default_traits_3<RepCls>                   PolyTraits;
typedef CGAL::Polyhedron_3< PolyTraits, HDS>                        Polyhedron;
typedef CGAL::Point_3<RepCls>                                       Point;

int
main()
{
  /* read points from file */
  std::ifstream F("Point_3_list.txt");
  std::vector<Point> V;
  std::copy( std::istream_iterator<Point>(F),
             std::istream_iterator<Point>(),
             std::back_inserter(V));

  /* define polyhedron to hold convex hull */
  Polyhedron P;

  /* compute convex hull */
  CGAL::convex_hull_3( V.begin(), V.end(), P);
  return 0;
}
