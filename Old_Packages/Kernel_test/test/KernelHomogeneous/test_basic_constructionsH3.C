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
// file          : test_basic_constructionsH3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>
#include <cassert>

#include <CGAL/Precise_numbers.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Point_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/basic_constructionsH3.h>
#include <CGAL/intersection_3.h>

int
main()
{
  typedef   CGAL::Homogeneous<Precise_integer>             H_Cls;
  typedef   CGAL::Point_3< H_Cls >                         Point;
  typedef   CGAL::Line_3< H_Cls >                          Line;
  typedef   CGAL::Plane_3< H_Cls >                         Plane;

  Point ep1( 8, 2, 4, 1);
  Point ep2( 9, 3,18, 3);
  Point ep3(12,24,-6, 2);
  Point p( 17, 12, 4, 1);

  Plane pl( ep1, ep2, ep3);
  CGAL::Object o = CGAL::intersection( pl, Line( p, pl.orthogonal_direction()));
  Point ip;
  assert( CGAL::assign(ip,o) );
  assert( pl.has_on(ip) );
  Point pp = CGAL::_projection(p, pl);
  assert( pl.has_on(pp) );
  assert( pp == ip );

  return 0;
}
