// ============================================================================
//
// Copyright (c) 1999,2002 The CGAL Consortium
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
// file          : test/Kernel/Simple_homogeneous.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra, Sylvain Pion
//
//
// coordinator   : MPI, Saarbruecken
// ============================================================================
 


#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Quotient.h>
#include <cassert>
#include <CGAL/intersection_3.h>

#include "CGAL/Precise_numbers.h"
#include "CGAL/_test_io.h"
#include "CGAL/_test_2.C"
#include "CGAL/_test_3.C"

#include "CGAL/_test_new_2.h"
#include "CGAL/_test_new_3.h"

#include "CGAL/_test_fct_points_implicit_sphere.h"
#include "CGAL/_test_orientation_and_bounded_side.h"
#include "CGAL/_test_fct_constructions_2.h"
#include "CGAL/_test_fct_constructions_3.h"
#include "CGAL/_test_fct_point_3.h"
#include "CGAL/_test_fct_coplanar_3.h"
#include "CGAL/_test_cls_iso_cuboid_3.h"
#include "CGAL/_test_angle.h"

#include "CGAL/_test_mf_plane_3_to_2d.h"

// This one should be merged with the global test-suite.
void
test_basic()
{
  typedef   CGAL::Simple_homogeneous<Precise_integer>             H_Cls;
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
}
 

int
main()
{
  test_basic();

  typedef   CGAL::Simple_homogeneous<double>     Clsd;
  std::cout << "Testing IO with Simple_homogeneous<double> :" << std::endl;
  _test_io( Clsd() );

  typedef   CGAL::Simple_homogeneous<Precise_integer>     Cls;
  std::cout << "Testing 2d with Simple_homogeneous<Precise_integer> :" << std::endl;
  _test_2( Cls() );

  std::cout << "Testing 3d with Simple_homogeneous<Precise_integer> :" << std::endl;
  _test_3( Cls() );

  test_new_2( Cls() );
  test_new_3( Cls() );

  std::cout << "Testing new parts with Simple_homogeneous<Precise_integer> :";
  std::cout << std::endl;
  _test_orientation_and_bounded_side( Cls() );
  _test_fct_points_implicit_sphere( Cls() );
  _test_fct_constructions_2( Cls() );
  _test_fct_constructions_3( Cls() );
  _test_fct_point_3( Cls() );
  _test_fct_coplanar_3( Cls() );
  _test_cls_iso_cuboid_3( Cls() );
  _test_angle( Cls() );

  std::cout << "Testing 3d-2d with Simple_homogeneous<Precise_integer> :";
  std::cout << std::endl;
  _test_mf_plane_3_to_2d( Cls() );
  return 0;
}
