// Copyright (c) 1999,2002  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra, Sylvain Pion


#include <CGAL/Homogeneous.h>
#include <cassert>
#include <CGAL/intersection_3.h>

#include "CGAL/Precise_numbers.h"
#define TEST_FILENAME "Test-Homogeneous-IO.out"
#include "CGAL/_test_io.h"
#include "CGAL/_test_2.h"
#include "CGAL/_test_3.h"

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
}
 

int
main()
{
  test_basic();

  typedef   CGAL::Homogeneous<double>     Clsd;
  std::cout << "Testing IO with Homogeneous<double> :" << std::endl;
  _test_io( Clsd() );

  typedef   CGAL::Homogeneous<Precise_integer>     Cls;
  std::cout << "Testing 2d with Homogeneous<Precise_integer> :" << std::endl;
  _test_2( Cls() );

  std::cout << "Testing 3d with Homogeneous<Precise_integer> :" << std::endl;
  _test_3( Cls() );

  test_new_2( Cls() );
  test_new_3( Cls() );

  std::cout << "Testing new parts with Homogeneous<Precise_integer> :";
  std::cout << std::endl;
  _test_orientation_and_bounded_side( Cls() );
  _test_fct_points_implicit_sphere( Cls() );
  _test_fct_constructions_2( Cls() );
  _test_fct_constructions_3( Cls() );
  _test_fct_point_3( Cls() );
  _test_fct_coplanar_3( Cls() );
  _test_cls_iso_cuboid_3( Cls() );
  _test_angle( Cls() );

  std::cout << "Testing 3d-2d with Homogeneous<Precise_integer> :";
  std::cout << std::endl;
  _test_mf_plane_3_to_2d( Cls() );
  return 0;
}
