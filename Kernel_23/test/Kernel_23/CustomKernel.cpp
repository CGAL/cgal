// Copyright (c) 2001,2002
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
//
//
// Author(s)     : Sylvain Pion


#include <CGAL/CustomKernel.h>
#include <cassert>

#include "CGAL/_test_cls_kernel.h"
#define TEST_FILENAME "Test-Cartesian-IO.out"
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
#include "CGAL/_test_cls_circle_3.h"

#include "CGAL/_test_mf_plane_3_to_2d.h"

#include <vector>
#include <CGAL/centroid.h>
#include <CGAL/bounding_box.h>

int
main()
{
  CGAL::Set_ieee_double_precision double_precision_guard;
  typedef   CGAL::Test::CustomKernel     Cls;
  typedef Cls::Point_3 Point_3;
  typedef Cls::Segment_3 Segment_3;
  Point_3 p(1,2,3);
  Point_3 q(4,5,6);
  std::cout << p  << std::endl;
  p.bbox();
  p.rep().x;
  Segment_3 seg(p,q);
  std::cout << seg.bbox() << std::endl;

  CGAL::Test::Vec3 vp = p;
  CGAL::Test::Vec3 vq = p;

  vp = p;
  p = vp;

  p.foo();
  std::cout << CGAL::midpoint(p,q) << std::endl;
  std::cout << CGAL::midpoint<Cls>(vp,vq) << std::endl;
  std::cout << Cls::Construct_midpoint_3()(vp,vq) << std::endl;
  std::vector<CGAL::Test::Vec3> vvec3;

  CGAL::bounding_box(vvec3.begin(), vvec3.end(), Cls());

  // CGAL::centroid(vvec3.begin(), vvec3.end(), Cls(), CGAL::Dimension_tag<0>());
#if 1
  Cls::Construct_point_3 cp3 = Cls().construct_point_3_object();
  p = cp3(4,5,6);
  p = cp3(CGAL::ORIGIN);

  std::cout << "Testing nested types with CustomKernel :" << std::endl;
  _test_kernel( Cls() );

  //std::cout << "Testing IO with CustomKernel  :" << std::endl;
  //_test_io( Cls() );


  std::cout << "Testing 2d with CustomKernel  :";
  std::cout << std::endl;
  _test_2( Cls() );

  std::cout << "Testing 3d with CustomKernel  :";
  std::cout << std::endl;
  _test_3( Cls() );

  std::cout << "Testing new 2d with CustomKernel  :";
  std::cout << std::endl;
  test_new_2( Cls() );

  std::cout << "Testing new 3d with CCustomKernel  :";
  std::cout << std::endl;
  test_new_3( Cls() );
  std::cout << "Testing new parts with CustomKernel :";
  std::cout << std::endl;
  _test_orientation_and_bounded_side( Cls() );
  _test_fct_points_implicit_sphere( Cls() );
  _test_fct_constructions_2( Cls() );
  _test_fct_constructions_3( Cls() );
  _test_fct_point_3( Cls() );
  _test_fct_coplanar_3( Cls() );
  _test_cls_iso_cuboid_3( Cls() );
  _test_angle( Cls() );

#elif 1
  std::cout << "Testing 3d-2d with CustomKernel  :";
  std::cout << std::endl;
  _test_mf_plane_3_to_2d( Cls() );
#endif
  std::cout << "All tests done" << std::endl;
  return 0;
}



