// Copyright (c) 2003
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

// This defines removes the operator/ from CGAL::Mpzf to check that functors not using
// the tag `Needs_FT<>` really only need a RT (ring type) without division.
#define CGAL_NO_MPZF_DIVISION_OPERATOR 1

#include <CGAL/Cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <cassert>

#include "CGAL/Precise_numbers.h"

#define TEST_FILENAME "Test-Filtered_cartesian-IO.out"
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

template <typename Cls>
void test();
int
main()
{
  CGAL::Set_ieee_double_precision double_precision_guard;
  typedef CGAL::Cartesian<double>                               Clsdb;
  typedef CGAL::Filtered_kernel<Clsdb>                          Clsd;

  // typedef CGAL::Cartesian<CGAL::Quotient<Precise_integer> >     Clsb;
  // typedef CGAL::Cartesian<CGAL::Quotient<CGAL::MP_Float> >      Clsb;
  // typedef CGAL::Filtered_kernel<Clsb>                           Cls;
  typedef CGAL::Exact_predicates_exact_constructions_kernel       Cls;

  std::cout <<
   // "Testing with Filtered_kernel<Cartesian<Quotient<Precise_integer>>>:"
   // "Testing with Filtered_kernel<Cartesian<Quotient<MP_Float>>>:"
   "Testing with Exact_predicates_exact_constructions_kernel:"
            << std::endl;
  std::cout << "Testing IO with F_k<Cartesian<double>>:" << std::endl;
  _test_io( Clsd() );

  std::cout << "Testing with Epeck:\n";
  test<Cls>();
  std::cout << "Testing with Double_precision_epick:\n";
  test<CGAL::Double_precision_epick>();

#  if defined(BOOST_MSVC)
#    pragma warning(push)
#    pragma warning(disable: 4244)
#  endif
  std::cout << "Testing with Simple_precision_epick:\n";
  test<CGAL::Single_precision_epick>();
#  if defined(BOOST_MSVC)
#    pragma warning(pop)
#  endif

  return 0;
}

template <typename Cls>
void test() {
  std::cout << "Testing 2d :";
  std::cout << std::endl;
  _test_2( Cls() );

  std::cout << "Testing 3d :";
  std::cout << std::endl;
  _test_3( Cls() );

  std::cout << "Testing new 2d :";
  std::cout << std::endl;
  test_new_2( Cls() );

  std::cout << "Testing new 3d :";
  std::cout << std::endl;
  test_new_3( Cls() );

  std::cout << "Testing new parts :";
  std::cout << std::endl;
  _test_orientation_and_bounded_side( Cls() );
  _test_fct_points_implicit_sphere( Cls() );
  _test_fct_constructions_2( Cls() );
  _test_fct_constructions_3( Cls() );
  _test_fct_point_3( Cls() );
  _test_fct_coplanar_3( Cls() );
  _test_cls_iso_cuboid_3( Cls() );
  _test_angle( Cls() );

  std::cout << "Testing 3d-2d :";
  std::cout << std::endl;
  _test_mf_plane_3_to_2d( Cls() );
}
