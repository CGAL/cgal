#define CGAL_PROFILE

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
// Author(s)     : Andreas Fabri

#include <CGAL/internal/Exact_type_selector.h>
#include <CGAL/Filtered_rational_kernel.h>

#define TEST_FILENAME "Test-Filtered_rational-IO.out"
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

#include <cassert>

#include <CGAL/gmpxx.h>
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#include <CGAL/CORE_Expr.h>

#if 1
typedef CORE::Expr ET;
#elif 0
typedef CGAL::internal::Exact_field_selector<void*> ET;
typedef leda_rational ET;
typedef CGAL::Gmpq ET;
#endif

typedef CGAL::Simple_cartesian<CGAL::Interval_nt<false> > AK;
typedef CGAL::Simple_cartesian<ET>                        EK;
typedef CGAL::Filtered_rational_kernel<AK, EK>            Cls;

int main()
{
  std::cout.precision(17);
  std::cerr.precision(17);

  std::cout << "ET: " << typeid(ET).name() << std::endl;
  std::cout << "Category: " << typeid(typename CGAL::Algebraic_structure_traits<ET>::Algebraic_category).name() << std::endl;
  std::cout << "FT: " << typeid(Cls::FT).name() << std::endl;
  std::cout << "Category: " << typeid(typename CGAL::Algebraic_structure_traits<Cls::FT>::Algebraic_category).name() << std::endl;

  std::cout << "Testing 2d with Filtered_rational_kernel :";
  std::cout << std::endl;
  _test_2( Cls() );

  std::cout << "Testing 3d with Filtered_rational_kernel :";
  std::cout << std::endl;
  _test_3( Cls() );

  std::cout << "Testing new 2d with Filtered_rational_kernel :";
  std::cout << std::endl;
  test_new_2( Cls() );
  std::cout << "Testing new 3d with Filtered_rational_kernel :";
  std::cout << std::endl;
  test_new_3( Cls() );

  std::cout << "Testing new parts with Filtered_rational_kernel :";
  std::cout << std::endl;
  _test_orientation_and_bounded_side( Cls() );
  _test_fct_points_implicit_sphere( Cls() );
  _test_fct_constructions_2( Cls() );
  _test_fct_constructions_3( Cls() );
  _test_fct_point_3( Cls() );
  _test_fct_coplanar_3( Cls() );
  _test_cls_iso_cuboid_3( Cls() );
  _test_angle( Cls() );

  std::cout << "Testing 3d-2d with Filtered_rational_kernel :";
  std::cout << std::endl;
  _test_mf_plane_3_to_2d( Cls() );

  std::cout << "Done!" << std::endl;

  return 0;
}
