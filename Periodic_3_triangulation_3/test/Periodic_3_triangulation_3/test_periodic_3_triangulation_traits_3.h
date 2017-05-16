// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// $Date$
//
//
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>

#include <cassert>

#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Lazy_exact_nt.h>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#endif
#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_real.h>
#endif
#ifdef CGAL_USE_CORE
#include <CGAL/CORE_Expr.h>
#endif


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Homogeneous.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/_test_periodic_3_static_filters.h>
#include <CGAL/_test_cls_periodic_3_triangulation_traits_3.h>

void test_periodic_3_triangulation_traits_3()
{
  std::cout<<"Statically filtered predicates:"<<std::endl;
  _test_periodic_3_static_filters();

#define test_traits _test_cls_periodic_3_triangulation_traits_3
  std::cout<<"\nTesting predicates:"<<std::endl;
  std::cout<<"  Predefined kernels...";std::cout.flush();
  test_traits<CGAL::Exact_predicates_inexact_constructions_kernel>();
  test_traits<CGAL::Exact_predicates_exact_constructions_kernel>();
  std::cout<<" done"<<std::endl;

  std::cout<<"  MP_Float...";std::cout.flush();

#define RT CGAL::MP_Float
#define FT CGAL::Quotient<RT>
  test_traits<KERNEL>();
#define LRT CGAL::Lazy_exact_nt< RT >
#define LFT CGAL::Lazy_exact_nt< FT >
  test_traits<LAZY_KERNEL>();
  std::cout<<" done"<<std::endl;
#undef RT
#undef FT
#undef LRT
#undef LFT

#ifdef CGAL_USE_GMP
  std::cout<<"  GMP...";std::cout.flush();
#define RT CGAL::Gmpz
#define FT CGAL::Gmpq
  test_traits<KERNEL>();
#define LRT CGAL::Lazy_exact_nt< RT >
#define LFT CGAL::Lazy_exact_nt< FT >
  test_traits<LAZY_KERNEL>();
#undef RT
#undef FT
#undef LRT
#undef LFT
#else
  std::cout<<"  GMP not available"<<std::endl;
#endif

#ifdef TEST_CARTESIAN
#ifdef CGAL_USE_LEDA
  std::cout<<"  LEDA...";std::cout.flush();
#define RT leda_integer
#define FT leda_rational
  test_traits<KERNEL>();

  std::cout<<" done"<<std::endl;
#undef RT
#undef FT
#else
  std::cout<<"  LEDA not available"<<std::endl;
#endif

#ifdef CGAL_USE_CORE
  std::cout<<"  CORE...";std::cout.flush();
#define RT CORE::Expr
#define FT CGAL::Quotient<RT>
  test_traits<KERNEL>();

  std::cout<<" done"<<std::endl;
#undef RT
#undef FT
#else
  std::cout<<"  CORE not available"<<std::endl;
#endif

#endif
}
