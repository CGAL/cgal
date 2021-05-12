// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// $Date$
//
//
// Author(s)     : Aymeric PELLE <aymeric.pelle@sophia.inria.fr>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Timer.h>

#include <CGAL/_test_cls_periodic_3_regular_triangulation_traits_3.h>

#include <iostream>

void test_periodic_3_regular_triangulation_traits_3()
{
  CGAL::Timer t;
  t.start();
  std::cout << "EPECK" << std::endl;
  _test_cls_periodic_3_regular_triangulation_traits_3_rational<CGAL::Exact_predicates_exact_constructions_kernel>();
  _test_cls_periodic_3_regular_triangulation_traits_3_irrational<CGAL::Exact_predicates_exact_constructions_kernel>();

  std::cout << "EPICK" << std::endl;
  _test_cls_periodic_3_regular_triangulation_traits_3_rational<CGAL::Exact_predicates_inexact_constructions_kernel>();

  std::cout << t.time() << " sec." << std::endl;
}
