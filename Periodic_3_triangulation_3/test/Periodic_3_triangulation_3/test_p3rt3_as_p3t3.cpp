// Copyright (c) 2015  INRIA Sophia-Antipolis (France).
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
//
//
// Author(s)     : Aymeric Pelle

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/_test_cls_periodic_3_triangulation_3.h>

#include <CGAL/Periodic_3_regular_triangulation_3.h>
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_triangulation_3.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel         Epeck;
typedef CGAL::Periodic_3_regular_triangulation_traits_3<Epeck>    PRTT_Exact;

typedef CGAL::Exact_predicates_inexact_constructions_kernel       Epick;
typedef CGAL::Periodic_3_regular_triangulation_traits_3<Epick>    PRTT_Inexact;

int main(int, char**)
{
  typedef CGAL::Periodic_3_regular_triangulation_3<PRTT_Exact>    P3RT3_Exact;
  _test_periodic_3_triangulation_3_constructors( P3RT3_Exact() );
  _test_cls_periodic_3_triangulation_3(P3RT3_Exact(),
                                       PRTT_Exact::Weighted_point_3(0.816982, 0.161518, -0.0942375),
                                       "data/P3RT3_covering_test_HOM.tri",
                                       "data/P3RT3_covering_test.tri",
                                       true);

  typedef CGAL::Periodic_3_regular_triangulation_3<PRTT_Inexact>  P3RT3_Inexact;
  _test_periodic_3_triangulation_3_constructors( P3RT3_Inexact() );
  _test_cls_periodic_3_triangulation_3(P3RT3_Inexact(),
                                       PRTT_Inexact::Weighted_point_3(0.816982, 0.161518, -0.0942375),
                                       "data/P3RT3_covering_test_HOM.tri",
                                       "data/P3RT3_covering_test.tri");

  return 0;
}
