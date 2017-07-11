// Copyright (c) 1998, 2015  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Francois Rebufat
//                 Manuel Caroli
//                 Aymeric Pelle

#include <iostream>

#include <CGAL/Timer.h>

#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_triangulation_3.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel     K1;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K1>    PDTT1;

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel       K2;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K2>    PDTT2;

#include <CGAL/MP_Float.h>
#include <CGAL/Simple_homogeneous.h>
typedef CGAL::Simple_homogeneous<CGAL::MP_Float>                K3;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K3>    PDTT3;
// Explicit instantiation of the whole class:
template class CGAL::Periodic_3_triangulation_3<PDTT3>;

#include <CGAL/_test_cls_periodic_3_triangulation_3.h>

// We use Periodic_3_Delaunay_triangulation_3 because
// Periodic_3_triangulation_3 does not have an insert function and
// thus we cannot construct non-trivial triangulations without using
// the insert from the periodic Delaunay triangulation.

int main(int, char**)
{
  CGAL::Timer timer;
  timer.start();

  typedef CGAL::Periodic_3_Delaunay_triangulation_3<PDTT1>   P3TD3_K1;
  _test_periodic_3_triangulation_3_constructors( P3TD3_K1() );
  _test_cls_periodic_3_triangulation_3(P3TD3_K1(),
                                       PDTT1::Point(0.711476,-0.0713565,-0.52312),
                                       "data/P3DT3_covering_test_HOM.tri",
                                       "data/P3DT3_covering_test.tri");

  typedef CGAL::Periodic_3_Delaunay_triangulation_3<PDTT2>   P3TD3_K2;
  _test_periodic_3_triangulation_3_constructors( P3TD3_K2() );
  _test_cls_periodic_3_triangulation_3(P3TD3_K2(),
                                       PDTT2::Point(0.711476,-0.0713565,-0.52312),
                                       "data/P3DT3_covering_test_HOM.tri",
                                       "data/P3DT3_covering_test.tri",
                                       true /*exact*/);

  typedef CGAL::Periodic_3_Delaunay_triangulation_3<PDTT3>   P3TD3_K3;
  _test_periodic_3_triangulation_3_constructors( P3TD3_K3() );

  // commented because it takes too long to test it
//  _test_cls_periodic_3_triangulation_3(P3TD3_K3(),
//                                       PDTT3::Point(0.711476,-0.0713565,-0.52312),
//                                       "data/P3DT3_covering_test_HOM.tri",
//                                       "data/P3DT3_covering_test.tri",
//                                       false /*not exact*/,
//                                       true /*homogeneous*/,
//                                       false /*do not test input/output*/);

  std::cerr << timer.time() << " sec." << std::endl;
  return 0;
}
