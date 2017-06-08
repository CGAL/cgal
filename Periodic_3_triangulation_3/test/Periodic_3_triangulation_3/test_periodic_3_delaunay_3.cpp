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
//
//
// Author(s)     : Nico Kruithof
//                 Manuel Caroli

#include <iostream>
#include <fstream>

#include <CGAL/Timer.h>

#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K1;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K1> PTT1;
// Explicit instantiation of the whole class :
template class CGAL::Periodic_3_Delaunay_triangulation_3<PTT1>;

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel K2;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K2> PTT2;
// Explicit instantiation of the whole class :
template class CGAL::Periodic_3_Delaunay_triangulation_3<PTT2>;

#include <CGAL/MP_Float.h>
#include <CGAL/Simple_homogeneous.h>
typedef CGAL::Simple_homogeneous<CGAL::MP_Float> K3;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K3> PTT3;
// Explicit instantiation of the whole class :
template class CGAL::Periodic_3_Delaunay_triangulation_3<PTT3>;

#include <CGAL/_test_cls_periodic_3_delaunay_3.h>
int main(int, char**)
{
  CGAL::Timer timer;
  timer.start();
  typedef CGAL::Periodic_3_Delaunay_triangulation_3< PTT1 > P3T3_1;
  _test_cls_periodic_3_delaunay_3( P3T3_1() );

  typedef CGAL::Periodic_3_Delaunay_triangulation_3< PTT2 > P3T3_2;
  _test_cls_periodic_3_delaunay_3( P3T3_2() );

  // typedef CGAL::Periodic_3_Delaunay_triangulation_3< PTT3 > P3T3_3;
  // this takes too much time for the test suite.
  //_test_cls_periodic_3_delaunay_3( P3T3_3(), true );

  std::cerr << timer.time() << " sec." << std::endl;
  return 0;
}
