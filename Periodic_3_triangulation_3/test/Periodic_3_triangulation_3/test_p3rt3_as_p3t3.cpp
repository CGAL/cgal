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
// Author(s)     : Francois Rebufat
//                 Manuel Caroli

#include <iostream>
#include <cassert>

#include <CGAL/Periodic_3_regular_triangulation_3.h>
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_triangulation_3.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include <CGAL/_test_cls_periodic_3_triangulation_3.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Epeck;
typedef CGAL::Regular_triangulation_euclidean_traits_3<Epeck> RT_Exact;
typedef CGAL::Periodic_3_regular_triangulation_traits_3<RT_Exact> PRTT_Exact;
// Explicit instantiation of the whole class:
template class CGAL::Periodic_3_regular_triangulation_3<PRTT_Exact>;

typedef CGAL::Exact_predicates_exact_constructions_kernel Epick;
typedef CGAL::Regular_triangulation_euclidean_traits_3<Epick> RT_Inexact;
typedef CGAL::Periodic_3_regular_triangulation_traits_3<RT_Inexact> PRTT_Inexact;
// Explicit instantiation of the whole class:
template class CGAL::Periodic_3_regular_triangulation_3<PRTT_Inexact>;


int main()
{
  typedef CGAL::Periodic_3_regular_triangulation_3<PRTT_Exact>            P3RT3_Exact;
  _test_periodic_3_triangulation_3_constructors( P3RT3_Exact() );
  _test_cls_periodic_3_triangulation_3( P3RT3_Exact() );

  typedef CGAL::Periodic_3_regular_triangulation_3<PRTT_Inexact>            P3RT3_Inexact;
  _test_periodic_3_triangulation_3_constructors( P3RT3_Inexact() );
  _test_cls_periodic_3_triangulation_3( P3RT3_Inexact(), true );

  return 0;
}
