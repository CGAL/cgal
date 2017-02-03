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
// Author(s)     : Aymeric PELLE <aymeric.pelle@sophia.inria.fr>

#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include <CGAL/Gmpz.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Lazy_exact_nt.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/_test_cls_periodic_3_regular_triangulation_traits_3.h>

#include <cassert>


void test_periodic_3_regular_triangulation_traits_3()
{
  std::cout << "EPECK" << std::endl;
  _test_cls_periodic_3_regular_triangulation_traits_3<CGAL::Exact_predicates_exact_constructions_kernel>();
  std::cout << "EPICK" << std::endl;
  _test_cls_periodic_3_regular_triangulation_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel>();
}
