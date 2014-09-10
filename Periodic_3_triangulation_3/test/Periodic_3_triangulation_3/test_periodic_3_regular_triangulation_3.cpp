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

#include <CGAL/Periodic_3_Regular_triangulation_traits_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include <CGAL/Gmpz.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Lazy_exact_nt.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Periodic_3_Regular_triangulation_3.h>

#include <cassert>


typedef CGAL::Epick K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Regular_traits;
typedef CGAL::Periodic_3_Regular_triangulation_traits_3<Regular_traits> Traits;
typedef typename Traits::Weighted_point Weighted_point;
typedef typename Traits::Bare_point Bare_point;
typedef typename Traits::Iso_cuboid_3 Iso_cuboid;

template class CGAL::Periodic_3_Regular_triangulation_3<Traits>;

typedef CGAL::Periodic_3_Regular_triangulation_3<Traits> P3RT3;

int main ()
{
  P3RT3 p3rt3;

  Weighted_point p(0,0,0);
  p3rt3.insert(p);

  std::cout << "EXIT SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
