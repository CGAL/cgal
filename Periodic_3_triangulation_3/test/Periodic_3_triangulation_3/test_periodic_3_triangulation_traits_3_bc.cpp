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
// Author(s)     : Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

//#include <CGAL/Ep>

#ifdef CGAL_NO_DEPRECATED_CODE
int main() { return 0;}
#else
#  define CGAL_NO_DEPRECATION_WARNINGS 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/_test_periodic_3_static_filters_bc.h>
#include <CGAL/_test_cls_periodic_3_triangulation_traits_3_bc.h>

int main()
{
  std::cout<<"Statically filtered predicates:"<<std::endl;
  _test_periodic_3_static_filters();

  std::cout<<"\nTesting predicates:"<<std::endl;
  std::cout<<"  Predefined kernels...";std::cout.flush();
  _test_cls_periodic_3_triangulation_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel>();
  _test_cls_periodic_3_triangulation_traits_3<CGAL::Exact_predicates_exact_constructions_kernel>();
  std::cout<<" done"<<std::endl;

  return 0;
}

#endif

