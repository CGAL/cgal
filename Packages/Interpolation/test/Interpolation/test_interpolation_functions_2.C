// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : /test/Interpolation/test_interpolation_functions_2.C
// package       : Interpolation
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Julia Floetotto <Julia.Flototto@sophia.inria.fr>
//
// coordinator   : 
// ============================================================================
#include <CGAL/basic.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/_test_interpolation_functions_2.C>

//////////////////////////////// 
struct K : CGAL::Exact_predicates_exact_constructions_kernel {};
typedef CGAL::Delaunay_triangulation_2<K>            Dt;

struct K2 : CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Delaunay_triangulation_2<K2>            Dt2;

//////////////////////////////// 
int main()
{
  
  std::cout << "Testing interpolation functions with 2D NN neighbors " 
	    << std::endl; 
  std::cout << " using Exact_predicates_exact_constructions_kernel: " 
	    << std::endl ;
  _test_interpolation_functions_2_delaunay( Dt(), K::FT(0));
  

  std::cout << "Testing interpolation functions with 2D NN neighbors " 
	    << std::endl; 
  std::cout << " using Exact_predicates_inexact_constructions_kernel: " 
	    << std::endl ;
  _test_interpolation_functions_2_delaunay( Dt2(), K2::FT(1e-10));

  return 0;
}
//end of file
