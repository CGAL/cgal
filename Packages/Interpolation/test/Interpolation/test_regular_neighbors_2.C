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
// file          : /test/Interpolation/test_regular_neighbors_2.C
// package       : Interpolation
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Julia Floetotto <Julia.Flototto@sophia.inria.fr>
//
// coordinator   : 
// ============================================================================
#include <CGAL/basic.h>
#include <CGAL/double.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Regular_neighbor_coordinates_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <CGAL/_test_regular_neighbors_2.C>


struct K : CGAL::Exact_predicates_exact_constructions_kernel {};
typedef double W;

typedef CGAL::Regular_neighbor_coordinates_traits_2<K,W>     Gt1;
typedef CGAL::Regular_triangulation_2<Gt1>                   Rt1;

int main()
{

  std::cout << "Testing NN_neighbors_2 " << std::endl; 
  std::cout << " with Exact_predicates_exact_constructions_kernel: " << std::endl ;
  _test_regular_neighbors_2( Rt1() );
  
  return 0;
};
