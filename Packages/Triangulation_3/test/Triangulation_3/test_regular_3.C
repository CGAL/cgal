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
// source        :
// file          : test_regular_3.C
// revision      : 
// revision_date : 
// author(s)     : Monique Teillaud (Monique.Teillaud@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/basic.h>
#include <iostream>
#include <cassert>
#include <list>
#include <vector>

#include <CGAL/triple.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer my_NT;
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz my_NT;
#else
#include <CGAL/double.h>
typedef double my_NT;
#endif
#endif

#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_HOMOGENEOUS_H
#include <CGAL/Homogeneous.h>
#endif // CGAL_HOMOGENEOUS_H

#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include <CGAL/Triangulation_data_structure_3.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Regular_triangulation_3.h>

#include <CGAL/_test_types.C>
#include <CGAL/_test_cls_regular_3.C>


typedef CGAL::Cartesian<my_NT> Test_rep_cartesian;
typedef CGAL::Homogeneous<my_NT> Test_rep_homogeneous;

bool del=true;

int main()
{

  std::cout << " with CGAL::Regular_triangulation_euclidean_traits_3: " << std::endl;
 
  typedef CGAL::Regular_triangulation_euclidean_traits_3<Test_rep_cartesian>  traits;
// works with both geom_traits
 // typedef CGAL::_Triangulation_test_traits_3                         traits;
  typedef CGAL::Triangulation_vertex_base_3<traits>                 Vb;
  typedef CGAL::Triangulation_cell_base_3<traits>                   Fb;
  typedef CGAL::Triangulation_data_structure_3<Vb,Fb>               Tds;
  typedef CGAL::Regular_triangulation_3<traits,Tds>                Cls;

  _test_cls_regular_3( Cls() );

  return 0;
}
