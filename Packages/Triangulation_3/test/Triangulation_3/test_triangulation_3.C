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
// file          : test_triangulation_3.C
// revision      : 
// revision_date : 
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/basic.h>
#include <cassert>

#include <list>
#include <vector>

#include <CGAL/_test_types.h>
#include <CGAL/triple.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Triangulation_data_structure_3.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/_test_cls_triangulation_3.C>

bool del = false;

int main()
{
  std::cout << " with Triangulation_test_traits_3 : " << std::endl;
  typedef CGAL::Triangulation_geom_traits_3<Test_rep_cartesian>     traits;
  // Using vertex_base_pointer induces a memory leak (not a bug, but the test
  // program is not adapted), so we use the normal vertex.
  // typedef CGAL::Triangulation_vertex_base_pointer_3<traits>        Vb;
  typedef CGAL::Triangulation_vertex_base_3<traits>                 Vb;
  typedef CGAL::Triangulation_cell_base_3<traits>                   Fb;
  typedef CGAL::Triangulation_data_structure_3<Vb,Fb>               Tds;
  typedef CGAL::Triangulation_3<traits,Tds>                         Cls3;

  _test_cls_triangulation_3( Cls3() );

  return 0;
}
