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
// file          : test_triangulation_tds.C
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
#include <CGAL/_test_cls_tds_3.C>

int main()
{
  std::cout << " with Triangulation_test_traits_3 : " << std::endl;
  typedef    CGAL::Triangulation_geom_traits_3<Test_rep_cartesian>       traits;
//  typedef   CGAL::_Triangulation_test_traits                 traits;
  typedef CGAL::Triangulation_vertex_base_3<traits>                 Vb;
  typedef CGAL::Triangulation_cell_base_3<traits>                   Fb;
  typedef CGAL::Triangulation_data_structure_3<Vb,Fb>               Tds;
  _test_cls_tds_3(Tds());
  return 0;
}
