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
// file          : test_triangulation_2.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <pair.h>
#include <list.h>
#include <vector.h>

#include <CGAL/_test_types.h>

#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_2.h>

#include <CGAL/_test_types.C>
#include <CGAL/_test_cls_triangulation_2.C>

int main()
{
  cout << "Testing Triangulation_2 with Gmpz coordinates"; 
  // cout << " with Euclidean_traits_2<Cartesian> : " << endl;
  // typedef CGAL_Triangulation_euclidean_traits_2<Test_rep_cartesian> Gt1;
  // typedef CGAL_Triangulation_vertex_base_2<Gt1>                     Vb1;
  // typedef CGAL_Triangulation_face_base_2<Gt1>                       Fb1;
  // typedef CGAL_Triangulation_default_data_structure_2<Gt1,Vb1,Fb1>  Tds1;
  // typedef CGAL_Triangulation_2<Gt1,Tds1>                            Cls1;

  // CGAL__test_cls_triangulation_2( Cls1() );

  cout << " with Triangulation_test_traits : " << endl;
  typedef CGAL__Triangulation_test_traits                           Gt2;
  typedef CGAL_Triangulation_vertex_base_2<Gt2>                     Vb2;
  typedef CGAL_Triangulation_face_base_2<Gt2>                       Fb2;
  typedef CGAL_Triangulation_default_data_structure_2<Gt2,Vb2,Fb2>  Tds2;
  typedef CGAL_Triangulation_2<Gt2,Tds2>                            Cls2;

  CGAL__test_cls_triangulation_2( Cls2() );

  return 0;
}
