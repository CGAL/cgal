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
// file          : test_delaunay_triangulation.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <utility>
#include <list>
#include <vector>

#include <CGAL/_test_types.h>

#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/_test_types.C>
#include <CGAL/_test_cls_delaunay_triangulation_2.C>

int main()
{
  cout << "Testing Delaunay Triangulation_2 "; 
  cout << " with Triangulation_test_traits : " << endl;
  cout << " this use double type coordinates " << endl;

  typedef CGAL::_Triangulation_test_traits                       Gt;
  typedef CGAL::Triangulation_vertex_base_2<Gt>                  Vb;
  typedef CGAL::Triangulation_face_base_2<Gt>                    Fb;
  typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
  typedef CGAL::Delaunay_triangulation_2<Gt,Tds>                 Cls;

  _test_cls_delaunay_triangulation_2( Cls() );

  return 0;
}
