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
// file          : /test/Triangulation/test_triangulation_2.C
// package       : Triangulation
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann, Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
// ============================================================================
#include <CGAL/basic.h>
#include <utility>
#include <list>
#include <vector>

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
  std::cout << "Testing Triangulation_2 " << std::endl; 
  std::cout << " with Euclidean_traits_2<Cartesian> : " << std::endl ;
  typedef CGAL::Triangulation_euclidean_traits_2<Test_rep_cartesian> Gt1;
  typedef CGAL::Triangulation_vertex_base_2<Gt1>                     Vb1;
  typedef CGAL::Triangulation_face_base_2<Gt1>                       Fb1;
  typedef CGAL::Triangulation_default_data_structure_2<Gt1,Vb1,Fb1>  Tds1;
  typedef CGAL::Triangulation_2<Gt1,Tds1>                            Cls1;

  _test_cls_triangulation_2( Cls1() );

  std::cout << std::endl << "Testing Triangulation_2 " <<std::endl; 
  std::cout << " with Euclidean_traits_2<Homogeneous> : " << std::endl ;
  typedef CGAL::Triangulation_euclidean_traits_2<Test_rep_homogeneous> Gt3;
  typedef CGAL::Triangulation_vertex_base_2<Gt3>                     Vb3;
  typedef CGAL::Triangulation_face_base_2<Gt3>                       Fb3;
  typedef CGAL::Triangulation_default_data_structure_2<Gt3,Vb3,Fb3>  Tds3;
  typedef CGAL::Triangulation_2<Gt3,Tds3>                            Cls3;

  _test_cls_triangulation_2( Cls3() );

  std::cout << std::endl << "Testing Triangulation_2" <<std::endl;
  std::cout << " with Triangulation_test_traits : " << std::endl;
  std::cout << " this use double type coordinates " << std::endl;
  typedef CGAL::_Triangulation_test_traits                           Gt2;
  typedef CGAL::Triangulation_vertex_base_2<Gt2>                     Vb2;
  typedef CGAL::Triangulation_face_base_2<Gt2>                       Fb2;
  typedef CGAL::Triangulation_default_data_structure_2<Gt2,Vb2,Fb2>  Tds2;
  typedef CGAL::Triangulation_2<Gt2,Tds2>    Cls2;

  _test_cls_triangulation_2( Cls2() );
 
  return 0;
}
