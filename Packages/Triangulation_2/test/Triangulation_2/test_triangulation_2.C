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

#include <CGAL/_test_types.h>

#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Triangulation_data_structure_using_list_2.h>
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

//   std::cout << std::endl << "Testing Triangulation_2 " <<std::endl; 
//   std::cout << " with Euclidean_traits_2<Homogeneous> : " << std::endl ;
//   typedef CGAL::Triangulation_euclidean_traits_2<Test_rep_homogeneous> Gt3;
//   typedef CGAL::Triangulation_vertex_base_2<Gt3>                     Vb3;
//   typedef CGAL::Triangulation_face_base_2<Gt3>                       Fb3;
//   typedef CGAL::Triangulation_default_data_structure_2<Gt3,Vb3,Fb3>  Tds3;
//   typedef CGAL::Triangulation_2<Gt3,Tds3>                            Cls3;

//   _test_cls_triangulation_2( Cls3() );

//   std::cout << std::endl << "Testing Triangulation_2" <<std::endl;
//   std::cout << " using Cartesaian Kernel traits : " << std::endl;
//   typedef CGAL::Cartesian<double>                                    Gt5;
//   typedef CGAL::Triangulation_vertex_base_2<Gt5>                     Vb5;
//   typedef CGAL::Triangulation_face_base_2<Gt5>                       Fb5;
//   typedef CGAL::Triangulation_default_data_structure_2<Gt5,Vb5,Fb5>  Tds5;
//   typedef CGAL::Triangulation_2<Gt5,Tds5>    Cls5;
//   _test_cls_triangulation_short_2( Cls5() );

  std::cout << std::endl << "Testing Triangulation_2" <<std::endl;
  std::cout << " using Homogeneous  Kernel traits : " << std::endl;
  typedef CGAL::Homogeneous<Rtype>                                     Gt6;
  typedef CGAL::Triangulation_vertex_base_2<Gt6>                     Vb6;
  typedef CGAL::Triangulation_face_base_2<Gt6>                       Fb6;
  typedef CGAL::Triangulation_default_data_structure_2<Gt6,Vb6,Fb6>  Tds6;
  typedef CGAL::Triangulation_2<Gt6,Tds6>    Cls6;
  _test_cls_triangulation_2( Cls6() );

  return 0;
}
