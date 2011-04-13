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
// source        : $RCSfile$
// file          : test_delaunay_triangulation.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================
#include <CGAL/basic.h>
#include <iostream>
//#include <vector>

#include <CGAL/_test_types.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Triangulation_data_structure_using_list_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/_test_types.C>
#include <CGAL/_test_cls_delaunay_triangulation_2.C>

int main()
{
  std::cout << "Testing Delaunay Triangulation_2 " << std::endl; 
  std::cout << " with Euclidean cartesian points : " << std::endl;
  typedef CGAL::Triangulation_euclidean_traits_2<Test_rep_cartesian> Gt1;
  typedef CGAL::Triangulation_vertex_base_2<Gt1>                  Vb1;
  typedef CGAL::Triangulation_face_base_2<Gt1>                    Fb1;
  typedef CGAL::Triangulation_default_data_structure_2<Gt1,Vb1,Fb1> Tds1;
  typedef CGAL::Delaunay_triangulation_2<Gt1,Tds1>                 Cls1;

  _test_cls_delaunay_triangulation_2( Cls1() );


//   std::cout << "Testing Delaunay Triangulation_2 "<< std::endl; 
//   std::cout << " with Triangulation_test_traits : " << std::endl;
//   typedef CGAL::_Triangulation_test_traits                       Gt;
//   typedef CGAL::Triangulation_vertex_base_2<Gt>                  Vb;
//   typedef CGAL::Triangulation_face_base_2<Gt>                    Fb;
//   typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
//   typedef CGAL::Delaunay_triangulation_2<Gt,Tds>                 Cls;

//   _test_cls_delaunay_triangulation_2( Cls() ); 

  std::cout << "Testing Delaunay Triangulation_2 " <<  std::endl;
  std::cout << " with Triangulation_data_structure_using_list_2 : " 
	    <<  std::endl << " and Cartesian<double>"
    	    << std::endl;
  std::cout << "this tests defaults setting" << std::endl;
  typedef CGAL::Cartesian<double>                       Gt3;
  typedef CGAL::Delaunay_triangulation_2<Gt3>           Cls3;

  _test_cls_delaunay_triangulation_2( Cls3() );

   std::cout << "Testing Delaunay Triangulation_2 " <<  std::endl;
  std::cout << " using Kernel_traits : " 	    << std::endl;
  std::cout << " and Homogeneous Points "       << std::endl;
  typedef CGAL::Homogeneous<double>                                 Gt4;
  typedef CGAL::Triangulation_vertex_base_2<Gt4>                    Vb4;
  typedef CGAL::Triangulation_face_base_2<Gt4>                      Fb4;
  typedef CGAL::Triangulation_default_data_structure_2<Gt4,Vb4,Fb4> Tds4;
  typedef CGAL::Delaunay_triangulation_2<Gt4,Tds4>                   Cls4;

  _test_cls_delaunay_triangulation_2( Cls4() );

  return 0;
}
