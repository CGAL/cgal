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
// source        : $URL$
// file          : test_delaunay_triangulation.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================
#include <CGAL/basic.h>


#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4661) // Explicit instantiation will not
                                //  instatiate template member functions
#endif


#include <iostream>

#include <CGAL/_test_types.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/_test_traits.h>
#include <CGAL/_test_cls_delaunay_triangulation_2.h>


// Explicit instantiation of the whole class :
template class CGAL::Delaunay_triangulation_2<TestK>;

int main()
{
  std::cout << "Testing Delaunay Triangulation_2 " << std::endl; 
  std::cout << " with Euclidean cartesian points : " << std::endl;
  typedef CGAL::Delaunay_triangulation_2<Test_rep_cartesian> Cls1;

  _test_cls_delaunay_triangulation_2( Cls1() );


  std::cout << "Testing Delaunay Triangulation_2 "<< std::endl; 
  std::cout << " with Triangulation_test_traits : " << std::endl;
  typedef CGAL::_Triangulation_test_traits                       Gt;
  typedef CGAL::Delaunay_triangulation_2<Gt>                 Cls;

  _test_cls_delaunay_triangulation_2( Cls() ); 

  std::cout << "Testing Delaunay Triangulation_2 " <<  std::endl;
  std::cout << " with Triangulation_data_structure_2 : " 
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
  typedef CGAL::Delaunay_triangulation_2<Gt4>                   Cls4;

  _test_cls_delaunay_triangulation_2( Cls4() );

  std::cout << "done" << std::endl;
  return 0;
}

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif
