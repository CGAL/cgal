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
// source        : $URL$
// revision      : $Id$
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
#include <CGAL/Triangulation_ds_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_2.h>

#include <CGAL/_test_traits.h>
#include <CGAL/_test_cls_triangulation_2.h>

// Explicit instantiation of the whole class :
template class CGAL::Triangulation_2<TestK>;

int main()
{

  std::cout << "Testing Triangulation_2 " << std::endl; 
  std::cout << " with Cartesian : " << std::endl ;
  typedef Test_rep_cartesian Gt1;
  typedef CGAL::Triangulation_vertex_base_2<Gt1>                     Vb1;
  typedef CGAL::Triangulation_face_base_2<Gt1>                       Fb1;
  typedef CGAL::Triangulation_data_structure_2<Vb1,Fb1> Tds1;
  typedef CGAL::Triangulation_2<Gt1,Tds1>                            Cls1;
   _test_cls_triangulation_2( Cls1() );


  std::cout << std::endl << "Testing Triangulation_2" << std::endl;
  std::cout << " using Homogeneous  Kernel traits : " << std::endl;
  std::cout << " and defaults setting " << std::endl;
  typedef CGAL::Homogeneous<Rtype>      Gt6;
  typedef CGAL::Triangulation_2<Gt6>    Cls6;
  _test_cls_triangulation_2( Cls6() );

  return 0;
}
