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

#include <CGAL/_test_types.h>

#include <CGAL/Triangulation_ds_vertex_base_2.h>
#include <CGAL/Triangulation_ds_face_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_2.h>

#include <CGAL/_test_traits.h>
#include <CGAL/_test_cls_triangulation_short_2.h>


int main()
{

  std::cout << std::endl << "Testing Triangulation_2" <<std::endl;
  std::cout << " with Triangulation_test_traits : " << std::endl;
  std::cout << " this use double type coordinates " << std::endl;
  typedef CGAL::_Triangulation_test_traits                           Gt2;
  typedef CGAL::Triangulation_vertex_base_2<Gt2>                     Vb2;
  typedef CGAL::Triangulation_face_base_2<Gt2>                       Fb2;
  typedef CGAL::Triangulation_data_structure_2<Vb2,Fb2>  Tds2;
  typedef CGAL::Triangulation_2<Gt2,Tds2>    Cls2;
  _test_cls_triangulation_short_2( Cls2() );

  std::cout << std::endl << "Testing Backward Compatibility"
	    <<std::endl;
  typedef CGAL::Triangulation_data_structure_using_list_2<Vb2,Fb2>  Tds3;
  typedef CGAL::Triangulation_2<Gt2,Tds3>    Cls3;
  _test_cls_triangulation_short_2( Cls3() );

  typedef CGAL::Triangulation_default_data_structure_2<Gt2,Vb2,Fb2>  Tds4;
  typedef CGAL::Triangulation_2<Gt2,Tds4>    Cls4;
  _test_cls_triangulation_short_2( Cls4() );

  return 0;
}
