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
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : Mariette Yvinec (Mariette.Yvinec@@sophia.inria.fr)
//
// ============================================================================

#include <CGAL/basic.h>
#include <CGAL/_test_types.h>

#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>

#include <CGAL/_test_types.C>
#include <CGAL/_test_cls_tds_2.C>

int
main()
{
  std::cout << "Testing Triangulation_defaut_data_structure_2" 
	    << std::endl;
  typedef CGAL::_Triangulation_test_traits Gt;
  typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
  typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
  typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Cls;
  _test_cls_tds_2( Cls(), Gt() );  
  return 0;
}
