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
// file          : Triangulation/test/Triangulation/test_triangulation_tdsul.C
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : Mariette Yvinec (Mariette.Yvinec@@sophia.inria.fr)
//
// ============================================================================

#include <CGAL/basic.h>
#include <CGAL/_test_types.h>

#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_using_list_2.h>

#include <CGAL/_test_types.C>
#include <CGAL/_test_cls_tds_2.C>

int
main()
{
  
  std::cout << "Testing Triangulation_data_structure_using_list_2" 
	    << std::endl;
  typedef CGAL::_Triangulation_test_traits Gt;
  typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
  typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
  typedef CGAL::Triangulation_data_structure_using_list_2<Vb,Fb> Cls;

  _test_cls_tds_2( Cls(), Gt() );
  
  return 0;
}
