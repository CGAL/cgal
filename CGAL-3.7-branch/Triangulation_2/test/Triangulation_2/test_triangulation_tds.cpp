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

#include <CGAL/Triangulation_data_structure_2.h>

#include <CGAL/_test_traits.h>
#include <CGAL/_test_cls_tds_2.h>



typedef CGAL::Triangulation_ds_vertex_base_2<>     Vb;
typedef CGAL::Triangulation_ds_face_base_2<>       Fb;

// Explicit instantiation :
// does not work because of off_file_input
// template class CGAL::Triangulation_data_structure_2<Vb,Fb>;

int main()
{
  std::cout << "Testing Triangulation_data_structure_2" 
	    << std::endl << std::endl;
    typedef CGAL::Triangulation_data_structure_2<> Cls1;
  _test_cls_tds_2( Cls1());  

  std::cout << "Testing bakward compatibility" << std::endl;
  std::cout << "Testing Triangulation_defaut_data_structure_2" 
	    << std::endl;
  typedef CGAL::_Triangulation_test_traits Gt;
  typedef CGAL::Triangulation_ds_vertex_base_2<> Vb;
  typedef CGAL::Triangulation_ds_face_base_2<>  Fb;
  typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Cls2;
  _test_cls_tds_2( Cls2());  

  std::cout << "Testing Triangulation_data_structure_using_list_2" 
	    << std::endl;
   typedef CGAL::Triangulation_data_structure_using_list_2<Vb,Fb> Cls3;

  _test_cls_tds_2( Cls3());
  return 0;
}
