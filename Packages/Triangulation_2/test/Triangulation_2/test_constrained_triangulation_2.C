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
// file          : test/Triangulation/test_constrained_triangulation.C
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <CGAL/_test_types.h>

#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Constrained_triangulation_2.h>

#include <CGAL/_test_types.C>
#include <CGAL/_test_cls_constrained_triangulation_2.C>

int main()
{
  std::cout << "Testing constrained_triangulation "<< endl;
  std::cout << " with Triangulation_test_traits : " << endl;
  std::cout << " this uses double type coordinates " << endl;
  typedef CGAL::_Triangulation_test_traits                       Gt;
  typedef CGAL::Triangulation_vertex_base_2<Gt>                  Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<Gt>        CFb;
  typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,CFb> Tds;
  typedef CGAL::Constrained_triangulation_2<Gt,Tds>              CCls;

  _test_cls_constrained_triangulation(CCls());
  return 0;
}
