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
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_data_structure_using_list_2.h>
#include <CGAL/Constrained_triangulation_2.h>

#include <CGAL/_test_cls_constrained_triangulation_2.C>

int main()
{
  std::cout << "Testing constrained_triangulation "<< std::endl;
  std::cout << " with Triangulation_euclidean_traits_2 : " << std::endl;
  std::cout << " and double coordinates " << std::endl;
  
  typedef double coord_type;
  typedef CGAL::Cartesian<coord_type>  Rep;
  typedef CGAL::Triangulation_euclidean_traits_2<Rep>            Gt;
//   typedef CGAL::Triangulation_vertex_base_2<Gt>                  Vb;
//   typedef CGAL::Constrained_triangulation_face_base_2<Gt>        Fb;
//   typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
//   typedef CGAL::Constrained_triangulation_2<Gt,Tds>              CCls;
  typedef CGAL::Constrained_triangulation_2<Gt>              CCls;
  _test_cls_constrained_triangulation(CCls());

  return 0;
}
