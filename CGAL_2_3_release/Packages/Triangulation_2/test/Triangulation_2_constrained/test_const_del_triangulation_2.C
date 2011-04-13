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
//               : Mariette Yvinec (Mariette.Yvinec@sophia.inria.fr)
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <CGAL/_test_types.h>

#include <CGAL/Triangulation_data_structure_using_list_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/_test_cls_const_Del_triangulation_2.C>

int main()
{


  std::cout << "Testing constrained_triangulation "<< std::endl;
  std::cout << " with Kernel traits Cartesian<double> : " << std::endl;
  typedef CGAL::Cartesian<double>                                    Gt2;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Gt2>        CDt2;

  _test_cls_const_Del_triangulation(CDt2());
 return 0;
}
