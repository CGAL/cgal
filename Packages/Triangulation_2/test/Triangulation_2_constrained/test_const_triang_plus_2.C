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
// file          : test/Triangulation_2_constrained/test_const_triang_plus_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec (Mariette.Yvinec@sophia.inria.fr)
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <CGAL/_test_types.h>

#include <CGAL/intersections.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <CGAL/_test_cls_hierarchy_2.C>
#include <CGAL/_test_cls_const_triang_plus_2.C>


int main()
{

  std::cout << "Testing constraint hirarchy" << std::endl;
  _test_cls_hierarchy_2();

  std::cout << "Testing constrained_triangulation_plus_2 "<<
    std::endl;
  typedef CGAL::Constrained_Delaunay_triangulation_2<TestK>   CDt;
  typedef CGAL::Constrained_triangulation_plus_2<CDt>   CDtplus;
  _test_cls_const_triang_plus_2(CDtplus());
  return 0;
}
