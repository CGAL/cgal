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
// source        : $URL$
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//               : Mariette Yvinec (Mariette.Yvinec@sophia.inria.fr)
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================
#include <CGAL/basic.h>
#include <CGAL/_test_types.h>

#include <CGAL/intersections.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/_test_cls_const_Del_triangulation_2.h>

// Explicit instantiation of the whole class :
template class CGAL::Constrained_Delaunay_triangulation_2<TestK>;

int main()
{

  std::cout << "Testing constrained_Delaunay_triangulation "<< std::endl;
  std::cout << " with No_intersection_tag : " << std::endl;
  typedef CGAL::Constrained_Delaunay_triangulation_2<TestK>        CDt2;

  _test_cls_const_Del_triangulation(CDt2());
  return 0;
}
