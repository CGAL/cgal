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
// source        : $URL$
// file          : test/Triangulation/test_triangulation_2.C
// revision      : $revision$
// revision_date : $Date$
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//                 Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/_test_types.h>

#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/_test_cls_regular_triangulation_2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>


typedef CGAL::Exact_predicates_exact_constructions_kernel  RGt;

// Explicit instantiation of the whole class :
template class CGAL::Regular_triangulation_2<RGt>;

int main()
{
  std::cout << "Testing Regular_triangulation_2" <<std::endl;
  std::cout << "using  Cartesian  points "   <<  std::endl;
  typedef CGAL::Regular_triangulation_2<RGt>                    RCls;
  _test_cls_regular_triangulation_2( RCls() );

  std::cout << "Testing Regular_triangulation_2" <<std::endl;
  std::cout << "using  Homogeneous points "<< std::endl;
  typedef Test_rep_homogeneous      RGt2;
  typedef CGAL::Regular_triangulation_2<RGt2>                    RCls2;
  _test_cls_regular_triangulation_2( RCls2() );

  std::cout << "done" << std::endl;
  return 0;
}

