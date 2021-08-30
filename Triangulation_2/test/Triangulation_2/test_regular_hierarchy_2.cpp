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
// file          : test/Triangulation/test_regular_hierarchy_2.C
// revision      : $revision$
// revision_date : $Date$
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//                 Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//                 Andreas Fabri    <Andreas.Fabri@geometryfactory.com>
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>

#include <CGAL/_test_cls_regular_triangulation_2.h>
#include <CGAL/_test_cls_regular_hierarchy_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel RGt;

typedef CGAL::Regular_triangulation_vertex_base_2<RGt> Vbb;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb> Vb;
typedef CGAL::Regular_triangulation_face_base_2<RGt>  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>  Tds;
typedef CGAL::Regular_triangulation_2<RGt,Tds>  Rt;
typedef CGAL::Triangulation_hierarchy_2<Rt> Regular_hierarchy_cartesian;

// Explicit instantiation of the whole class  :
template class CGAL::Triangulation_hierarchy_2<Rt>;


int main()
{
  std::cout << "Testing Triangulation_hierarchy_2<Regular_triangulation_2>"
            <<std::endl;
  std::cout << "using  Cartesian  points "   <<  std::endl;
  //_test_cls_regular_hierarchy_2( Rt());
  _test_cls_regular_hierarchy_2( Regular_hierarchy_cartesian());
  std::cout << "done" << std::endl;
  return 0;
}
