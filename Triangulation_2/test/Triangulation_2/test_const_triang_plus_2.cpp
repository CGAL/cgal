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
// source        : $URL$
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Mariette Yvinec (Mariette.Yvinec@sophia.inria.fr)
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/_test_types.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <CGAL/_test_cls_hierarchy_2.h>
#include <CGAL/_test_cls_const_triang_plus_2.h>


typedef CGAL::Constrained_Delaunay_triangulation_2<TestK> CDt_ni;
// Explicit instantiation of the whole class :
template class CGAL::Constrained_triangulation_plus_2<CDt_ni>;

int main()
{

  std::cout << "Testing constraint hierarchy" << std::endl;
  _test_cls_hierarchy_2();

  std::cout << "Testing constrained_triangulation_plus_2 "<< std::endl;
  std::cout << " with Exact_predicates_tag : " << std::endl;
  typedef CGAL::Triangulation_vertex_base_2<TestK>               Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<TestK>     Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>            TDS;
  typedef CGAL::Exact_predicates_tag                             Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<TestK,TDS,Itag> CDt;
  typedef CGAL::Constrained_triangulation_plus_2<CDt>   CDtplus;



  _test_cls_const_triang_plus_2(CDtplus());

  std::cout << "Testing constrained_triangulation_plus_2 "<<   std::endl;
  std::cout << " with Exact_intersections_tag : " << std::endl;
  typedef CGAL::Triangulation_vertex_base_2<EK>                 Vbb;
  typedef CGAL::Constrained_triangulation_face_base_2<EK>       Fbb;
  typedef CGAL::Triangulation_data_structure_2<Vbb,Fbb>         TDSS;
  typedef CGAL::Exact_intersections_tag                         EItag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<EK,TDSS,EItag>  CDtei;
  typedef CGAL::Constrained_triangulation_plus_2<CDtei>   CDteiplus;
  _test_cls_const_triang_plus_2(CDteiplus());

  return 0;
}
