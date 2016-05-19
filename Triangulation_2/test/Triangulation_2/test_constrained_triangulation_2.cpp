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
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/basic.h>
#include <CGAL/_test_types.h>

#include <CGAL/intersections.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/_test_cls_constrained_triangulation_2.h>

// Explicit instantiation of the whole class :
template class CGAL::Constrained_triangulation_2<TestK>;

int main()
{
  std::cout << "Testing constrained_triangulation "<< std::endl;
  std::cout << " with No_intersection_tag : " << std::endl;
  typedef CGAL::Constrained_triangulation_2<TestK>         Ct;
  _test_cls_constrained_triangulation(Ct());

  std::cout << "Testing constrained_triangulation "<< std::endl;
  std::cout << " with Exact_predicates_tag : " << std::endl;
  typedef CGAL::Triangulation_vertex_base_2<TestK>                 Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<TestK>       Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
  typedef CGAL::Exact_predicates_tag                               Itag;
  typedef CGAL::Constrained_triangulation_2<TestK,TDS,Itag>        Ctwi;
  _test_cls_constrained_triangulation(Ctwi());
  
  std::cout << "Testing constrained_triangulation "<< std::endl;
  std::cout << " with Exact_intersections_tag : " << std::endl;
  typedef CGAL::Triangulation_vertex_base_2<EK>                 Vbb;
  typedef CGAL::Constrained_triangulation_face_base_2<EK>         Fbb;
  typedef CGAL::Triangulation_data_structure_2<Vbb,Fbb>         TDSS;
  typedef CGAL::Exact_intersections_tag                         EItag;
  typedef CGAL::Constrained_triangulation_2<EK,TDSS,EItag>      Ctwei;
  _test_cls_constrained_triangulation(Ctwei());
  
  return 0;
}
