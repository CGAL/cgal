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

#include <CGAL/basic.h>
#include <CGAL/_test_types.h>

#include <CGAL/intersections.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/_test_cls_constrained_triangulation_2.C>

typedef CGAL::Quotient<Ftype>                       Exact_type;
typedef CGAL::Simple_cartesian<Exact_type>          Exact_kernel;
struct EK : public Exact_kernel {};

int main()
{
  std::cout << "Testing constrained_triangulation "<< std::endl;
  std::cout << " with No_intersection_tag : " << std::endl;
  typedef CGAL::Triangulation_euclidean_traits_2<TestK>            Gt;
  typedef CGAL::Triangulation_vertex_base_2<TestK>                 Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<TestK>       Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
  typedef CGAL::No_intersection_tag                                Itag;
  typedef CGAL::Constrained_triangulation_2<TestK,TDS,Itag>        Ct;
  _test_cls_constrained_triangulation(Ct());

  std::cout << "Testing constrained_triangulation "<< std::endl;
  std::cout << " with Exact_predicates_tag : " << std::endl;
  typedef CGAL::Constrained_triangulation_2<TestK>        Ctwi;
  _test_cls_constrained_triangulation(Ctwi());
  
  std::cout << "Testing constrained_triangulation "<< std::endl;
  std::cout << " with Exact_intersections_tag : " << std::endl;
  

  typedef CGAL::Triangulation_vertex_base_2<EK>                 Vbb;
  typedef CGAL::Constrained_triangulation_face_base_2<EK>       Fbb;
  typedef CGAL::Triangulation_data_structure_2<Vbb,Fbb>         TDSS;
  typedef CGAL::Exact_intersections_tag                         EItag;
  typedef CGAL::Constrained_triangulation_2<EK,TDSS,EItag>      Ctwei;
  _test_cls_constrained_triangulation(Ctwei());
  
  return 0;
}
