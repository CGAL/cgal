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
// source        :
// file          : test_triangulation_2.C
// revision      : 
// revision_date : 
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================



#include <CGAL/_test_types.h>

//#include <CGAL/_test_types.C>



#include <CGAL/Weighted_point_2.h>

#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/_test_cls_regular_triangulation_2.C>

int main()
{
  cout << "Testing Regular_triangulation_2" <<endl;
  cout << " with CGAL::Regular_triangulation_euclidean_traits_2 : "<<endl;
  cout << "using  double type coordinates " << endl;
  typedef CGAL::Regular_triangulation_euclidean_traits_2
    <Test_rep_cartesian, double>           RGt2;
  typedef CGAL::Triangulation_vertex_base_2<RGt2>                     RVb2;
  typedef CGAL::Regular_triangulation_face_base_2<RGt2>               RFb2;
  typedef CGAL::Triangulation_default_data_structure_2<RGt2,RVb2,RFb2>  RTds2;
  typedef CGAL::Regular_triangulation_2<RGt2,RTds2>                    RCls2;

    _test_cls_reg_triangulation_2( RCls2() );

  return 0;
}

