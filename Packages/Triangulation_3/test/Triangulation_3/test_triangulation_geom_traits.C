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
// file          : test_triangulation_geom_traits_2.C
// revision      : 
// revision_date : 
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/basic.h>

#include <CGAL/_test_types.h>

#include <CGAL/Triangulation_geom_traits_3.h>

#include <CGAL/_test_types.C>
#include <CGAL/_test_cls_geom_traits.C>

int
main()
{
  std::cout << "Testing Triangulation_geom_traits" << std::endl;
  std::cout << "   with Cartesian" << std::endl;

  typedef CGAL::Triangulation_geom_traits_3<Test_rep_cartesian> Cls1;

typedef Cls1::Point  Pt1;
  Pt1 p1[55] = {
	Pt1(1,3,0), Pt1(3,5,2), Pt1(3,8,0), 
	Pt1(1,3,1), Pt1(3,5,0), Pt1(7,5,0), 
	Pt1(4,4,0), Pt1(4,4,0), Pt1(4,4,1),
        Pt1(2,6,0), Pt1(1,1,0), Pt1(2,2,0), Pt1(1,3,1), 
        Pt1(0,0,0), Pt1(0,0,1), Pt1(1,0,1), Pt1(1,1,1),

        Pt1(1,2,0), Pt1(1,3,3), Pt1(1,1,0), Pt1(3,2,2),
        Pt1(2,6,0), Pt1(1,1,0), Pt1(2,2,0), Pt1(1,3,-1), 
        Pt1(0,0,0), Pt1(0,0,1), Pt1(1,0,1), Pt1(1,-1,1),

        Pt1(1,2,0), Pt1(1,3,3), Pt1(1,1,0), Pt1(-3,2,2),
       	Pt1(1,1,1), Pt1(2,3,1), Pt1(3,4,1), Pt1(5,6,1),
        Pt1(1,2,3), Pt1(0,2,5), Pt1(-5,2,4), Pt1(2,2,2),
        Pt1(3,0,0), Pt1(3,1,2), Pt1(3,5,2), Pt1(3,2,7),
        Pt1(0,0,0), Pt1(2,0,2), Pt1(0,2,2), Pt1(1,1,2),
        Pt1(1,3,0), Pt1(0,0,0), Pt1(1,0,0), Pt1(2,0,0),
        Pt1(1,-3,0), Pt1(3,0,0)
  };

   _test_cls_geom_traits(p1, Cls1());


  std::cout << "   with Homogeneous" << std::endl;
  typedef CGAL::Triangulation_geom_traits_3<Test_rep_homogeneous> Cls2;
  typedef Cls2::Point  Pt2;

  Pt2 p2[55] = {
	Pt2(1,3,0,1), Pt2(3,5,2,1), Pt2(6,16,0,2), 
	Pt2(1,3,1,1), Pt2(6,10,0,2), Pt2(7,5,0,1), 
	Pt2(4,4,0,1), Pt2(4,4,0,1), Pt2(4,4,1,1),
        Pt2(2,6,0,1), Pt2(1,1,0,1), Pt2(2,2,0,1), Pt2(1,3,1,1), 
        Pt2(0,0,0,10), Pt2(0,0,1,10), Pt2(1,0,1,10), Pt2(1,1,1,10),

        Pt2(1,2,0,125), Pt2(1,3,3,125), Pt2(1,1,0,125), Pt2(3,2,2,125),
        Pt2(2,6,0,1), Pt2(10,10,0,10), Pt2(2,2,0,1), Pt2(1,3,-1,1), 
        Pt2(0,0,0,1), Pt2(0,0,1,1), Pt2(1,0,1,1), Pt2(1,-1,1,1),

        Pt2(1,2,0,125), Pt2(1,3,3,125), Pt2(1,1,0,125), Pt2(-3,2,2,125),
       	Pt2(10,10,10,10), Pt2(4,6,2,2), Pt2(6,8,2,2), Pt2(5,6,1,1),
        Pt2(1,2,3,1), Pt2(0,2,5,1), Pt2(-5,2,4,1), Pt2(2,2,2,1),
        Pt2(3,0,0,1), Pt2(3,1,2,1), Pt2(3,5,2,1), Pt2(3,2,7,1),
        Pt2(0,0,0,1), Pt2(2,0,2,1), Pt2(0,2,2,1), Pt2(1,1,2,1),
        Pt2(1,3,0,1), Pt2(0,0,0,1), Pt2(1,0,0,1), Pt2(2,0,0,1),
        Pt2(1,-3,0,1), Pt2(3,0,0,1)
  };

   _test_cls_geom_traits(p2, Cls2());


   std::cout << "   Testing Triangulation_test_traits for the requirements" << std::endl;
   typedef _Triangulation_test_traits_3 Cls3;

   typedef Cls3::Point  Pt3;

   Pt3 p3[55] = {
	Pt3(1,3,0), Pt3(3,5,2), Pt3(3,8,0), 
	Pt3(1,3,1), Pt3(3,5,0), Pt3(7,5,0), 
	Pt3(4,4,0), Pt3(4,4,0), Pt3(4,4,1),
        Pt3(2,6,0), Pt3(1,1,0), Pt3(2,2,0), Pt3(1,3,1), 
        Pt3(0,0,0), Pt3(0,0,1), Pt3(1,0,1), Pt3(1,1,1),

        Pt3(1,2,0), Pt3(1,3,3), Pt3(1,1,0), Pt3(3,2,2),
        Pt3(2,6,0), Pt3(1,1,0), Pt3(2,2,0), Pt3(1,3,-1), 
        Pt3(0,0,0), Pt3(0,0,1), Pt3(1,0,1), Pt3(1,-1,1),

        Pt3(1,2,0), Pt3(1,3,3), Pt3(1,1,0), Pt3(-3,2,2),
       	Pt3(1,1,1), Pt3(2,3,1), Pt3(3,4,1), Pt3(5,6,1),
        Pt3(1,2,3), Pt3(0,2,5), Pt3(-5,2,4), Pt3(2,2,2),
        Pt3(3,0,0), Pt3(3,1,2), Pt3(3,5,2), Pt3(3,2,7),
        Pt3(0,0,0), Pt3(2,0,2), Pt3(0,2,2), Pt3(1,1,2),
        Pt3(1,3,0), Pt3(0,0,0), Pt3(1,0,0), Pt3(2,0,0),
        Pt3(1,-3,0), Pt3(3,0,0)
  };

   _test_cls_geom_traits(p3, Cls3());
}

