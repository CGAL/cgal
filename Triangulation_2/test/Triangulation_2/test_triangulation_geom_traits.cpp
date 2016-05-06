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
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/_test_types.h>

#include <CGAL/internal/Projection_traits_3.h>
#include <CGAL/_test_traits.h>
#include <CGAL/_test_cls_geom_traits.h>

// Just define our own traits for the tests here to prevent including
// the deprecated headers. Otherwise they trigger spurious warnings in
// the test suite.

namespace CGAL {
template < class R >
class Triangulation_euclidean_traits_yz_3
  : public CGAL::internal::Projection_traits_3<R,0>
{};
template < class R >
class Triangulation_euclidean_traits_xz_3
  : public CGAL::internal::Projection_traits_3<R,1>
{};
template < class R >
class Triangulation_euclidean_traits_xy_3
  : public CGAL::internal::Projection_traits_3<R,2>
{};
} // CGAL

int
main()
{
  std::cout << "Testing with Cartesian" << std::endl;
  typedef Test_rep_cartesian Cls1;
  typedef Cls1::Point_2  Pt1;
  Pt1 p1[34] = {
	Pt1(1,3), Pt1(3,5), Pt1(3,8), 
	Pt1(1,3), Pt1(3,5), Pt1(7,5), 
	Pt1(4,4), Pt1(4,4), Pt1(2,6), 
	Pt1(1,1), Pt1(2,2), Pt1(1,3), 
	Pt1(1,1), Pt1(2,3), Pt1(5,0), 
	Pt1(6,0), Pt1(0,6), Pt1(3,3),
	Pt1(91,312), Pt1(125,300), Pt1(204,253), Pt1(324,0),
	Pt1(91,312), Pt1(125,300), Pt1(204,253), Pt1(326,0),
	Pt1(91,312), Pt1(125,300), Pt1(204,253), Pt1(325,0),
	Pt1(91,312), Pt1(125,300), Pt1(204,253), Pt1(0,0)
  };
  _test_cls_delaunay_geom_traits( p1, Cls1() );

  std::cout << "   with Kernel traits" << std::endl;
  typedef TestK Cls1b;
  typedef Cls1b::Point_2  Pt1b;
  Pt1b p1b[34] = {
	Pt1b(1,3), Pt1b(3,5), Pt1b(3,8), 
	Pt1b(1,3), Pt1b(3,5), Pt1b(7,5), 
	Pt1b(4,4), Pt1b(4,4), Pt1b(2,6), 
	Pt1b(1,1), Pt1b(2,2), Pt1b(1,3), 
	Pt1b(1,1), Pt1b(2,3), Pt1b(5,0), 
	Pt1b(6,0), Pt1b(0,6), Pt1b(3,3),
	Pt1b(91,312), Pt1b(125,300), Pt1b(204,253), Pt1b(324,0),
	Pt1b(91,312), Pt1b(125,300), Pt1b(204,253), Pt1b(326,0),
	Pt1b(91,312), Pt1b(125,300), Pt1b(204,253), Pt1b(325,0),
	Pt1b(91,312), Pt1b(125,300), Pt1b(204,253), Pt1b(0,0)
  };
  _test_cls_delaunay_geom_traits( p1b, Cls1b() );

  std::cout << "   with Homogeneous" << std::endl;
  typedef Test_rep_homogeneous Cls2;
  typedef Cls2::Point_2  Pt2;
  Pt2 p2[34] = {
	Pt2(1,3,1), Pt2(6,10,2), Pt2(3,8,1), 
	Pt2(1,3,1), Pt2(3,5,1), Pt2(14,10,2), 
	Pt2(4,4,1), Pt2(8,8,2), Pt2(2,6,1), 
	Pt2(1,1,1), Pt2(2,2,1), Pt2(1,3,1), 
	Pt2(1,1,1), Pt2(2,3,1), Pt2(5,0,1), 
	Pt2(6,0,1), Pt2(0,6,1), Pt2(6,6,2),
	Pt2(91,312,325), Pt2(125,300,325), Pt2(204,253,325), Pt2(324,0,325),
	Pt2(91,312,325), Pt2(125,300,325), Pt2(204,253,325), Pt2(326,0,325),
	Pt2(91,312,325), Pt2(125,300,325), Pt2(204,253,325), Pt2(1,0,1),
	Pt2(91,312,325), Pt2(125,300,325), Pt2(204,253,325), Pt2(0,0,43)
  };
  _test_cls_delaunay_geom_traits( p2, Cls2() );

  std::cout << "Testing Triangulation_euclidean_traits_xy_3" << std::endl;
  std::cout << "   with Cartesian" << std::endl;
  typedef CGAL::Triangulation_euclidean_traits_xy_3<Test_rep_cartesian> Cls3;
  typedef Cls3::Point_2  Pt3;
  Pt3 p3[34] = {
	Pt3(1,3,5), Pt3(3,5,4), Pt3(3,8,1), 
	Pt3(1,3,6), Pt3(3,5,3), Pt3(7,5,7), 
	Pt3(4,4,9), Pt3(4,4,4), Pt3(2,6,5), 
	Pt3(1,1,3), Pt3(2,2,7), Pt3(1,3,2), 
	Pt3(1,1,4), Pt3(2,3,2), Pt3(5,0,3), 
	Pt3(6,0,6), Pt3(0,6,1), Pt3(3,3,4),
	Pt3(91,312,5), Pt3(125,300,5), Pt3(204,253,5), Pt3(324,0,5),
	Pt3(91,312,5), Pt3(125,300,5), Pt3(204,253,5), Pt3(326,0,5),
	Pt3(91,312,5), Pt3(125,300,5), Pt3(204,253,5), Pt3(325,0,5),
	Pt3(91,312,5), Pt3(125,300,5), Pt3(204,253,5), Pt3(0,0,5)
  };
  _test_cls_geom_traits( p3, Cls3() );

  std::cout << "   with Homogeneous" << std::endl;
  typedef CGAL::Triangulation_euclidean_traits_xy_3<Test_rep_homogeneous> Cls4;
  typedef Cls4::Point_2  Pt4;
  Pt4 p4[34] = {
	Pt4(1,3,6,1), Pt4(6,10,6,2), Pt4(3,8,6,1), 
	Pt4(1,3,6,1), Pt4(3,5,6,1), Pt4(14,10,6,2), 
	Pt4(4,4,6,1), Pt4(8,8,6,2), Pt4(2,6,6,1), 
	Pt4(1,1,6,1), Pt4(2,2,6,1), Pt4(1,3,6,1), 
	Pt4(1,1,6,1), Pt4(2,3,6,1), Pt4(5,0,6,1), 
	Pt4(6,0,6,1), Pt4(0,6,6,1), Pt4(6,6,6,2),
	Pt4(91,312,234,325), Pt4(125,300,534,325), Pt4(204,253,241,325), Pt4(324,0,214,325),
	Pt4(91,312,21,325),  Pt4(125,300,534,325), Pt4(204,253,52,325), Pt4(326,0,421,325),
	Pt4(91,312,42,325),  Pt4(125,300,42,325), Pt4(204,253,542,325), Pt4(1,0,441,1),
	Pt4(91,312,534,325), Pt4(125,300,21,325), Pt4(204,253,123,325), Pt4(0,423,0,43)
  };
  _test_cls_geom_traits( p4, Cls4() );

  std::cout << "Testing Triangulation_euclidean_traits_yz_3" << std::endl;

  std::cout << "   with Cartesian" << std::endl;
  typedef CGAL::Triangulation_euclidean_traits_yz_3<Test_rep_cartesian> Cls5;
  typedef Cls5::Point_2  Pt5;
  Pt5 p5[34] = {
	Pt5(1,1,3), Pt5(3,3,5), Pt5(6,3,8), 
	Pt5(4,1,3), Pt5(8,3,5), Pt5(8,7,5), 
	Pt5(5,4,4), Pt5(1,4,4), Pt5(4,2,6), 
	Pt5(2,1,1), Pt5(6,2,2), Pt5(2,1,3), 
	Pt5(4,1,1), Pt5(3,2,3), Pt5(5,5,0), 
	Pt5(7,6,0), Pt5(2,0,6), Pt5(8,3,3),
	Pt5(1,91,312), Pt5(1,125,300), Pt5(1,204,253), Pt5(1,324,0),
	Pt5(1,91,312), Pt5(1,125,300), Pt5(1,204,253), Pt5(1,326,0),
	Pt5(1,91,312), Pt5(1,125,300), Pt5(1,204,253), Pt5(1,325,0),
	Pt5(1,91,312), Pt5(1,125,300), Pt5(1,204,253), Pt5(1,0,0)
  };
  _test_cls_geom_traits( p5, Cls5() );

  std::cout << "   with Homogeneous" << std::endl;
  typedef CGAL::Triangulation_euclidean_traits_yz_3<Test_rep_homogeneous> Cls6;
  typedef Cls6::Point_2  Pt6;
  Pt6 p6[34] = {
	Pt6(5,1,3,1), Pt6(5,6,10,2), Pt6(5,3,8,1), 
	Pt6(5,1,3,1), Pt6(5,3,5,1), Pt6(5,14,10,2), 
	Pt6(5,4,4,1), Pt6(5,12,12,3), Pt6(5,2,6,1), 
	Pt6(5,1,1,1), Pt6(5,2,2,1), Pt6(5,1,3,1), 
	Pt6(5,1,1,1), Pt6(5,2,3,1), Pt6(5,5,0,1), 
	Pt6(5,6,0,1), Pt6(5,0,6,1), Pt6(5,6,6,2),
	Pt6(325,91,312,325), Pt6(123,125,300,325), Pt6(12,204,253,325), Pt6(342,324,0,325),
	Pt6(765,91,312,325), Pt6(534,125,300,325), Pt6(542,204,253,325), Pt6(352,326,0,325),
	Pt6(876,91,312,325), Pt6(722,125,300,325), Pt6(513,204,253,325), Pt6(5,1,0,1),
	Pt6(533,91,312,325), Pt6(453,125,300,325), Pt6(32,204,253,325), Pt6(532,0,0,43)
  };
  _test_cls_geom_traits( p6, Cls6() );

  std::cout << "Testing Triangulation_euclidean_traits_xz_3" << std::endl;

  std::cout << "   with Cartesian" << std::endl;
  typedef CGAL::Triangulation_euclidean_traits_xz_3<Test_rep_cartesian> Cls7;
  typedef Cls7::Point_2  Pt7;
  Pt7 p7[34] = {
	Pt7(1,1,3), Pt7(3,3,5), Pt7(3,6,8), 
	Pt7(1,4,3), Pt7(3,8,5), Pt7(7,8,5), 
	Pt7(4,5,4), Pt7(4,1,4), Pt7(2,4,6), 
	Pt7(1,2,1), Pt7(2,6,2), Pt7(1,2,3), 
	Pt7(1,4,1), Pt7(2,3,3), Pt7(5,5,0), 
	Pt7(6,7,0), Pt7(0,2,6), Pt7(3,8,3),
	Pt7(91,3,312), Pt7(125,3,300), Pt7(204,3,253), Pt7(324,42,0),
	Pt7(91,3,312), Pt7(125,3,300), Pt7(204,3,253), Pt7(326,42,0),
	Pt7(91,3,312), Pt7(125,3,300), Pt7(204,3,253), Pt7(325,42,0),
	Pt7(91,3,312), Pt7(125,3,300), Pt7(204,3,253), Pt7(0,23,0)
  };

  _test_cls_geom_traits( p7, Cls7() );
  std::cout << "   with Homogeneous" << std::endl;
  typedef CGAL::Triangulation_euclidean_traits_xz_3<Test_rep_homogeneous> Cls8;
  typedef Cls8::Point_2  Pt8;
  Pt8 p8[34] = {
	Pt8(1,7,3,1), Pt8(6,7,10,2), Pt8(3,7,8,1), 
	Pt8(1,7,3,1), Pt8(3,7,5,1), Pt8(14,7,10,2), 
	Pt8(4,7,4,1), Pt8(12,7,12,3), Pt8(2,7,6,1), 
	Pt8(1,7,1,1), Pt8(2,7,2,1), Pt8(1,7,3,1), 
	Pt8(1,7,1,1), Pt8(2,7,3,1), Pt8(5,7,0,1), 
	Pt8(6,7,0,1), Pt8(0,7,6,1), Pt8(6,7,6,2),
	Pt8(91,653,312,325), Pt8(125,632,300,325), Pt8(204,632,253,325), Pt8(324,632,0,325),
	Pt8(91,632,312,325), Pt8(125,632,300,325), Pt8(204,632,253,325), Pt8(326,632,0,325),
	Pt8(91,632,312,325), Pt8(125,632,300,325), Pt8(204,632,253,325), Pt8(1,632,0,1),
	Pt8(91,632,312,325), Pt8(125,632,300,325), Pt8(204,632,253,325), Pt8(0,632,0,43)
  };
  _test_cls_geom_traits( p8, Cls8() );

  std::cout << "Testing Triangulation_test_traits for the requirements" 
	    << std::endl;
  typedef CGAL::_Triangulation_test_traits Cls9;
  typedef Cls9::Point_2  Pt9;
  Pt9 p9[34] = {
	Pt9(1,3), Pt9(3,5), Pt9(3,8), 
	Pt9(1,3), Pt9(3,5), Pt9(7,5), 
	Pt9(4,4), Pt9(4,4), Pt9(2,6), 
	Pt9(1,1), Pt9(2,2), Pt9(1,3), 
	Pt9(1,1), Pt9(2,3), Pt9(5,0), 
	Pt9(6,0), Pt9(0,6), Pt9(3,3),
	Pt9(91,312), Pt9(125,300), Pt9(204,253), Pt9(324,0),
	Pt9(91,312), Pt9(125,300), Pt9(204,253), Pt9(326,0),
	Pt9(91,312), Pt9(125,300), Pt9(204,253), Pt9(325,0),
	Pt9(91,312), Pt9(125,300), Pt9(204,253), Pt9(0,0)
  };
  _test_cls_delaunay_geom_traits( p9, Cls9() );
    
  return 0;
}
