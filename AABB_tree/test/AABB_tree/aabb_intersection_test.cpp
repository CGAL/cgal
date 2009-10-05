// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#include <string>

#include <CGAL/AABB_intersections.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>


template <class Triangle, class Query, class Result>
bool test_aux(const Triangle t,
              const Query& q,
              const std::string& name,
              const Result& expected)
{
  CGAL::Object object = CGAL::intersection(t,q);
  const Result* pr = CGAL::object_cast<Result>(&object);
  
  if ( (NULL != pr) && (expected == *pr) )
  {
    return true;
  }
  else
  {
    std::cout << "ERROR: intersection(" << name
              << ") did not answer the expected result !";  
    
    if ( NULL != pr )
      std::cout << " (answer: ["<< *pr << "])";
    
    std::cout << std::endl;
  }
  
  return false;
}

template <class K>
bool test()
{
	// types
  typedef typename K::FT FT;
  typedef typename K::Line_3 Line;
  typedef typename K::Point_3 Point;
  typedef typename K::Segment_3 Segment;
  typedef typename K::Ray_3 Ray;
  typedef typename K::Line_3 Line;
  typedef typename K::Triangle_3 Triangle;

  /* -------------------------------------
  // Test data is something like that (in t supporting plane)
  // Triangle is (p1,p2,p3)
  //
  //       +E          +1
  //                 /   \
  //        +C     6+  +8  +4      +B
  //              /   +7    \
  //            3+-------+5--+2
  //     
  //         +F        +A      
  ------------------------------------- */
  
  Point p1(FT(1.), FT(0.), FT(0.));
  Point p2(FT(0.), FT(1.), FT(0.));
  Point p3(FT(0.), FT(0.), FT(1.));
  
  Triangle t(p1,p2,p3);
  
  // Edges of t 
  Segment s12(p1,p2);
  Segment s13(p1,p3);
  Segment s23(p2,p3);
  
  bool b = test_aux(t,s12,"t-s12",s12);
  b &= test_aux(t,s13,"t-s13",s13);
  b &= test_aux(t,s23,"t-s23",s23);

  // Inside points
  Point p4(FT(0.5), FT(0.5), FT(0.));
  Point p5(FT(0.), FT(0.75), FT(0.25));
  Point p6(FT(0.5), FT(0.), FT(0.5));
  Point p7(FT(0.25), FT(0.625), FT(0.125));
  Point p8(FT(0.5), FT(0.25), FT(0.25));
  
  Segment s14(p1,p4);
  Segment s41(p4,p1);
  Segment s24(p2,p4);
  Segment s42(p4,p2);
  Segment s15(p1,p5);
  Segment s25(p2,p5);
  Segment s34(p3,p4);
  Segment s35(p3,p5);
  Segment s36(p3,p6);
  Segment s45(p4,p5);
  Segment s16(p1,p6);
  Segment s26(p2,p6);
  Segment s62(p6,p2);
  Segment s46(p4,p6);
  Segment s48(p4,p8);
  Segment s56(p5,p6);
  Segment s17(p1,p7);
  Segment s67(p6,p7);
  Segment s68(p6,p8);
  Segment s86(p8,p6);
  Segment s78(p7,p8);
  Segment s87(p8,p7);
  
  b &= test_aux(t,s14,"t-s14",s14);
  b &= test_aux(t,s41,"t-s41",s41);
  b &= test_aux(t,s24,"t-s24",s24);
  b &= test_aux(t,s42,"t-s42",s42);
  b &= test_aux(t,s15,"t-s15",s15);
  b &= test_aux(t,s25,"t-s25",s25);
  b &= test_aux(t,s34,"t-s34",s34);
  b &= test_aux(t,s35,"t-s35",s35);
  b &= test_aux(t,s36,"t-s36",s36);
  b &= test_aux(t,s45,"t-s45",s45);
  b &= test_aux(t,s16,"t-s16",s16);
  b &= test_aux(t,s26,"t-s26",s26);
  b &= test_aux(t,s62,"t-s62",s62);
  b &= test_aux(t,s46,"t-s46",s46);
  b &= test_aux(t,s48,"t-s48",s48);
  b &= test_aux(t,s56,"t-s56",s56);
  b &= test_aux(t,s17,"t-t17",s17);
  b &= test_aux(t,s67,"t-t67",s67);
  b &= test_aux(t,s68,"t-s68",s68);
  b &= test_aux(t,s86,"t-s86",s86);
  b &= test_aux(t,s78,"t-t78",s78);
  b &= test_aux(t,s87,"t-t87",s87);
  
  // Outside points (in triangle plane)
  Point pA(FT(-0.5), FT(1.), FT(0.5));
  Point pB(FT(0.5), FT(1.), FT(-0.5));
  Point pC(FT(0.5), FT(-0.5), FT(1.));
  Point pE(FT(1.), FT(-1.), FT(1.));
  Point pF(FT(-1.), FT(0.), FT(2.));
  
  Segment sAB(pA,pB);
  Segment sBC(pB,pC);
  Segment s2E(p2,pE);
  Segment sE2(pE,p2);
  Segment s2A(p2,pA);
  Segment s6E(p6,pE);
  Segment sB8(pB,p8);
  Segment sC8(pC,p8);
  Segment s8C(p8,pC);
  Segment s1F(p1,pF);
  Segment sF6(pF,p6);
  
  b &= test_aux(t,sAB,"t-sAB",p2);
  b &= test_aux(t,sBC,"t-sBC",s46);
  b &= test_aux(t,s2E,"t-s2E",s26);
  b &= test_aux(t,sE2,"t-sE2",s62);
  b &= test_aux(t,s2A,"t-s2A",p2);
  b &= test_aux(t,s6E,"t-s6E",p6);
  b &= test_aux(t,sB8,"t-sB8",s48);
  b &= test_aux(t,sC8,"t-sC8",s68);
  b &= test_aux(t,s8C,"t-s8C",s86);
  b &= test_aux(t,s1F,"t-s1F",s13);
  b &= test_aux(t,sF6,"t-sF6",s36);
  
  // Outside triangle plane
  Point pa(FT(0.), FT(0.), FT(0.));
  Point pb(FT(2.), FT(0.), FT(0.));
  Point pc(FT(1.), FT(0.), FT(1.));
  Point pe(FT(1.), FT(0.5), FT(0.5));
  
  Segment sab(pa,pb);
  Segment sac(pa,pc);
  Segment sae(pa,pe);
  Segment sa8(pa,p8);
  Segment sb2(pb,p2);
  
  b &= test_aux(t,sab,"t-sab",p1);
  b &= test_aux(t,sac,"t-sac",p6);
  b &= test_aux(t,sae,"t-sae",p8);
  b &= test_aux(t,sa8,"t-sa8",p8);
  b &= test_aux(t,sb2,"t-sb2",p2);
  
	return b;
}

int main()
{
  std::cout << "Testing with Simple_cartesian<float>..." << std::endl ;
  bool b = test<CGAL::Simple_cartesian<float> >();
  
  std::cout << "Testing with Simple_cartesian<double>..." << std::endl ;
	b &= test<CGAL::Simple_cartesian<double> >();
  
  std::cout << "Testing with Cartesian<float>..." << std::endl ;
	b &= test<CGAL::Cartesian<float> >();
  
  std::cout << "Testing with Cartesian<double>..." << std::endl ;
	b &= test<CGAL::Cartesian<double> >();
  
  std::cout << "Testing with Exact_predicates_inexact_constructions_kernel..." << std::endl ;
  b &= test<CGAL::Exact_predicates_inexact_constructions_kernel>();
	
  std::cout << "Testing with Exact_predicates_exact_constructions_kernel..." << std::endl ;
  b &= test<CGAL::Exact_predicates_exact_constructions_kernel>();
  
  if ( b )
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}
