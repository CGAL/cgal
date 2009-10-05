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


template <class T>
bool test_aux(const T& t,
              const std::string& name,
              const CGAL::Bbox_3& bbox,
              bool expected)
{
  bool b = CGAL::do_intersect(t,bbox);
  
  if ( b != expected )
    std::cout << "ERROR: do_intersect(" << name
              << ") did not answer the expected result !" << std::endl;
  
  return (b == expected);
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

  CGAL::Bbox_3 bbox(1.0,1.0,1.0,10.0,50.0,100.0);
  
  Point p1(FT(0.), FT(0.), FT(0.));
  Point p2(FT(0.), FT(100.), FT(100.));
  Point p3(FT(100.), FT(0.), FT(100.));
  Point p4(FT(100.), FT(100.), FT(0.));
  Point p5(FT(0.), FT(0.), FT(100.));
  Point p6(FT(100.), FT(100.), FT(100.));
  Point p7(FT(-1.), FT(-1.), FT(-1.));
  
  Point pA(FT(5.), FT(0.), FT(0.));
  Point pB(FT(5.), FT(10.), FT(100.));
  Point pC(FT(5.), FT(10.), FT(-1.));
  Point pE(FT(5.), FT(10.), FT(0.));
  
  Segment s12(p1,p2);
  Segment s13(p1,p3);
  Segment s14(p1,p4);
  Segment s15(p1,p5);
  Segment s16(p1,p6);
  Segment s17(p1,p7);
  Segment s71(p7,p1);
  
  Segment sAB(pA,pB);
  Segment sBA(pB,pA);
  Segment sBC(pB,pC);
  Segment sCE(pC,pE);
  Segment sEC(pE,pC);
  
  bool b = test_aux(s12,"s12",bbox,false);
  b &= test_aux(s13,"s13",bbox,false);
  b &= test_aux(s14,"s14",bbox,false);
  b &= test_aux(s15,"s15",bbox,false);
  b &= test_aux(s16,"s16",bbox,true);
  b &= test_aux(s17,"s17",bbox,false);
  b &= test_aux(s71,"s71",bbox,false);
  
  b &= test_aux(sAB,"sAB",bbox,true);
  b &= test_aux(sBA,"sBA",bbox,true);
  b &= test_aux(sBC,"sBC",bbox,true);
  b &= test_aux(sCE,"sCE",bbox,false);
  b &= test_aux(sEC,"sEC",bbox,false);
  
  
  Ray r12(p1,p2);
  Ray r13(p1,p3);
  Ray r14(p1,p4);
  Ray r15(p1,p5);
  Ray r16(p1,p6);
  Ray r17(p1,p7);
  Ray r71(p7,p1);
  
  Ray rAB(pA,pB);
  Ray rBA(pB,pA);
  Ray rBC(pB,pC);
  Ray rCE(pC,pE);
  Ray rEC(pE,pC);
  
  b &= test_aux(r12,"r12",bbox,false);
  b &= test_aux(r13,"r13",bbox,false);
  b &= test_aux(r14,"r14",bbox,false);
  b &= test_aux(r15,"r15",bbox,false);
  b &= test_aux(r16,"r16",bbox,true);
  b &= test_aux(r17,"r17",bbox,false);
  b &= test_aux(r71,"r71",bbox,true);
  
  b &= test_aux(rAB,"rAB",bbox,true);
  b &= test_aux(rBA,"rBA",bbox,true);
  b &= test_aux(rBC,"rBC",bbox,true);
  b &= test_aux(rCE,"rCE",bbox,true);
  b &= test_aux(rEC,"rEC",bbox,false);
    
  Line l12(p1,p2);
  Line l13(p1,p3);
  Line l14(p1,p4);
  Line l15(p1,p5);
  Line l16(p1,p6);
  Line l17(p1,p7);
  Line l71(p7,p1);
  
  Line lAB(pA,pB);
  Line lBA(pB,pA);
  Line lBC(pB,pC);
  Line lCE(pC,pE);
  Line lEC(pE,pC);
  
  b &= test_aux(l12,"l12",bbox,false);
  b &= test_aux(l13,"l13",bbox,false);
  b &= test_aux(l14,"l14",bbox,false);
  b &= test_aux(l15,"l15",bbox,false);
  b &= test_aux(l16,"l16",bbox,true);
  b &= test_aux(l17,"l17",bbox,true);
  b &= test_aux(l71,"l71",bbox,true);
  
  b &= test_aux(lAB,"lAB",bbox,true);
  b &= test_aux(lBA,"lBA",bbox,true);
  b &= test_aux(lBC,"lBC",bbox,true);
  b &= test_aux(lCE,"lCE",bbox,true);
  b &= test_aux(lEC,"lEC",bbox,true);

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
	
  if ( b )
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}
