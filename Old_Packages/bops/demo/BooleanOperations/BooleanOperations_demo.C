//  -*- Mode: c++ -*-
// ============================================================================
// 
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-1.0 $
// release_date  : $CGAL_Date: 1998/09/12 $
//
// file          : demo/BooleanOperations/BooleanOperations_demo.C
// source        : demo/BooleanOperations/BooleanOperations_demo.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                        Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

/*
  This is an example program for the usage of boolean opertations:
  

  The number type it uses is a non-exact one:    Quotient<int>
                               ---------
  The coordinates are usually cartesian and if one wants to use
  homogeneous coordinates the compiler option -DCGAL_BOPS_HOMOGENEOUS
  has to be added.




  

*/

#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/basic.h>

#include <list>
#include <CGAL/boolean_operations_2.h>
#include <iostream>

using CGAL::Object;
using std::list;
using std::vector;
using std::cout;
using std::endl;

//#include <CGAL/Rational.h>
//typedef Rational TestNum;  // only with Leda usage
//typedef double TestNum;         // non exact
typedef CGAL::Quotient<long int>       TestNum;  // exact, but very finite

#ifdef BOPS_HOMOGENEOUS
  typedef CGAL::Homogeneous<TestNum>   R_type;
#else
  typedef CGAL::Cartesian<TestNum>     R_type;
#endif

typedef CGAL::Point_2<R_type>          Point_2;
typedef CGAL::Segment_2<R_type>        Segment_2;
                                    // Polygon_2
typedef list< Point_2 >               Container;
typedef CGAL::Polygon_traits_2<R_type> Polygon_traits_2;
typedef CGAL::Polygon_2< Polygon_traits_2, Container >  Polygon_2;

#include "CGAL/example_io.h"  /* I/O to cin/cout */



int test_intersection(void) {
  vector<Point_2> vA(6), vB(4);
  test_input( vA, vB);
  Polygon_2 A(vA.begin(), vA.end()), B(vB.begin(),vB.end());

  cout << "INTERSECTION" << endl;
  cout << "============" << endl;
  cout << "polygon A: " << A << endl;
  cout << "polygon B: " << B << endl;

  list<Object> result;
  CGAL::intersection(A,B, back_inserter(result));
  test_result_output(result);
  
  return 0;
}



int test_difference(void) {
  vector<Point_2> vA(6), vB(4);
  test_input( vA, vB);
  Polygon_2 A(vA.begin(), vA.end()), B(vB.begin(),vB.end());

  cout << "DIFFERENCE" << endl;
  cout << "==========" << endl;
  cout << "polygon A: " << A << endl;
  cout << "polygon B: " << B << endl;

  list<Object> result;
  CGAL::difference(A,B, back_inserter(result) );
  test_result_output(result);

  return 0;
}



int test_union(void) {
  vector<Point_2> vA(6), vB(4);
  test_input( vA, vB);
  Polygon_2 A(vA.begin(), vA.end()), B(vB.begin(),vB.end());

  cout << "UNION" << endl;
  cout << "=====" << endl;
  cout << "polygon A: " << A << endl;
  cout << "polygon B: " << B << endl;

  list<Object> result;
  CGAL::Union(A,B, back_inserter(result) );
  test_result_output(result);

  return 0;
}


int main(int argc, char *argv[])
{
  int s= -1;
  if( argc>1 ) s= atoi(argv[1]);

  switch(s) {
    case 1: test_intersection(); break;
    case 2: test_difference(); break;
    case 3: test_union(); break;
    default:
      test_intersection();
      test_difference(); 
      test_union(); 
  }
 
  return 0;
}
