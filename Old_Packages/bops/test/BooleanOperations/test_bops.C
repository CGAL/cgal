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
// file          : test/BooleanOperations/test_bops.C
// source        : test/BooleanOperations/test_bops.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :            Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

#include <CGAL/test_bops.h>
#include <iostream>


using namespace std;
using CGAL::Object;
using CGAL::intersection;
using CGAL::difference;
using CGAL::Union;
using std::cout;


int test_intersection(void) {
  vector<Point> vA(6), vB(4);
  test_input( vA, vB);
  Polygon A(vA.begin(), vA.end()), B(vB.begin(),vB.end());
  list<Object> result, correct_result;
  intersection(A,B, back_inserter(result));
  test_input_intersection_result( correct_result );
  return test_result_compare(result, correct_result);
}

int test_difference(void) {
  vector<Point> vA(6), vB(4);
  test_input( vA, vB);
  Polygon A(vA.begin(), vA.end()), B(vB.begin(),vB.end());
  list<Object> result, correct_result;
  difference(A,B, back_inserter(result) );
  test_input_difference_result( correct_result );
  return test_result_compare(result, correct_result);
}

int test_union(void) {
  vector<Point> vA(6), vB(4);
  test_input( vA, vB);
  Polygon A(vA.begin(), vA.end()), B(vB.begin(),vB.end());
  list<Object> result, correct_result;
  Union(A,B, back_inserter(result) );
  test_input_union_result( correct_result );
  return test_result_compare(result, correct_result);
}


int main(int argc, char *argv[])
{
  if( argc > 1 ) {
    cout << test_intersection() ? "intersection false" : "intersection OK";
    cout << endl;
    cout << test_difference() ? "difference false" : "difference OK";
    cout << endl;
    cout << test_union() ? "union false" : "union OK";
    cout << endl;
  }
  return test_intersection() || test_difference() || test_union(); 
}

