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
// file          : demo/BooleanOperations/test_bops.C
// source        : demo/BooleanOperations/test_bops.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                        Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

#include <CGAL/test_bops.h>

using namespace CGAL;
using namespace std;

int test_iso_rectangles(void) {
  Iso_rectangle_2 A(Point_2(0,5), Point_2(6,2));
  Iso_rectangle_2 B(Point_2(4,3), Point_2(7,0));
  list<Object> result;
  cout << "iso_rectangle A: " << A << endl;
  cout << "iso_rectangle B: " << B << endl;
  cout << "INTERSECTION(A,B)" << " - ";
  intersection(A, B, back_inserter(result));
  test_result_output(result);
  result= list<Object>();
  cout << "UNION(A,B)" << " - ";
  Union(A, B, back_inserter(result));
  test_result_output(result);
  result= list<Object>();
  cout << "DIFFERENCE(A,B)" << " - ";
  difference(A, B, back_inserter(result));
  test_result_output(result);
  result= list<Object>();
  cout << "DIFFERENCE(B,A)" << " - ";
  difference(B, A, back_inserter(result));
  test_result_output(result);
  return 0;
}

int test_triangles(void) {
  Triangle_2 A(Point_2(0,3), Point_2(2,0), Point_2(4,1));
  Triangle_2 B(Point_2(-1,1), Point_2(1,0), Point_2(3,3));
  list<Object> result;
  cout << "triangle A: " << A << endl;
  cout << "triangle B: " << B << endl;
  cout << "INTERSECTION" << endl;
  //intersection(A, B, back_inserter(result));
  //test_result_output(result);
  cout << "UNION" << endl;
  union(A, B, back_inserter(result));
  test_result_output(result);
  result= list<Object>();
  cout << "DIFFERENCE" << endl;
  difference(A, B, back_inserter(result));
  test_result_output(result);
  return 0;
}

int test_intersection(void) {
  vector<Point_2> vA(6), vB(4);
  test_input( vA, vB);
  Polygon_2 A(vA.begin(), vA.end()), B(vB.begin(),vB.end());

  cout << "INTERSECTION" << endl;
  cout << "============" << endl;
  cout << "polygon A: " << A << endl;
  cout << "polygon B: " << B << endl;

  list<Object> result;
  /*
  insert_iterator< list<Object> > result_iterator(result, result.begin());
  intersection(A,B, result_iterator);
  */
  intersection(A,B, back_inserter(result));
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
  difference(A,B, back_inserter(result) );
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
  union(A,B, back_inserter(result) );
  test_result_output(result);

  return 0;
}


int main(int argc, char *argv[])
{
  int s= -1;
  if( argc>1 ) s= atoi(argv[1]);

  switch(s) {
    case 0: 
      cout << "usage: test_bops [number]" << endl;
      cout << "   - tests boolean operations on two simple polygons," << endl;
      cout << "   - where 'number' can be one of the following:" << endl;
      cout << "                  1 ... intersection" << endl;
      cout << "                  2 ... difference" << endl;
      cout << "                  3 ... union " << endl;
      cout << "                  4 ... iso_rectangles" << endl;
      cout << "                  5 ... triangles" << endl;
      cout << "            default ... all tests" << endl;
      break;
    case 1: test_intersection(); break;
    case 2: test_difference(); break;
    case 3: test_union(); break;
    case 4: test_iso_rectangles(); break;
    case 5: test_triangles(); break;
    default:
      test_iso_rectangles();
      test_triangles();
      test_intersection();
      test_difference(); 
      test_union(); 
  }
 
  return 0;
}
