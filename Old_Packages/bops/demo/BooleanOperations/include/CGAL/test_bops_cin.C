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
// file          : demo/BooleanOperations/include/CGAL/test_bops_cin.C
// source        : demo/BooleanOperations/include/CGAL/test_bops_cin.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                        Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

int test_intersection(const Polygon_2& A, const Polygon_2& B) {
  list<Object> result;
  intersection(A,B, back_inserter(result));
  test_result_output(result);
  return 0;
}


int test_difference(const Polygon_2& A, const Polygon_2& B) {
  list<Object> result;
  difference(A,B, back_inserter(result));
  test_result_output(result);
  return 0;
}


int test_union(const Polygon_2& A, const Polygon_2& B) {
  list<Object> result;
  union(A,B, back_inserter(result));
  test_result_output(result);
  return 0;
}


int main( int argc, char *argv[] )
{
  int arg= 0;
  if( argc > 1 ) {
    arg= atoi(argv[1]);
  }

  Polygon_2 A,B;

  if (!read_pgn(A)) {
	cerr << "Polygon A is wrong!" << endl;
	return 1;
  }
  // check counterclockwise orientation
  if ( A.orientation() == CLOCKWISE) {
    A.reverse_orientation();
    cout << "orientation reversed" << endl;
  }
  cout << "polygon A: " << A << endl;

  if (!read_pgn(B)) {
	cerr << "Polygon B is wrong!" << endl;
	return 1;
  }
  // check counterclockwise orientation
  if ( B.orientation() == CLOCKWISE) {
    B.reverse_orientation();
    cout << "orientation reversed" << endl;
  }
  cout << "polygon B: " << B << endl;

  switch( arg ) {
  case 1:
    cout << "INTERSECTION: " << endl;
    test_intersection(A,B);
    break;
  case 2:
    cout << "UNION: A + B " << endl;
    test_union(A, B);
    break;
  case 3:
    cout << "DIFFERENCE:  A - B" << endl;
    test_difference(A, B);
    break;
  case 4:
    cout << "DIFFERENCE:  B - A" << endl;
    test_difference(B, A);
    break;
  default:

    cout << "INTERSECTION: " << endl;
    test_intersection(A,B);

    cout << "UNION: A + B " << endl;
    test_union(A, B);

    cout << "DIFFERENCE:  A - B" << endl;
    test_difference(A, B);

    cout << "DIFFERENCE:  B - A" << endl;
    test_difference(B, A);
  }
  return 0;
}

