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
// file          : test/BooleanOperations/include/CGAL/test_bops_data.h
// source        : test/BooleanOperations/include/CGAL/test_bops_data.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                        Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

#ifndef CGAL_TEST_BOPS_DATA_H
#define CGAL_TEST_BOPS_DATA_H

using namespace std;

void test_input(vector<Point>& vA, vector<Point>& vB) {
  vA[0]= Point(2,4);
  vA[1]= Point(0,3);
  vA[2]= Point(1,1);
  vA[3]= Point(2,3);
  vA[4]= Point(3,1);
  vA[5]= Point(4,3);
  vB[0]= Point(0,2);
  vB[1]= Point(0,0);
  vB[2]= Point(5,0);
  vB[3]= Point(5,2);
}

void test_input_intersection_result( list<CGAL::Object>& result ) {
  //RESULT: (size=2)
  //pgon convex, n=3, CCW, 1 1, 1.5 2, 0.5 2
  //pgon convex, n=3, CCW, 2.5 2, 3 1, 3.5 2
  Container co1, co2;
  co1.push_back(Point(1,1,1));
  co1.push_back(Point(3,4,2));
  co1.push_back(Point(1,4,2));

  co2.push_back(Point(5,4,2));
  co2.push_back(Point(3,1,1));
  co2.push_back(Point(7,4,2));
  
  Polygon pgon1(co1.begin(), co1.end()), pgon2(co2.begin(), co2.end());
  result.push_back(CGAL::make_object(pgon1) );
  result.push_back(CGAL::make_object(pgon2) );

  return;
}

void test_input_difference_result( list<CGAL::Object>& result ) {
//RESULT: (size=1)
//pgon simple, n=8, CCW, 0 3, 0.5 2, 1.5 2, 2 3, 2.5 2, 3.5 2, 4 3, 2 4
  Container co1;
  co1.push_back(Point(0,3,1));
  co1.push_back(Point(1,4,2));
  co1.push_back(Point(3,4,2));
  co1.push_back(Point(2,3,1));
  co1.push_back(Point(5,4,2));
  co1.push_back(Point(7,4,2));
  co1.push_back(Point(4,3,1));
  co1.push_back(Point(2,4,1));

  Polygon pgon1(co1.begin(), co1.end());
  result.push_back(CGAL::make_object(pgon1) );
  return;

}

void test_input_union_result( list<CGAL::Object>& result ) {
//RESULT: (size=2)
//pgon simple, n=9, CCW, 0 3, 0.5 2, 0 2, 0 0, 5 0, 5 2, 3.5 2, 4 3, 2 4
//pgon convex, n=3, CLW, 1.5 2, 2 3, 2.5 2
  Container co1, co2;
  co1.push_back(Point(0,3,1));
  co1.push_back(Point(1,4,2));
  co1.push_back(Point(0,2,1));
  co1.push_back(Point(0,0,1));
  co1.push_back(Point(5,0,1));
  co1.push_back(Point(5,2,1));
  co1.push_back(Point(7,4,2));
  co1.push_back(Point(4,3,1));
  co1.push_back(Point(2,4,1));

  co2.push_back(Point(3,4,2));
  co2.push_back(Point(2,3,1));
  co2.push_back(Point(5,4,2));
  
  Polygon pgon1(co1.begin(), co1.end()), pgon2(co2.begin(), co2.end());
  result.push_back(CGAL::make_object(pgon1) );
  result.push_back(CGAL::make_object(pgon2) );

  return;
}

int test_result_compare( const list<CGAL::Object>& result ,
                          const list<CGAL::Object>& original ) {

  if( result.size() != original.size() ) return 1;
  list<CGAL::Object>::const_iterator it1, it2;

  Point pt1, pt2;
  Segment seg1, seg2;
  Polygon pgon1, pgon2;
  Iso_rectangle irect1, irect2;
  
  for( it1= result.begin(), it2= original.begin();
       it1 != result.end();
       it1++, it2++) {

    if( CGAL::assign( pgon1, *it1) && CGAL::assign( pgon2, *it2)) { /* polygon */
      if( pgon1.size() != pgon2.size() ) return 1;
      //if( pgon1 != pgon2 ) return 1;
    }
    else if( CGAL::assign( seg1, *it1) && CGAL::assign( seg2, *it2)) { /* segment */
      if( seg1 != seg2 ) return 1;
    }
    else if( CGAL::assign( pt1, *it1) && CGAL::assign( pt2, *it2))  { /* point */
      if( pt1 != pt2 ) return 1;
    }
    else if( CGAL::assign( irect1, *it1) && CGAL::assign( irect2, *it2))  { 
      if( irect1 != irect2 ) return 1;
    }
    else { /* nothing */
      return 1;
    }
  } 
  
  return 0;
}

#endif // CGAL_TEST_BOPS_DATA_H
