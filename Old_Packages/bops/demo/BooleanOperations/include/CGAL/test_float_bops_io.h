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
// file          : demo/BooleanOperations/include/CGAL/test_float_bops_io.h
// source        : demo/BooleanOperations/include/CGAL/test_float_bops_io.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                        Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

#ifndef TEST_FLOAT_BOPS_IO_H
#define TEST_FLOAT_BOPS_IO_H

ostream& operator<<(ostream& o, Orientation orient) {
  if(orient == CLOCKWISE )
      o << "CLW";
  else if(orient == COUNTERCLOCKWISE )
      o << "CCW";
  else /* COLLINEAR */
      o << "COL";
  return o;
}

ostream& operator<<(ostream& o, Polygon_2 p) {
  Polygon_2::Vertex_const_iterator it;
  if( p.is_convex() )      o << "convex,";
  /* else if( p.is_simple() ) o << "simple,";
  else                     o << "NOT simple,"; */
  else o << "simple,";
  o << " n=" << p.size() << ", " << p.orientation();
  for( it= p.vertices_begin(); it != p.vertices_end(); it++) 
    o  << ", " << *it;
  return o;
}

void test_input(vector<Point_2>& vA, vector<Point_2>& vB) {
/* 6
4.15 4.15
-0.2 6.55
-5.05 1.15
-4.65 -5.8
3.05 -5.15
5.85 0.3
*/
  vA.reserve(6);
  vA[0]= Point_2(4.15,4.15);
  vA[1]= Point_2(-0.2,6.55);
  vA[2]= Point_2(-5.05,1.15);
  vA[3]= Point_2(-4.65,-5.8);
  vA[4]= Point_2(3.05,-5.15);
  vA[5]= Point_2(5.85,0.3);

/* 6
8.1 -2.55
1.75 1.55
-3.5 4.6
-6.2 -2.15
1.8 -5.95
6.45 -4
*/
  vB.reserve(6);
  vB[0]= Point_2(8.1, -2.55);
  vB[1]= Point_2(1.75, 1.55);
  vB[2]= Point_2(-3.5, 4.6);
  vB[3]= Point_2(-6.2, -2.15);
  vB[4]= Point_2(1.8, -5.95);
  vB[5]= Point_2(6.45, -4);

}

void test_result_output( const list<Object>& result ) {
  Point_2 pt;
  Segment_2 seg;
  Polygon_2 pgon;
  Iso_rectangle_2 irect;

  list<Object>::const_iterator it;
  cout << endl << "RESULT: (size=" << result.size() << ")" << endl;

  for( it= result.begin(); it != result.end(); it++) {
    if( assign( pgon, *it) ) { /* polygon */
      //cout << "POLYGON" << endl;
      cout << "pgon " << pgon << endl;
    }
    else if( assign( seg, *it) ) { /* segment */
      cout << "seg  " << seg << endl;
    }
    else if( assign( pt, *it) )  { /* point */
      cout << "pt   " << pt << endl;
    }
    else if( assign( irect, *it) )  { /* point */
      cout << "irect   " << irect << endl;
    }
    else { /* nothing */
      cout << "undefined object " << endl;
    }
  }
  cout << endl << endl << endl;
  
  return;
}

bool read_point(Point_2 &pt)
{
    double x, y;
    cin >> x >> y;
    if (!cin.good())
	return false;
    pt = Point_2(TestNum(x), TestNum(y));
    return true;
}

bool read_pgn(Polygon_2 &pgn)
{
    int n, i;
    cin >> n;
    if (!cin.good())
	return false;
    if (n < 3) {
	cin.clear(ios::failbit);
	return false;
    }
    vector<Point_2> points(n);
    for (i=0; i<n; i++) {
	if (!read_point(points[i]))
	    return false;
    }
    pgn = Polygon_2(points.begin(), points.end());
    return true;
}

#endif // TEST_FLOAT_BOPS_IO_H
