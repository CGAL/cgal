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
// file          : demo/BooleanOperations/include/CGAL/example_io.h
// source        : demo/BooleanOperations/include/CGAL/example_io.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                        Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

#ifndef EXAMPLE_IO_H
#define EXAMPLE_IO_H

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
  vA[0]= Point_2(2,4);
  vA[1]= Point_2(0,3);
  vA[2]= Point_2(1,1);
  vA[3]= Point_2(2,3);
  vA[4]= Point_2(3,1);
  vA[5]= Point_2(4,3);
  vB[0]= Point_2(0,2);
  vB[1]= Point_2(0,0);
  vB[2]= Point_2(5,0);
  vB[3]= Point_2(5,2);
}

void test_result_output( const list<Object>& result ) {
  Point_2 pt;
  Segment_2 seg;
  Polygon_2 pgon;

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
    else { /* nothing */
      cout << "undefined object " << endl;
    }
  }
  cout << endl << endl << endl;
  
  return;
}

bool read_point(Point_2 &pt)
{
    long x, y, w;
    cin >> x >> y >> w;
    if (!cin.good())
	return false;
    pt = Point_2(TestNum(x), TestNum(y), TestNum(w));
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

#endif // EXAMPLE_IO_H
