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
// file          : demo/BooleanOperations/test_bops_V2E.C
// source        : demo/BooleanOperations/test_bops_V2E.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                  Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

#ifdef __GNUC__
#include <typeinfo>
#endif
#include <cstdlib>

#define DEBUG_ON
#include <CGAL/bops_V2E_rep.h>
#include <CGAL/bops_V2E_rep_out.h>

typedef vector<char>::iterator vertex;
typedef vector<pair<int,int> >::iterator edge;

  int ary_x[]= { 0, 0, 1, 1};
  int ary_y[]= { 0, 1, 0, 1};

struct bops_compare_vertices {
  vertex v0;
  bops_compare_vertices() {}
  bops_compare_vertices(vertex _v0) : v0(_v0) {}


  int xcoord(char c) const { return ary_x[c-'a']; }
  int ycoord(char c) const { return ary_y[c-'a']; }

  bool compare(vertex v1, vertex v2) const {
    cout << " v0= " << (*v0)
         << " v1= " << (*v1)
         << " v2= " << (*v2);

    int result= 
      xcoord(*v0) * ycoord(*v1) - ycoord(*v0) * xcoord(*v1) +
      ycoord(*v0) * xcoord(*v2) - xcoord(*v0) * ycoord(*v2) +
      xcoord(*v1) * ycoord(*v2) - ycoord(*v1) * xcoord(*v2);
    cout << "   : " << result << endl;

    return result >0; 
  }
  bool operator()(vertex v1, vertex v2) const {
    return compare(v1, v2);
  }
};


int test0( bool four_points ) {

  vector< char > V;
  vector< pair<int,int> > E;

  /* vertices: 0, 1, 2, 3
     edges:    0(0,1), 1(1,2), 2(2,3), 3(3,0), 4(0,2)
  int four_points= 0;
  if( argc > 1) {
    cout << "FOUR-POINTS !!" << endl << endl;
    four_points= atoi(argv[1]);
  }
  */

  if( four_points ) {
    V.push_back('a');
    V.push_back('b');
    V.push_back('c');
    V.push_back('d');

    E.push_back( make_pair(0, 1) );
    E.push_back( make_pair(1, 2) );
    E.push_back( make_pair(2, 3) );
    E.push_back( make_pair(3, 0) );
    E.push_back( make_pair(0, 2) ); /* diagonale */
  }
  else {
    V.push_back('a');
    V.push_back('b');
    V.push_back('c');
    V.push_back('d');
    V.push_back('e');

    E.push_back( make_pair(0, 1) );
    E.push_back( make_pair(1, 2) );
    E.push_back( make_pair(2, 3) );
    E.push_back( make_pair(3, 4) );
    E.push_back( make_pair(4, 0) );
    E.push_back( make_pair(4, 1) ); /* diagonale */
  }
  
  int n= V.size();
  int m= E.size();
  _V2E_rep_type< vertex, edge, bops_compare_vertices> v2e(n,m);

  edge e= E.begin();
  vertex v1, v2;
  v1= V.begin(); v2= v1+1;
  v2e.insert( v1-V.begin(), v1, v2-V.begin(), v2, e);
  e++; v1++, v2++;
  v2e.insert( v1-V.begin(), v1, v2-V.begin(), v2, e);
  e++; v1++, v2++;
  v2e.insert( v1-V.begin(), v1, v2-V.begin(), v2, e);
  e++; v1++, v2++;
  if( v2 == V.end() ) v2= V.begin();
  v2e.insert( v1-V.begin(), v1, v2-V.begin(), v2, e);
  e++; v1= V.begin()+0; v2= V.begin()+2;
  v2e.insert( v1-V.begin(), v1, v2-V.begin(), v2, e);

  cout << (_V2E_rep_base_type<vertex, edge>)v2e;
  v2e.sort_vertices_CCW();
  cout << (_V2E_rep_base_type<vertex, edge>)v2e;

  return 0;
}

template<class NT>
struct point_2 {
  point_2() {}
  point_2(const NT& cx, const NT& cy) : _x(cx), _y(cy) {}
  point_2(const point_2& pt) : _x(pt.x()), _y(pt.y()) {}
  point_2& operator=(const point_2& pt) {
    _x= pt._x; _y= pt._y; return *this;
  } 
  const NT& x() const { return _x; }
  const NT& y() const { return _y; }
  NT& x() { return _x; }
  NT& y() { return _y; }
  NT _x, _y;
};

template<class NT>
ostream& operator <<( ostream& o, const point_2<NT>& pt) {
  o << '(' << pt.x() << ", " << pt.y() << ')';
  return o;
}

typedef point_2<double> point;
typedef vector<point>::iterator   pt_vertex;


struct bops_compare_points {
  bops_compare_points() {}
  bops_compare_points(pt_vertex v) : v0(v) {}
  pt_vertex v0;
  double xcoord( pt_vertex v ) const { return (*v).x(); }
  double ycoord( pt_vertex v ) const { return (*v).y(); }
 
  bool compare( pt_vertex v1, pt_vertex v2) const {
    /*  cout << " v0= " << (*v0) << " v1= " << (*v1)
         << " v2= " << (*v2);
	 */

    double result= 
      xcoord(v0) * ycoord(v1) - ycoord(v0) * xcoord(v1) +
      ycoord(v0) * xcoord(v2) - xcoord(v0) * ycoord(v2) +
      xcoord(v1) * ycoord(v2) - ycoord(v1) * xcoord(v2);
    cout << "   : " << result << endl;

    return result > 0.0; 
  }

  bool operator()(pt_vertex v1, pt_vertex v2) const {
    return compare(v1,v2);
  }
};

int test1( void ) {

  vector< point > V(16);
  vector< pair<int,int> > E(24);

  V[0]= point(4.15, 4.15);
  V[1]= point(-0.2, 6.55);
  V[2]= point(-5.05, 1.15);
  V[3]= point(-4.65, -5.8);
  V[4]= point(3.05, -5.15);
  V[5]= point(5.85, 0.3);
  V[6]= point(8.1, -2.55);
  V[7]= point(1.75, 1.55);
  V[8]= point(-3.5, 4.6);
  V[9]= point(-6.2, -2.15);
  V[10]= point(-2.48237, 4.00881); 
  V[11]= point(-4.74349, 1.49126); 
  V[12]= point(-5.02862, 0.778459);
  V[13]= point(-4.85791, -2.18754);
  V[14]= point(5.31096, -0.749203);
  V[15]= point(4.4384, -2.44758);

  E[0]=  make_pair(0, 1);
  E[1]=  make_pair(3, 4);
  E[2]=  make_pair(0, 5);
  E[3]=  make_pair(2, 11);
  E[4]=  make_pair(10, 11);
  E[5]=  make_pair(1, 10);
  E[6]=  make_pair(8, 10);
  E[7]=  make_pair(7, 10);
  E[8]=  make_pair(9, 12);
  E[9]=  make_pair(11, 12);
  E[10]=  make_pair(8, 11);
  E[11]=  make_pair(2, 12);
  E[12]=  make_pair(12, 13);
  E[13]=  make_pair(3, 13);
  E[14]=  make_pair(9, 13);
  E[15]=  make_pair(13, 15);
  E[16]=  make_pair(6, 15);
  E[17]=  make_pair(4, 15);
  E[18]=  make_pair(14, 15);
  E[19]=  make_pair(5, 14);
  E[20]=  make_pair(7, 14);
 
  
  int n= V.size();
  int m= E.size();
  _V2E_rep_type< pt_vertex, edge, bops_compare_points> v2e(n,m);

  for( edge e= E.begin(); e != E.end(); e++) {
    int i1= (*e).first;
    int i2= (*e).second;
    v2e.insert( i1, &V[i1], i2, &V[i2], e);
  }

  cout << (_V2E_rep_base_type<pt_vertex, edge>)v2e;
  v2e.sort_vertices_CCW();
  cout << (_V2E_rep_base_type<pt_vertex, edge>)v2e;

  return 0;
}

int main( int argc , char *argv[] ) {

  test0( ((bool)argc > 1) );
  test1();

  return 0;
}
