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
// file          : demo/BooleanOperations/test_float_bops_wout.C
// source        : demo/BooleanOperations/test_float_bops_wout.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                        Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

//#define CGAL__BOPS_DEBUG_ON
//#define  CGAL__DCEL_DEBUG_ON

#include <CGAL/test_float_bops.h>

//#include <CGAL/IO/ostream.h>
//#include <fstream.h>
//#include <CGAL/IO/new_iostream.h>
#define red leda_red
#define black leda_black
#include <CGAL/IO/Window_stream.h>
#undef red 
#undef black

using namespace CGAL;
using namespace std;

Window_stream& operator <<(Window_stream& w, const Polygon_2& pgon);
Window_stream& operator <<(Window_stream& w, const Object& obj);

template <class T>
class Window_ostream_iterator : public output_iterator
{
protected:
    Window_stream* stream;
public:
    Window_ostream_iterator(Window_stream& s) : stream(&s) {}
    Window_ostream_iterator<T>& operator=(const T& value)
    {
        *stream << value;
        return *this;
    }
    Window_ostream_iterator<T>& operator*() { return *this; }
    Window_ostream_iterator<T>& operator++() { return *this; }
    Window_ostream_iterator<T>& operator++(int) { return *this; }
};


Window_stream& operator <<(Window_stream& w, const Polygon_2& pgon) {
  Window_ostream_iterator< Segment_2 > w_seg(w);
  copy(pgon.edges_begin(), pgon.edges_end(), w_seg);
  return w;
}


Window_stream& operator <<(Window_stream& w, const Object& obj) {
  Point_2 pt;
  Segment_2 seg;
  Polygon_2 pgon;
  if( assign( pgon, obj) )     { /* polygon */ w << pgon; }
  else if( assign( seg, obj) ) { /* segment */ w << seg; }
  else if( assign( pt, obj) )  { /* point */   w << pt; }
  return w;
}





void click_to_continue(Window_stream& W)
{
  double x, y;
  W.read_mouse(x,y);
}

#define wout (*W_global_ptr)
Window_stream *W_global_ptr;
Window_ostream_iterator< Object >   winout_obj(wout);
Window_ostream_iterator< Segment_2 > winout_seg(wout);


void test_result_win_output( const list<Object>& result ) {
  Point_2 pt;
  Segment_2 seg;
  Polygon_2 pgon;
 
  list<Object>::const_iterator it;
  wout << RED;
  for( it= result.begin(); it != result.end(); it++) {
    if( assign( pgon, *it) )     { /* polygon */ wout << pgon; }
    else if( assign( seg, *it) ) { /* segment */ wout << seg; }
    else if( assign( pt, *it) )  { /* point */   wout << pt; }
  }
 
  return;
}



int test_intersection(const Polygon_2& A, const Polygon_2& B) {
  cout << "INTERSECTION" << endl;
  wout << BLACK << A << BLUE << B;

  list<Object> result;
  intersection(A,B, back_inserter(result));
  test_result_output(result);
  //copy(result.begin(), result.end(), winout_obj);
  test_result_win_output( result );
  
  return 0;
}


int test_difference(const Polygon_2& A, const Polygon_2& B) {
  cout << "DIFFERENCE" << endl;
  wout << BLACK << A << BLUE << B;

  list<Object> result;
  difference(A,B, back_inserter(result));
  test_result_output(result);
  //copy(result.begin(), result.end(), winout_obj);
  test_result_win_output( result );
  
  return 0;
}


int test_union(const Polygon_2& A, const Polygon_2& B) {

  cout << "UNION" << endl;
  wout << BLACK << A << BLUE << B;

  list<Object> result;
  Union(A,B, back_inserter(result));
  test_result_output(result);
  //copy(result.begin(), result.end(), winout_obj);
  test_result_win_output( result );
  
  return 0;
}


int main( void )
{
  vector<Point_2> vA(6), vB(4);
  test_input( vA, vB);
  Polygon_2 A(vA.begin(), vA.end()), B(vB.begin(),vB.end());

  cout << "polygon A: " << A << endl;
  cout << "polygon B: " << B << endl;
 
  Window_stream W(600,600);
  W_global_ptr = &W;
  Bbox_2 box_A= A.bbox();
  Bbox_2 box_B= B.bbox();
  double xmin= min(box_A.xmin(), box_B.xmin());
  double xmax= max(box_A.xmax(), box_B.xmax());
  double ymin= min(box_A.ymin(), box_B.ymin());
  double ymax= max(box_A.ymax(), box_B.ymax());
  double dx=(xmax-xmin);
  double dy=(ymax-ymin);
  double d= max(dx,dy);
  dx *= 0.025;
  dy *= 0.025;
  W.init(xmin-dx,xmin+d+dx,ymin-dy);
  W.clear();

  test_intersection(A,B);
  click_to_continue(W);
  test_union(A, B);
  click_to_continue(W);
  test_difference(A, B);
  click_to_continue(W);
  test_difference(B, A);
  click_to_continue(W);
  return 0;
}

