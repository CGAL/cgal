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
// file          : demo/BooleanOperations/include/CGAL/test_bops_cin_wout.C
// source        : demo/BooleanOperations/include/CGAL/test_bops_cin_wout.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                        Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

//#include <CGAL/test_bops.h>

//#include <CGAL/IO/ostream.h>
//#include <fstream.h>
//#include <CGAL/IO/new_iostream.h>
#define red leda_red
#define black leda_black
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/polygon_Window_stream.h>
#undef red 
#undef black


//Window_stream& operator <<(Window_stream& w, const Polygon_2& pgon);
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


Window_stream& operator <<(
		Window_stream& w,
		const Polygon_2& pgon)
{
	Window_ostream_iterator< Segment_2 > w_seg(w);
	copy(pgon.edges_begin(), pgon.edges_end(), w_seg);
	return w;
}

void print_stream(Window_stream& w, const Polygon_2& pgon)
{
	Window_ostream_iterator< Segment_2 > w_seg(w);
	copy(pgon.edges_begin(), pgon.edges_end(), w_seg);
	return;
}                     

Window_stream& operator <<(Window_stream& w, const Object& obj) {
  Point_2 pt;
  Segment_2 seg;
  Polygon_2 pgon;
  Iso_rectangle_2 irect;
  Triangle_2 triangle;
  if( assign( pgon, obj) )          {
	/* polygon */
	print_stream(w,pgon);
	//w << pgon;
  }
  else if( assign( irect, obj) )    { /* iso-rectangle */ w << irect; }
  else if( assign( triangle, obj) ) { /* triangle */ w << triangle; }
  else if( assign( seg, obj) )      { /* segment */ w << seg; }
  else if( assign( pt, obj) )       { /* point */   w << pt; }
  return w;
}

#define wout (*W_global_ptr)
Window_stream *W_global_ptr;
Window_ostream_iterator< Object >   winout_obj(wout);
Window_ostream_iterator< Segment_2 > winout_seg(wout);

int click_to_continue(Window_stream& W= wout)
{
  double x, y;
  return W.read_mouse(x,y);
}


void test_result_win_output_base( const list<Object>& result ) {
  Point_2 pt;
  Segment_2 seg;
  Polygon_2 pgon;
  Iso_rectangle_2 irect;
  Triangle_2 triangle;
 
  list<Object>::const_iterator it;
  wout << RED;
  for( it= result.begin(); it != result.end(); it++) {
    if( assign( pgon, *it) )          { /* polygon */ wout << pgon; }
    else if( assign( irect, *it) )    { /* iso-rectangle */ wout << irect; }
    else if( assign( triangle, *it) ) { /* triangle */ wout << triangle; }
    else if( assign( seg, *it) )      { /* segment */ wout << seg; }
    else if( assign( pt, *it) )       { /* point */   wout << pt; }
  }
 
  return;
}

template< class _polygon>
void test_result_win_output(
   const _polygon& A,
   const _polygon& B,
   const list<Object>& result
)
{
  bool loop_on= true;
  wout.clear();
  wout << BLACK << A << BLUE << B;
  test_result_win_output_base(result);
  while( loop_on == true ) {
    switch( click_to_continue() ) {
      case -2: //middle button
        wout.clear();
        wout << BLACK << A << BLUE << B;
        break;
      case -3: //right button
        wout.clear();
        wout << BLACK << A << BLUE << B;
        test_result_win_output_base(result);
        break;
      case -1: //left button
      default:
        loop_on= false;
        //cout << "button= " << button << endl;
    }
  }
  return;
}


/*
    --------------------------------------------------------------------
       Polygons
    --------------------------------------------------------------------
*/

template< class _polygon>
void test_intersection(const _polygon& A, const _polygon& B) {
  wout << BLACK << A << BLUE << B;
  list<Object> result;
  intersection(A,B, back_inserter(result));
  test_result_output(result);
  test_result_win_output(A,B,result);
}


template< class _polygon>
void test_difference(const _polygon& A, const _polygon& B) {
  wout << BLACK << A << BLUE << B;
  list<Object> result;
  difference(A,B, back_inserter(result));
  test_result_output(result);
  test_result_win_output(A,B,result);
}


template< class _polygon>
void test_union(const _polygon& A, const _polygon& B) {
  wout << BLACK << A << BLUE << B;
  list<Object> result;
  union(A,B, back_inserter(result));
  test_result_output(result);
  test_result_win_output(A,B,result);
}



/*
    --------------------------------------------------------------------
       MAIN
    --------------------------------------------------------------------
*/

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
  cout << "polygon A: (BLACK) " << A << endl;

  if (!read_pgn(B)) {
	cerr << "Polygon B is wrong!" << endl;
	return 1;
  }
  // check counterclockwise orientation
  if ( B.orientation() == CLOCKWISE) {
    B.reverse_orientation();
    cout << "orientation reversed" << endl;
  }
  cout << "polygon B: (BLUE) " << B << endl;
  
  // WINDOW-INIT
  float w = 600.0, h = 600.0;
  Window_stream W(w,h);
  cgalize(W);
  //double x_extension= 1.0;
  //double y_extension = x_extension * h / w;
  //W.init(-x_extension, x_extension, -y_extension);//(-1.0, 1.0, -1.0);
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
  W.display();
  W.clear();
  wout << BLACK << A << BLUE << B;
  click_to_continue(W);

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
  case 5:  // iso-rectangles
    //Iso_rectangle_2 rA, rB;
    if( A.size() != 4 )
      cout << "size of polygon A is NOT 4 (=" << A.size() << ")" <<endl;
    else if( B.size() != 4 )
      cout << "size of polygon B is NOT 4 (=" << B.size() << ")" <<endl;
    else {
      Polygon_2::Vertex_const_iterator it;
      it= A.vertices_begin();
      Iso_rectangle_2 rA(*it++, *++it);
      it= B.vertices_begin();
      Iso_rectangle_2 rB(*it++, *++it);
      cout << "ISO-RECT. INTERSECTION: " << endl;
      test_intersection(rA,rB);
      cout << "ISO-RECT. UNION: " << endl;
      test_union(rA,rB);
      cout << "ISO-RECT. DIFFERENCE: A-B" << endl;
      test_difference(rA,rB);
      cout << "ISO-RECT. DIFFERENCE: B-A" << endl;
      test_difference(rB,rA);
    }
    
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
