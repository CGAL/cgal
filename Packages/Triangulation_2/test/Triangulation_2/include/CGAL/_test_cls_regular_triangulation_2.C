// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        : 
// file          : include/CGAL/_test_cls_regular_triangulation_2.C
// revision      : 
// revision_date : 

// author(s)     : Francois Rebufat, Mariette Yvinec

// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <iostream>
#include <fstream>

//#include <vector>

#include <CGAL/_test_fct_is_infinite.C>
#include <CGAL/_test_triangulation_iterators.C>
#include <CGAL/_test_triangulation_circulators.C>



template < class Triangulation, class Point, class Face_handle >
bool
_test_is_to_the_left( const Triangulation &T,
		      const Point &p,
		      const Face_handle &f,
		      const int li)
{
  return( T.orientation(f->vertex(f->ccw(li))->point(),
			f->vertex(f->cw(li))->point(),
			p)  == CGAL::LEFTTURN );
}

template <class Triangulation>
void 
_test_cls_reg_triangulation_2( const Triangulation & )
{
  typedef Triangulation                      Cls;

  // We assume the traits class has been tested already
  // actually, any traits is good if it has been tested
  typedef typename Cls::Geom_traits          Gt;

  typedef typename Cls::Point                Point;
  typedef typename Cls::Weighted_point       WPoint;
  typedef typename Cls::Segment              Segment;
  typedef typename Cls::Triangle             Triangle;

  typedef typename Cls::Vertex               Vertex;
  typedef typename Cls::Face                 Face;

  typedef typename Cls::Vertex_handle        Vertex_handle;
  typedef typename Cls::Face_handle          Face_handle;

  typedef std::pair<Face_handle,int>              Edge;

  typedef typename Cls::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Cls::Finite_faces_iterator    Finite_faces_iterator;
  typedef typename Cls::Finite_edges_iterator    Finite_edges_iterator;

  typedef typename Cls::Vertex_circulator    Vertex_circulator;
  typedef typename Cls::Face_circulator      Face_circulator;
  typedef typename Cls::Edge_circulator      Edge_circulator;
  typedef typename Cls::Line_face_circulator Line_face_circulator;

  typedef typename Cls::Locate_type          Locate_type;

  // Build a few objects
  // p1,p3,p2,p9,p8 aligned in this order
  // p0,p5,p7 also

  Point p0(5,6,1);
  Point p1(1,9,1);
  Point p2(6,14,1);
  Point p3(4,12,1);
  Point p4(3,29,1);
  Point p5(6,7,1);
  Point p6(6,39,1);
  Point p7(8,9,1);
  Point p8(10,18,1);
  Point p9(75625,155625,10000); // intersection of p2,p8 and p6,p7
  Point p10(10,50,2);
  Point p11(6,15,2); // midpoint p1,p0
  Point p12(6,16,2); // slightly above, in face
  Point p13(10,11,1);
  Point p14(10,40,1);
  Point p15(60,-10,1);

  int px=1, py=1;
  int qx=-1, qy=2;

  std::list<Point> l; l.push_back(p0);
  l.push_back(p1); l.push_back(p2); l.push_back(p3);
  l.push_back(p4); l.push_back(p5); l.push_back(p6);
  l.push_back(p7); l.push_back(p8); l.push_back(p9);
   
  std::vector<Point> v; v.push_back(p0);
  v.push_back(p1); v.push_back(p2); v.push_back(p3);
  v.push_back(p4); v.push_back(p5); v.push_back(p6);
  v.push_back(p7); v.push_back(p8); v.push_back(p9);
  
  
  WPoint wp0(p0,1);
  WPoint wp1(p1,1);
  WPoint wp2(p2,1);
  WPoint wp3(p3,5);
  WPoint wp4(p4,1);
  WPoint wp5(p5,8);
  WPoint wp6(p6,1);
  WPoint wp7(p7,1);
  WPoint wp8(p8,1);
  WPoint wp9(p9,3); // intersection of p2,p8 and p6,p7
  WPoint wp10(p10,2);
  WPoint wp11(p11,2); // midpoint p1,p0
  WPoint wp12(p12,2); // slightly above, in face
  WPoint wp13(p13,1);
  WPoint wp14(p14,1);
  WPoint wp15(p15,1);


  std::list<WPoint> lw; lw.push_back(wp0);
  lw.push_back(wp1); lw.push_back(wp2); lw.push_back(wp3);
 
   
  std::vector<WPoint> vw; vw.push_back(wp0);
  vw.push_back(wp1); vw.push_back(wp2); vw.push_back(wp3);
 

  /*****************************/
  /***** CONSTRUCTORS (1) ******/
  std::cout << "    constructors(1)" << std::endl;

  Cls T1;
  assert( T1.dimension() == -1 ); 
  assert( T1.number_of_vertices() == 0 );

  
  Cls T3(T1);
  Cls T4 = T1;
  T3.swap(T1);


  std::cout << "    insertions 0-dim" << std::endl;
  
  Cls T0_0;
  assert( T0_0.dimension() == -1 );
  assert( T0_0.number_of_vertices() == 0 );
  assert( T0_0.is_valid() );

  Cls T0_1; 
  Vertex_handle v0_1_0 = T0_1.insert(p0); assert( !v0_1_0.is_null() );
  assert( T0_1.dimension() == 0 );
  assert( T0_1.number_of_vertices() == 1 );
  assert( T0_1.is_valid() );

  Cls T0_2; 
  T0_2.insert_first(wp0);
  assert( T0_2.dimension() == 0 );
  assert( T0_2.number_of_vertices() == 1 );
  assert( T0_2.is_valid() );


  /******** 1-dimensional triangulations ******/
  // T1_n denotes a 1-dimensional triangulation with n vertices
  // when there are several, we use T1_n_p
  std::cout << "    insertions 1-dim" << std::endl;
  
  Cls T1_5;
  Vertex_handle v1_5_1 = T1_5.insert(wp1);
  Vertex_handle v1_5_2 = T1_5.insert(wp2);
  Vertex_handle v1_5_3 = T1_5.insert(wp3);
  Vertex_handle v1_5_8 = T1_5.insert(wp8);
  Vertex_handle v1_5_9 = T1_5.insert(wp9);
  assert( T1_5.dimension() == 1 );
  assert( T1_5.number_of_vertices() == 5 );
  // assert( T1_5.number_of_faces() == 0 );
  assert( T1_5.is_valid() );

  // test insert_second()
  Cls T1_6 = T0_2; 
  T1_6.insert_second( wp3);
  assert( T1_6.dimension() == 1 );
  assert( T1_6.number_of_vertices() == 2 );
  assert( T1_6.is_valid() ); 
  
  /******** 2-dimensional triangulations ******/ 
  std::cout << "    insertions 2-dim" << std::endl;
  Cls T2_3;
  Locate_type lt;
  Face_handle loc;
  int li;
  Vertex_handle v2_3_1 = T2_3.insert(wp1);
  Vertex_handle v2_3_2 = T2_3.insert(wp2);
  Vertex_handle v2_3_3 = T2_3.insert(wp3);
  Vertex_handle v2_3_8 = T2_3.insert(wp8);
  Vertex_handle v2_3_9 = T2_3.insert(wp9);
  assert( T2_3.dimension() == 1 );
  Vertex_handle v2_3_4 = T2_3.insert(wp4);
  assert( T2_3.dimension() == 2 );
  Vertex_handle v2_3_6 = T2_3.insert(wp6, T2_3.finite_faces_begin());
  Vertex_handle v2_3_0 = T2_3.insert(wp0, ++T2_3.finite_faces_begin());
  Vertex_handle v2_3_5 = T2_3.insert(wp5);
  Vertex_handle v2_3_7 = T2_3.insert(wp7);
  loc = T2_3.locate(wp10,lt,li);
  Vertex_handle v2_3_10 = T2_3.insert(wp10, lt, loc,li);
  assert( T2_3.dimension() == 2 );
  assert( T2_3.number_of_vertices() == 11 );
  assert( T2_3.is_valid() );
  
 
  // test list iterator insert
  Cls T2_5;
  assert( T2_5.insert(lw.begin(), lw.end()) == 4);
  assert( T2_5.dimension() == 2 );
  assert( T2_5.number_of_vertices() == 4 );
  assert( T2_5.is_valid() );

  // test list iterator insert
  Cls T2_6;
  assert( T2_6.insert(vw.begin(), vw.end()) == 4 );
  assert( T2_6.dimension() == 2 );
  assert( T2_6.number_of_vertices() == 4 );
  assert( T2_6.is_valid() );
  
  // test grid insert and make sure push_back exists
  Cls T2_7;
  int m, p;
  for (m=0; m<3; m++)
    for (p=0; p<3; p++)
      T2_7.push_back( WPoint(Point(m*px+p*qx, m*py+p*qy), 1) );
  assert( T2_7.number_of_vertices() == m*p );
  assert( T2_7.is_valid() );


 std::cout << "    constructors (2)" << std::endl;

  // test copy_constructor with non-empty 0-triangulation
  Cls T0_1_1( T0_1 );
  assert( T0_1_1.dimension() == 0 );
  assert( T0_1_1.number_of_vertices() == 1 );
  assert( T0_1_1.is_valid() );

  // test assignement
  Cls T0_1_2;
  T0_1_2 = T0_1;
  assert( T0_1_2.dimension() == 0 );
  assert( T0_1_2.number_of_vertices() == 1 );
  assert( T0_1_2.is_valid() );
  
  // test copy_constructor with non-empty 1-triangulation
  Cls T1_5_1( T1_5 );
  assert( T1_5_1.dimension() == 1 );
  assert( T1_5_1.number_of_vertices() == 5 );
  assert( T1_5_1.is_valid() );

  // Test assignment operator
  Cls T1_5_2 = T1_5;
  assert( T1_5_2.dimension() == 1 );
  assert( T1_5_2.number_of_vertices() == 5 );
  assert( T1_5_2.is_valid() );

   // test copy_constructor with non-empty 2-triangulation
  Cls T2_3_1( T2_3 );
  assert( T2_3_1.dimension() == 2 );
  assert( T2_3_1.number_of_vertices() == 11 );
  assert( T2_3_1.is_valid() );

  // test assignment operator
  Cls T2_3_4 = T2_3;
  assert( T2_3_4.dimension() == 2 );
  assert( T2_3_4.number_of_vertices() == 11 );
  assert( T2_3_4.is_valid() );
  
  /*********************************************/
  /****** FINITE/INFINITE VERTICES/FACES *******/

  std::cout << "    finite/infinite vertices/faces" << std::endl;
  _test_fct_is_infinite( T0_0 );
  _test_fct_is_infinite( T0_1 );
  _test_fct_is_infinite( T1_5 );
  _test_fct_is_infinite( T2_3 );




  /*************************************/
  /******** POINT LOCATIONS ************/

  // Check point location in 0-dimensional triangulations
  // No need because of precondition (at least two vertices)
  
  // Check point location in 1-dimensional triangulations
  std::cout << "    point locations 1-dim" << std::endl;
  Cls T1_3_2;
  T1_3_2.insert(p1);
  T1_3_2.insert(p2);
  T1_3_2.insert(p9); 
  loc = T1_3_2.locate(p1,lt,li); assert( lt == Cls::VERTEX );
  assert( T1_3_2.xy_equal(loc->vertex(li)->point(), p1) );
  loc = T1_3_2.locate(p2,lt,li); assert( lt == Cls::VERTEX );
  assert( T1_3_2.xy_equal(loc->vertex(li)->point(), p2) );
  loc = T1_3_2.locate(p9,lt,li); assert( lt == Cls::VERTEX );
  assert( T1_3_2.xy_equal(loc->vertex(li)->point(), p9) );
  loc = T1_3_2.locate(p3,lt,li); assert( lt == Cls::EDGE );
  assert( (T1_3_2.xy_equal(loc->vertex(loc->ccw(li))->point(), p1)
        && T1_3_2.xy_equal(loc->vertex(loc->cw(li))->point(), p2))
       || (T1_3_2.xy_equal(loc->vertex(loc->ccw(li))->point(), p2)
        && T1_3_2.xy_equal(loc->vertex(loc->cw(li))->point(), p1)));
  loc = T1_3_2.locate(p8,lt,li); assert( lt == Cls::OUTSIDE_CONVEX_HULL );
  loc = T1_3_2.locate(p7,lt,li); assert( lt == Cls::OUTSIDE_AFFINE_HULL );
  loc = T1_3_2.locate(p5,lt,li); assert( lt == Cls::OUTSIDE_AFFINE_HULL );
  loc = T1_3_2.locate(p4,lt,li); assert( lt == Cls::OUTSIDE_AFFINE_HULL );
  loc = T1_3_2.locate(p6,lt,li); assert( lt == Cls::OUTSIDE_AFFINE_HULL);
 

  // Check point location in 2-dimensional triangulations
  std::cout << "    point locations 2-dim" << std::endl;
  loc = T2_3.locate(p0,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), p0) );
  loc = T2_3.locate(p1,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), p1) );
  loc = T2_3.locate(p2,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), p2) );
  loc = T2_3.locate(p3,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), p3) );
  loc = T2_3.locate(p4,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), p4) );
  loc = T2_3.locate(p5,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), p5) );
  loc = T2_3.locate(p6,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), p6) );
  loc = T2_3.locate(p7,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), p7) );
  loc = T2_3.locate(p8,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), p8) );
  loc = T2_3.locate(p9,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), p9) );
  loc = T2_3.locate(p10,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), p10) );

  loc = T2_3.locate(p11,lt,li); assert( lt == Cls::EDGE);
  assert( (T2_3.xy_equal(loc->vertex(loc->ccw(li))->point(), p1)
        && T2_3.xy_equal(loc->vertex(loc->cw(li))->point(), p0))
       || (T2_3.xy_equal(loc->vertex(loc->ccw(li))->point(), p0)
        && T2_3.xy_equal(loc->vertex(loc->cw(li))->point(), p1)));
  loc = T2_3.locate(p12,lt,li); assert( lt == Cls::FACE );
  assert( T2_3.oriented_side(loc,p12) == CGAL::ON_POSITIVE_SIDE );
  loc = T2_3.locate(p13,lt,li,loc); assert( lt == Cls::OUTSIDE_CONVEX_HULL );
  li = loc->index(T2_3.infinite_vertex());
  assert( _test_is_to_the_left(T2_3,p13,loc,li) );
  loc = T2_3.locate(p14,lt,li); assert( lt == Cls::OUTSIDE_CONVEX_HULL );
  li = loc->index(T2_3.infinite_vertex());
  assert( _test_is_to_the_left(T2_3,p14,loc,li) );
  loc = T2_3.locate(p15,lt,li); assert( lt == Cls::OUTSIDE_CONVEX_HULL );
  li = loc->index(T2_3.infinite_vertex());
  assert( _test_is_to_the_left(T2_3,p15,loc,li) );


  /*************************/
  /******* Iterators *******/
  std::cout << "    iterators" << std::endl;
  // _test_iterators(T0_0);
  // _test_iterators(T0_1);
  _test_iterators(T1_5);
  _test_iterators(T2_3);


  /***************************/
  /******* Circulators *******/
  std::cout << "    circulators" << std::endl;
  _test_circulators(T1_5);
  _test_circulators(T2_3);
 


  // Line_face_circulator
  std::cout << "    line face circulator  " << std::endl;
  typedef typename Cls::Line_face_circulator LFC;
  // here == operator needed for Point!
  // testing with the grid triangulation
  LFC fc= T2_3.line_walk(p1,p10);
  assert(fc.ptr()!=NULL);
  assert(!fc.is_empty());
  LFC fc2=fc;
  assert(fc==fc2);
  fc++;
  fc--;
  ++fc;
  --fc;
  //Point pp(0,0.5,1);
  Point pp(0,1,2);
  loc = T2_3.locate(pp,lt,li);
  // assert(lt==Cls::VERTEX); no longer true, I've changed the triangulation
  fc= T2_3.line_walk(pp,p10,loc);
  fc2=fc;
  assert(fc==fc2);
  fc++;
  fc--;
  ++fc;
  --fc;
  // testing with dummy triangulations
  Cls T2_8;
     T2_8.insert(Point(0,0,1));
     T2_8.insert(Point(1,0,1));
     T2_8.insert(Point(1,1,1));
     T2_8.insert(Point(0,1,1));
  int n=0;
  do {fc2++ ; n = n+1;} while (fc2 != fc);
  assert(T2_8.number_of_vertices()>=2);
  assert(T2_8.is_valid());
  //fc= T2_8.line_walk(Point(0.5,0.4),Point(5,5));
  fc= T2_8.line_walk(Point(5,4,10),Point(5,5));
  fc2=fc;
  n=0;
  assert(fc==fc2);
  do {fc2++ ; n = n+1;} while (fc2 != fc);
  assert(n==4);
  // the two point are vertices of the triangulation.
  Cls TT;
  TT.insert(Point(0,0)); TT.insert(Point(1,0));
  TT.insert(Point(0,1)); TT.insert(Point(1,1));
  assert(TT.dimension()==2);
  assert(TT.is_valid());
  assert(TT.number_of_vertices()==4);
  loc = TT.locate(Point(0,0));
  fc = TT.line_walk(Point(0,0),Point(1,1));
  fc2 = TT.line_walk(Point(0,0),Point(1,1),loc);
  assert(fc==fc2);
  n=0;
  do {fc2++ ; n = n+1;} while (fc2 != fc);
  assert(n==3);

  /*****************************/
  /******** Miscellaneaous *****/
  std::cout << "    misc." << std::endl;
  assert( T0_0.ccw(0) == 1 );
  assert( T0_0.ccw(1) == 2 );
  assert( T0_0.ccw(2) == 0 );
  assert( T0_0.cw(0) == 2 );
  assert( T0_0.cw(1) == 0 );
  assert( T0_0.cw(2) == 1 );

  // the assert() are to avoid compiler warnings about unused variables
  loc = T2_3.locate(p12,lt,li); // from section locate above
  Triangle t = T2_3.triangle(loc); assert( &t == &t );
  Segment  s = T2_3.segment(loc,0); assert( &s == &s );
  s = T2_3.segment(Edge(loc,1)); assert( &s == &s );
  s = T2_3.segment(v2_3_6->incident_edges()); assert( &s == &s );
  s = T2_3.segment(T2_3.finite_edges_begin()); assert( &s == &s );

  /********************/
  /***** Duality ******/
  std::cout << "    duality" << std::endl;
   _test_regular_duality(T1_5);
   _test_regular_duality(T2_3);


  /********************/
  /******** I/O *******/
  std::cout << "    output to a file" << std::endl;
  std::ofstream of1_5("T15.triangulation");
  CGAL::set_ascii_mode(of1_5); 
  of1_5 << T1_5; of1_5.close();
  std::ofstream of2_3("T23.triangulation");
  CGAL::set_ascii_mode(of2_3); 
  of2_3 << T2_3; of2_3.close();
  

  std::cout << "    input from a file" << std::endl;
  
  std::ifstream if1_5("T15.triangulation"); CGAL::set_ascii_mode(if1_5); 
  Cls T1_5_copy; if1_5 >> T1_5_copy;
  assert( T1_5_copy.is_valid() &&
	  T1_5_copy.number_of_vertices() == T1_5.number_of_vertices() );
  
  std::ifstream if2_3("T23.triangulation"); CGAL::set_ascii_mode(if2_3);
  Cls T2_3_copy; if2_3 >> T2_3_copy;
  assert( T2_3_copy.is_valid() &&
	  T2_3_copy.number_of_vertices() == T2_3.number_of_vertices() );
  
  /**********************/
  /***** REMOVALS *******/ 
  std::cout << "    removals" << std::endl;

  // test remove_first()
  T0_1.remove_first(T0_1.finite_vertex());
  assert( T0_1.number_of_vertices() == 0 );

   // remove from 1-dimensional triangulations
  T1_5.remove(v1_5_1);
  T1_5.remove(v1_5_2);
  T1_5.remove(v1_5_3);
  T1_5.remove(v1_5_8);
  T1_5.remove(v1_5_9);
  assert( T1_5.number_of_vertices() == 0 );

  // remove from 2-dimensional triangulations
  T2_3.remove(v2_3_0);
  T2_3.remove(v2_3_1);
  T2_3.remove(v2_3_9);
  T2_3.remove(v2_3_8);
  T2_3.remove(v2_3_5);
  T2_3.remove(v2_3_3);
  T2_3.remove(v2_3_4);
  T2_3.remove(v2_3_2);
  T2_3.remove(v2_3_6);
  T2_3.remove(v2_3_7);
  T2_3.remove(v2_3_10);
  assert( T2_3.number_of_vertices() == 0 );

  int i;
  T2_5.clear();
  assert( T2_5.number_of_vertices() == 0 );

  for (i=T2_6.number_of_vertices(); i>0; i--)
    T2_6.remove(T2_6.finite_vertex());
  assert( T2_6.number_of_vertices() == 0 );
  


  // test destructors and return
  std::cout << "    test destructors and return" << std::endl;
}


template <class Del>
void
_test_regular_duality( const Del &T )
{
  typedef typename Del::Geom_traits          Gt;
  typedef typename Del::Finite_faces_iterator        Face_iterator;
  typedef typename Del::Finite_edges_iterator        Edge_iterator;
  typedef typename Del::Edge_circulator              Edge_circulator;

  // Test dual(face iterator)
  //dual of faces is tested via dual of edges
  
  // Test dual(edge iterator)
  Edge_iterator eit;
  for (eit =  T.finite_edges_begin(); eit !=  T.finite_edges_end(); ++eit)
    {
      CGAL::Object o = T.dual(eit);
      typename Gt::Ray_2 r;
      typename Gt::Segment_2 s;
      typename Gt::Line_2 l;
      if ( CGAL::assign(s,o) ) {
        assert(  ! T.is_infinite((*eit).first) );
	assert( ! T.is_infinite(((*eit).first)->neighbor((*eit).second )) );
      } 
      else if ( CGAL::assign(l,o) ) {
        assert( T.dimension() == 1 );
      } 
      else {
        assert( CGAL::assign(r,o) );
      }
    }

  // Test dual(edge circulator)
  Edge_circulator ec=T.finite_vertices_begin()->incident_edges(), done(ec);
  if ( !ec.is_empty() ) 
  do  
    {
      if (! T.is_infinite(ec)){
	CGAL::Object o = T.dual(ec);
	typename Gt::Ray_2 r;
        typename Gt::Segment_2 s;
	typename Gt::Line_2 l;
	assert( CGAL::assign(s,o) || CGAL::assign(r,o) || CGAL::assign(l,o) );
      }
      ++ec;
    } while ( ec == done);
}
