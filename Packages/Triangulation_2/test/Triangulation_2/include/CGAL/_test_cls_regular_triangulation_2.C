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
// file          : include/CGAL/_test_cls_triangulation_2.C
// revision      : 
// revision_date : 

// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)

// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <iostream>
#include <fstream>

//#include <pair.h>

#include <vector.h>

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
  return( T.geom_traits().orientation
              (f->vertex(f->ccw(li))->point(),
	       f->vertex(f->cw(li))->point(),
	       p)
	  == CGAL::LEFTTURN );
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

  typedef typename Cls::Vertex_iterator      Vertex_iterator;
  typedef typename Cls::Face_iterator        Face_iterator;
  typedef typename Cls::Edge_iterator        Edge_iterator;

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
  assert( Gt().compare(p11,p11) );
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
  assert( Gt().compare(wp11,wp11) );
  WPoint wp12(p12,2); // slightly above, in face
  WPoint wp13(p13,1);
  WPoint wp14(p14,1);
  WPoint wp15(p15,1);


  std::list<WPoint> lw; lw.push_back(wp0);
  lw.push_back(wp1); lw.push_back(wp2); lw.push_back(wp3);
  lw.push_back(wp4); lw.push_back(wp5); lw.push_back(wp6);
  lw.push_back(wp7); lw.push_back(wp8); lw.push_back(wp9);
   
  std::vector<WPoint> vw; vw.push_back(wp0);
  vw.push_back(wp1); vw.push_back(wp2); vw.push_back(wp3);
  vw.push_back(wp4); vw.push_back(wp5); vw.push_back(wp6);
  vw.push_back(wp7); vw.push_back(wp8); vw.push_back(wp9);


  /*****************************/
  /***** CONSTRUCTORS (1) ******/
  cout << "    constructors(1)" << endl;

  Cls T1;
  assert( T1.dimension() == 0 ); // should be -1, but we'll discuss this
  assert( T1.number_of_vertices() == 0 );

  // Cls T2(Gt()), T3(T2);
  Cls T3(T1);
  // assert( T3.geom_traits() == T2.geom_traits() ); // assert identity of traits

  Cls T4 = T1;
  // assert( T4.geom_traits() == T1.geom_traits() ); // assert identity of traits
  
  T3.swap(T1);


  cout << "    insertions 0-dim" << endl;
  
  Cls T0_0;
  assert( T0_0.dimension() == 0 );
  assert( T0_0.number_of_vertices() == 0 );
  assert( T0_0.is_valid() );

  Cls T0_1; 
  Vertex_handle v0_1_0 = T0_1.insert(p0); assert( !v0_1_0.is_null() );
  assert( T0_1.dimension() == 0 );
  assert( T0_1.number_of_vertices() == 1 );
  assert( T0_1.is_valid() );

Cls T0_2; 
  // T0_2.insert_first(Vertex(p0).handle());
  // this  statement cause a segmentation fault on Linux
  // when the whole procedure is leaved
  Vertex_handle v0_2_0(new Vertex(wp0));
  T0_2.insert_first(v0_2_0);

  assert( T0_2.dimension() == 0 );
  assert( T0_2.number_of_vertices() == 1 );
  assert( T0_2.is_valid() );


  /******** 1-dimensional triangulations ******/
  // T1_n denotes a 1-dimensional triangulation with n vertices
  // when there are several, we use T1_n_p
  cout << "    insertions 1-dim" << endl;
  
  Cls T1_2;
  Vertex_handle v1_2_1 = T1_2.insert(wp1);
  Vertex_handle v1_2_2 = T1_2.insert(wp2);
  assert( T1_2.dimension() == 1 );
  assert( T1_2.number_of_vertices() == 2 );
  // assert( T1_2.number_of_faces() == 0 );
  assert( T1_2.is_valid() );
  
  // p1,p2,p3  [endpoints first]
  Cls T1_3_0;
  Vertex_handle v1_3_0_1 = T1_3_0.insert(wp1); assert( !v1_3_0_1.is_null() );
  Vertex_handle v1_3_0_3 = T1_3_0.insert(wp3); assert( !v1_3_0_3.is_null() );
  Vertex_handle v1_3_0_2 = T1_3_0.insert(wp2); assert( !v1_3_0_2.is_null() );
  assert( T1_3_0.dimension() == 1 );
  assert( T1_3_0.number_of_vertices() == 3 );
  // assert( T1_3_0.number_of_faces() == 0 );
  assert( T1_3_0.is_valid() );
  
  // p1,p2,p3  [middle point first]
  Cls T1_3_1;
  Vertex_handle v1_3_1_1 = T1_3_1.insert(wp1); assert( !v1_3_1_1.is_null() );
  Vertex_handle v1_3_1_3 = T1_3_1.insert(wp3); assert( !v1_3_1_3.is_null() );
  Vertex_handle v1_3_1_2 = T1_3_1.insert(wp2); assert( !v1_3_1_2.is_null() );
  assert( T1_3_1.dimension() == 1 );
  assert( T1_3_1.number_of_vertices() == 3 );
  // assert( T1_3_1.number_of_faces() == 0 );
  // assert( T1_3_0 == T1_3_1 ); // operator== is not defined!
  assert( T1_3_1.is_valid() );

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
  //T1_6.insert_second(Vertex(p3).handle());
  // the following statement cause a segmentation fault on Linux
  // when the whole procedure is leaved
  T1_6.insert_second( Vertex_handle(new Vertex(wp3)));
  assert( T1_6.dimension() == 1 );
  assert( T1_6.number_of_vertices() == 2 );
  assert( T1_6.is_valid() ); 
  
  /******** 2-dimensional triangulations ******/ 
  cout << "    insertions 2-dim" << endl;
  
  Cls T2_1;
  Vertex_handle v2_1_0 = T2_1.insert(wp0);
  Vertex_handle v2_1_1 = T2_1.insert(wp1);
  Vertex_handle v2_1_2 = T2_1.insert(wp2);
  Vertex_handle v2_1_3 = T2_1.insert(wp3);    // on the edge p1,p2, on the convex hull
  Vertex_handle v2_1_4 = T2_1.insert(wp4);    // outside, with two visible collinear edges
  Vertex_handle v2_1_5 = T2_1.insert(wp5); 
  Vertex_handle v2_1_6 = T2_1.insert(wp6);    // outside, collinear with p2,p5
  Vertex_handle v2_1_7 = T2_1.insert(wp7);    // outside with two visible collinear edges
                      		             // but also collinear with and extending p0,p5 
  Vertex_handle v2_1_8 = T2_1.insert(wp8); 
  Vertex_handle v2_1_9 = T2_1.insert(wp9);    // inside, on the edge p6,p7
  Vertex_handle v2_1_10 = T2_1.insert(wp10);  // inside the face p2,p4,p6
  assert( T2_1.dimension() == 2 );
  assert( T2_1.number_of_vertices() == 11 );
  
  // test is_valid for 2-triangulations
  assert( T2_1.is_valid() );

  // we now test the other insert functions
  // more vicious, we insert all the points on a single line first
  Cls T2_3;
  Locate_type lt;
  Vertex_handle v2_3_1 = T2_3.insert(wp1);
  Vertex_handle v2_3_2 = T2_3.insert(wp2);
  Vertex_handle v2_3_3 = T2_3.insert(wp3, lt);
  assert( lt == Cls::EDGE );
  Vertex_handle v2_3_8 = T2_3.insert(wp8, lt);
  assert( lt == Cls::COLLINEAR_OUTSIDE );
  Vertex_handle v2_3_9 = T2_3.insert(wp9);
  assert( T2_3.dimension() == 1 );
  Vertex_handle v2_3_4 = T2_3.insert(wp4);
  assert( T2_3.dimension() == 2 );
  Vertex_handle v2_3_6 = T2_3.insert(wp6, T2_3.faces_begin());
  Vertex_handle v2_3_0 = T2_3.insert(wp0, lt, ++T2_3.faces_begin());
  assert( lt == Cls::OUTSIDE );
  Vertex_handle v2_3_5 = T2_3.insert(wp5);
  Vertex_handle v2_3_7 = T2_3.insert(wp7);
  Vertex_handle v2_3_10 = T2_3.insert(wp10, lt, ++(++(T2_3.faces_begin())));
  assert( lt == Cls::FACE );
  assert( T2_3.dimension() == 2 );
  assert( T2_3.number_of_vertices() == 11 );
  // assert( T2_3.number_of_faces() == 13 );
  assert( T2_3.is_valid() );
  
  // make sure inserting on a previous point does not insert it again
  assert( T2_3.insert(wp10, lt) == v2_3_10 );


  //  assert( lt == Cls::VERTEX );


  assert( T2_3.number_of_vertices() == 11 );

  // make sure push_back exists and does the same thing as insert
  assert( T2_3.push_back(wp10) == v2_3_10 );
  assert( T2_3.number_of_vertices() == 11 );


  // test list iterator insert
  Cls T2_5;
  assert( T2_5.insert(lw.begin(), lw.end()) == 10 );
  assert( T2_5.dimension() == 2 );
  assert( T2_5.number_of_vertices() == 10 );
  // assert( T2_5.number_of_faces() == 13 );
  assert( T2_5.is_valid() );

  // test list iterator insert
  Cls T2_6;
  assert( T2_6.insert(vw.begin(), vw.end()) == 10 );
  assert( T2_6.dimension() == 2 );
  assert( T2_6.number_of_vertices() == 10 );
  // assert( T2_6.number_of_faces() == 13 );
  assert( T2_6.is_valid() );
  
  // test grid insert
  Cls T2_7;
  int m, p;
  for (m=0; m<20; m++)
    for (p=0; p<20; p++)
      T2_7.insert( WPoint(Point(m*px+p*qx, m*py+p*qy), 1) );
  assert( T2_7.number_of_vertices() == m*p );
  assert( T2_7.is_valid() );


 cout << "    constructors (2)" << endl;

  // test copy_constructor with non-empty 0-triangulation
  Cls T0_1_1( T0_1 );
  assert( T0_1_1.dimension() == 0 );
  assert( T0_1_1.number_of_vertices() == 1 );
  assert( T0_1_1.is_valid() );

  // test constructor that takes a vertex handle
  // this must(!) be the vertex at infinity
  Cls T0_1_2( T0_1.infinite_vertex() , T0_1.geom_traits());
  assert( T0_1_2.dimension() == 0 );
  assert( T0_1_2.number_of_vertices() == 1 );
  assert( T0_1_2.is_valid() );
  // copy the triangulation to avoid having two triangulations
  // with the same set of vertices and faces
  // which causes a segmentation fault when the 2d one is deleted
  T0_1_2 = T0_1;
  
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

  // test constructor that takes a vertex handle
  // this must(!) be the vertex at infinity
  Cls T1_5_3( T1_5_2.infinite_vertex() , T1_5_2.geom_traits());
  assert( T1_5_3.dimension() == 1 );
  assert( T1_5_3.number_of_vertices() == 5 );
  assert( T1_5_3.is_valid() );
  // copy the triangulation to avoid having two triangulations
  // with the same set of vertices and faces
  // which causes a segmentation fault when the 2d one is deleted
  T1_5_3 = T1_5_2;
  
  // test constructor that takes a vertex handle and a geom_traits
  Cls T1_5_4( T1_5_2.infinite_vertex(), T2_1.geom_traits() );
  assert( T1_5_4.dimension() == 1 );
  assert( T1_5_4.number_of_vertices() == 5 );
  assert( T1_5_4.is_valid() );
  // Exchange the triangulation -- same reason
  T1_5_4 = T1_5_2;
 
  // test copy_constructor with non-empty 2-triangulation
  Cls T2_1_1( T2_1 );
  assert( T2_1_1.dimension() == 2 );
  assert( T2_1_1.number_of_vertices() == 11 );
  assert( T2_1_1.is_valid() );
  
  // test constructor that takes a vertex handle
  Cls T2_1_2( T2_1.infinite_vertex(), T2_1.geom_traits() );
  assert( T2_1_2.dimension() == 2 );
  assert( T2_1_2.number_of_vertices() == 11 );
  assert( T2_1_2.is_valid() );
  // copy the triangulation to avoid having two triangulations
  // with the same set of vertices and faces
  // which causes a segmentation fault when the 2d one is deleted
  T2_1_2 = T2_1;

  // test constructor that takes a vertex handle and a geom_traits
  Cls T2_1_3( T2_1.infinite_vertex(), T2_1.geom_traits() );
  assert( T2_1_3.dimension() == 2 );
  assert( T2_1_3.number_of_vertices() == 11 );
  assert( T2_1_3.is_valid() );
  // Copy the  triangulation, see T1_5_3
  T2_1_3 = T2_1;

  // test assignment operator
  Cls T2_1_4 = T2_1;
  assert( T2_1_4.dimension() == 2 );
  assert( T2_1_4.number_of_vertices() == 11 );
  assert( T2_1_4.is_valid() );
  
  /*********************************************/
  /****** FINITE/INFINITE VERTICES/FACES *******/

  cout << "    finite/infinite vertices/faces" << endl;
  _test_fct_is_infinite( T0_0 );
  _test_fct_is_infinite( T0_1 );
  _test_fct_is_infinite( T1_2 );
  _test_fct_is_infinite( T1_5 );
  _test_fct_is_infinite( T2_1 );
  _test_fct_is_infinite( T2_3 );
  _test_fct_is_infinite( T2_5 );
  _test_fct_is_infinite( T2_6 );



  /*************************************/
  /******** POINT LOCATIONS ************/

  // Locate_type lt; // see above
  int            li;
  Face_handle    f;

  // Check point location in 0-dimensional triangulations
  // No need because of precondition (at least two vertices)
  
  // Check point location in 1-dimensional triangulations
  cout << "    point locations 1-dim" << endl;
  Cls T1_3_2;
  T1_3_2.insert(p1);
  T1_3_2.insert(p2);
  T1_3_2.insert(p9); 
  f = T1_3_2.locate(p1,lt,li); assert( lt == Cls::VERTEX );
  assert( T1_3_2.geom_traits().compare(f->vertex(li)->point(), p1) );
  f = T1_3_2.locate(p2,lt,li); assert( lt == Cls::VERTEX );
  assert( T1_3_2.geom_traits().compare(f->vertex(li)->point(), p2) );
  f = T1_3_2.locate(p9,lt,li); assert( lt == Cls::VERTEX );
  assert( T1_3_2.geom_traits().compare(f->vertex(li)->point(), p9) );
  f = T1_3_2.locate(p3,lt,li); assert( lt == Cls::EDGE );
  assert( (T1_3_2.geom_traits().compare(f->vertex(f->ccw(li))->point(), p1)
        && T1_3_2.geom_traits().compare(f->vertex(f->cw(li))->point(), p2))
       || (T1_3_2.geom_traits().compare(f->vertex(f->ccw(li))->point(), p2)
        && T1_3_2.geom_traits().compare(f->vertex(f->cw(li))->point(), p1)));
  f = T1_3_2.locate(p8,lt,li); assert( lt == Cls::COLLINEAR_OUTSIDE );
  assert( T1_3_2.geom_traits().compare(f->vertex(li)->point(), p9) );
  f = T1_3_2.locate(p0,lt,li); assert( lt == Cls::OUTSIDE );
  li = f->index(T1_3_2.infinite_vertex());
  assert( _test_is_to_the_left(T1_3_2,p0,f,li) );
  f = T1_3_2.locate(p7,lt,li); assert( lt == Cls::OUTSIDE );
  li = f->index(T1_3_2.infinite_vertex());
  assert( _test_is_to_the_left(T1_3_2,p7,f,li) );
  f = T1_3_2.locate(p5,lt,li); assert( lt == Cls::OUTSIDE );
  li = f->index(T1_3_2.infinite_vertex());
  assert( _test_is_to_the_left(T1_3_2,p5,f,li) );
  f = T1_3_2.locate(p4,lt,li); assert( lt == Cls::OUTSIDE );
  li = f->index(T1_3_2.infinite_vertex());
  assert( _test_is_to_the_left(T1_3_2,p4,f,li) );
  f = T1_3_2.locate(p6,lt,li); assert( lt == Cls::OUTSIDE );
  li = f->index(T1_3_2.infinite_vertex());
  assert( _test_is_to_the_left(T1_3_2,p6,f,li) );

  // Check point location in 2-dimensional triangulations
  cout << "    point locations 2-dim" << endl;
  f = T2_1.locate(p0,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_1.geom_traits().compare(f->vertex(li)->point(), p0) );
  f = T2_1.locate(p1,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_1.geom_traits().compare(f->vertex(li)->point(), p1) );
  f = T2_1.locate(p2,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_1.geom_traits().compare(f->vertex(li)->point(), p2) );
  f = T2_1.locate(p3,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_1.geom_traits().compare(f->vertex(li)->point(), p3) );
  f = T2_1.locate(p4,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_1.geom_traits().compare(f->vertex(li)->point(), p4) );
  f = T2_1.locate(p5,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_1.geom_traits().compare(f->vertex(li)->point(), p5) );
  f = T2_1.locate(p6,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_1.geom_traits().compare(f->vertex(li)->point(), p6) );
  f = T2_1.locate(p7,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_1.geom_traits().compare(f->vertex(li)->point(), p7) );
  f = T2_1.locate(p8,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_1.geom_traits().compare(f->vertex(li)->point(), p8) );
  f = T2_1.locate(p9,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_1.geom_traits().compare(f->vertex(li)->point(), p9) );
  f = T2_1.locate(p10,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_1.geom_traits().compare(f->vertex(li)->point(), p10) );

  f = T2_1.locate(p11,lt,li); assert( lt == Cls::EDGE);
  assert( (T2_1.geom_traits().compare(f->vertex(f->ccw(li))->point(), p1)
        && T2_1.geom_traits().compare(f->vertex(f->cw(li))->point(), p0))
       || (T2_1.geom_traits().compare(f->vertex(f->ccw(li))->point(), p0)
        && T2_1.geom_traits().compare(f->vertex(f->cw(li))->point(), p1)));
  f = T2_1.locate(p12,lt,li); assert( lt == Cls::FACE );
  assert( T2_1.oriented_side(f,p12) == CGAL::ON_POSITIVE_SIDE );
  f = T2_1.locate(p13,lt,li,f); assert( lt == Cls::OUTSIDE );
  li = f->index(T2_1.infinite_vertex());
  assert( _test_is_to_the_left(T2_1,p13,f,li) );
  f = T2_1.locate(p14,lt,li); assert( lt == Cls::OUTSIDE );
  li = f->index(T2_1.infinite_vertex());
  assert( _test_is_to_the_left(T2_1,p14,f,li) );
  f = T2_1.locate(p15,lt,li); assert( lt == Cls::OUTSIDE );
  li = f->index(T2_1.infinite_vertex());
  assert( _test_is_to_the_left(T2_1,p15,f,li) );


  /*************************/
  /******* Iterators *******/
  cout << "    iterators" << endl;
  // _test_iterators(T0_0);
  // _test_iterators(T0_1);
  _test_iterators(T1_2);
  _test_iterators(T1_3_0);
  _test_iterators(T1_3_1);
  _test_iterators(T1_5);
  _test_iterators(T1_6);
  _test_iterators(T2_1);
  _test_iterators(T2_3);
  _test_iterators(T2_5);
  _test_iterators(T2_6);
  _test_iterators(T2_7);

  /***************************/
  /******* Circulators *******/
  cout << "    circulators" << endl;
  // _test_circulators(T0_0);
  // _test_circulators(T0_1);
  _test_circulators(T1_2);
  _test_circulators(T1_3_0);
  _test_circulators(T1_3_1);
  _test_circulators(T1_5);
  _test_circulators(T1_6);
  _test_circulators(T2_1);
  _test_circulators(T2_3);
  _test_circulators (T2_5);
  _test_circulators(T2_6);
  _test_circulators(T2_7);


  // Line_face_circulator
  cout << "    line face circulator  " << endl;
  typedef typename Cls::Line_face_circulator LFC;
  // here == operator needed for Point!
  // testing with the grid triangulation
  LFC fc= T2_7.line_walk(p1,p10);
  assert(fc.ptr()!=NULL);
  assert(!fc.is_empty());
  LFC fc2=fc;
  assert(fc==fc2);
  fc++;
  fc--;
  ++fc;
  --fc;
  Point pp(0,0.5,1);
  f = T2_7.locate(pp,lt,li);
  assert(lt==Cls::VERTEX);
  fc= T2_7.line_walk(pp,p10,f);
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
  fc= T2_8.line_walk(Point(0.5,0.4),Point(5,5));
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
  f = TT.locate(Point(0,0));
  fc = TT.line_walk(Point(0,0),Point(1,1));
  fc2 = TT.line_walk(Point(0,0),Point(1,1),f);
  assert(fc==fc2);
  n=0;
  do {fc2++ ; n = n+1;} while (fc2 != fc);
  assert(n==3);

  /*****************************/
  /******** Miscellaneaous *****/
  cout << "    misc." << endl;
  assert( T0_0.ccw(0) == 1 );
  assert( T0_0.ccw(1) == 2 );
  assert( T0_0.ccw(2) == 0 );
  assert( T0_0.cw(0) == 2 );
  assert( T0_0.cw(1) == 0 );
  assert( T0_0.cw(2) == 1 );

  // the assert() are to avoid compiler warnings about unused variables
  f = T2_1.locate(p12,lt,li); // from section locate above
  Triangle t = T2_1.triangle(f); assert( &t == &t );
  Segment  s = T2_1.segment(f,li); assert( &s == &s );
  s = T2_1.segment(Edge(f,li)); assert( &s == &s );
  s = T2_1.segment(v2_1_6->incident_edges()); assert( &s == &s );
  s = T2_1.segment(T2_1.edges_begin()); assert( &s == &s );


  /********************/
  /******** I/O *******/
  cout << "    output to a file" << endl;
  ofstream of0_0("T00.triangulation", ios::out);
  CGAL::set_ascii_mode(of0_0); 
  of0_0 << T0_1; of0_0.close();
  ofstream of0_1("T01.triangulation");
  CGAL::set_ascii_mode(of0_1); 
  of0_1 << T0_1; of0_1.close();
  ofstream of1_2("T12.triangulation");
  CGAL::set_ascii_mode(of1_2); 
  of1_2 << T1_2; of1_2.close();
  ofstream of1_5("T15.triangulation");
  CGAL::set_ascii_mode(of1_5); 
  of1_5 << T1_5; of1_5.close();
  ofstream of1_6("T16.triangulation");
  CGAL::set_ascii_mode(of1_6); 
  of1_6 << T1_6; of1_6.close();
  ofstream of2_1("T21.triangulation");
  CGAL::set_ascii_mode(of2_1); 
  of2_1 << T2_1; of2_1.close();
  ofstream of2_3("T23.triangulation");
  CGAL::set_ascii_mode(of2_3); 
  of2_3 << T2_3; of2_3.close();
  ofstream of2_5("T25.triangulation");
  CGAL::set_ascii_mode(of2_5); 
  of2_5 << T2_5; of2_5.close();
  ofstream of2_6("T26.triangulation");
  CGAL::set_ascii_mode(of2_6); 
  of2_6 << T2_6; of2_6.close();

  cout << "    input from a file" << endl;
  ifstream if0_0("T00.triangulation"); CGAL::set_ascii_mode(if0_0);
  Cls T0_0_copy;   if0_0 >> T0_0_copy;
  // assert( T0_0_copy.number_of_vertices() == T0_0.number_of_vertices() );
  ifstream if0_1("T01.triangulation"); CGAL::set_ascii_mode(if0_1);
  Cls T0_1_copy; if0_1 >> T0_1_copy;
  // assert( T0_1_copy.number_of_vertices() == T0_1.number_of_vertices() );
  ifstream if1_2("T12.triangulation"); CGAL::set_ascii_mode(if1_2); 
  Cls T1_2_copy; if1_2 >> T1_2_copy;
  // assert( T1_2_copy.number_of_vertices() == T1_2.number_of_vertices() );
  ifstream if1_5("T15.triangulation"); CGAL::set_ascii_mode(if1_5); 
  Cls T1_5_copy; if1_5 >> T1_5_copy;
  // assert( T1_5_copy.number_of_vertices() == T1_5.number_of_vertices() );
  ifstream if1_6("T16.triangulation"); CGAL::set_ascii_mode(if1_6);
  Cls T1_6_copy; if1_6 >> T1_6_copy;
  // assert( T1_6_copy.number_of_vertices() == T1_6.number_of_vertices() );
  ifstream if2_1("T21.triangulation"); CGAL::set_ascii_mode(if2_1);
  Cls T2_1_copy; if2_1 >> T2_1_copy;
  // assert( T2_1_copy.number_of_vertices() == T2_1.number_of_vertices() );
  ifstream if2_3("T23.triangulation"); CGAL::set_ascii_mode(if2_3);
  Cls T2_3_copy; if2_3 >> T2_3_copy;
  // assert( T2_3_copy.number_of_vertices() == T2_3.number_of_vertices() );
  ifstream if2_5("T25.triangulation"); CGAL::set_ascii_mode(if2_5); 
  Cls T2_5_copy; if2_5 >> T2_5_copy;
  // assert( T2_5_copy.number_of_vertices() == T2_5.number_of_vertices() );
  ifstream if2_6("T26.triangulation"); CGAL::set_ascii_mode(if2_6);
  Cls T2_6_copy; if2_6 >> T2_6_copy;
  // assert( T2_6_copy.number_of_vertices() == T2_6.number_of_vertices() );

  /**********************/
  /***** REMOVALS *******/ 
  cout << "    removals" << endl;

//   // test remove_first()
//   T0_1.remove_first(T0_1.finite_vertex());
//   assert( T0_1.number_of_vertices() == 0 );

  // test remove_second()
  T1_6.remove_second(T1_6.finite_vertex());
  assert( T1_6.number_of_vertices() == 1 );

  // remove from 1-dimensional triangulations
  T1_2.remove(v1_2_1);
  T1_2.remove(v1_2_2);
  assert( T1_2.number_of_vertices() == 0 );

  T1_5.remove(v1_5_1);
  T1_5.remove(v1_5_2);
  T1_5.remove(v1_5_3);
  T1_5.remove(v1_5_8);
  T1_5.remove(v1_5_9);
  assert( T1_5.number_of_vertices() == 0 );

  // remove from 2-dimensional triangulations
  T2_1.remove(v2_1_0);
  T2_1.remove(v2_1_2);
  T2_1.remove(v2_1_1);
  T2_1.remove(v2_1_6);
  T2_1.remove(v2_1_5);
  T2_1.remove(v2_1_4);
  T2_1.remove(v2_1_3);
  T2_1.remove(v2_1_9);
  T2_1.remove(v2_1_8);
  T2_1.remove(v2_1_7);
  T2_1.remove(v2_1_10);
  assert( T2_1.number_of_vertices() == 0 );

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
  
  for (i=T2_7.number_of_vertices(); i>0; i--)
    T2_7.remove(T2_7.finite_vertex());
  assert( T2_7.number_of_vertices() == 0 );

  // test destructors and return
  cout << "    test destructors and return" << endl;
}
