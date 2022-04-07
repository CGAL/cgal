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
// file          : include/CGAL/_test_cls_regular_triangulation_2.h
// revision      :
// revision_date :

// author(s)     : Francois Rebufat, Mariette Yvinec

// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#ifndef TEST_CLS_REGULAR_TRIANGULATION_C
#define TEST_CLS_REGULAR_TRIANGULATION_C

#include <iostream>
#include <fstream>


#include <CGAL/_test_fct_is_infinite.h>
#include <CGAL/_test_triangulation_iterators.h>
#include <CGAL/_test_triangulation_circulators.h>


template <class Del>
void
_test_regular_duality( const Del &T );


template < class Triangulation, class Point, class Face_handle >
bool
_test_is_to_the_left( const Triangulation &T,
                      const Point &p,
                      const Face_handle &f,
                      const int li)
{
  typename Triangulation::Geom_traits::Construct_weighted_point_2 p2wp =
      T.geom_traits().construct_weighted_point_2_object();

  return( T.orientation(f->vertex(f->ccw(li))->point(),
                        f->vertex(f->cw(li))->point(),
                        p2wp(p)) == CGAL::LEFT_TURN );
}

template <class Triangulation>
void
_test_cls_regular_triangulation_2( const Triangulation & )
{
  typedef Triangulation                      Cls;

  // We assume the traits class has been tested already
  // actually, any traits is good if it has been tested
  typedef typename Cls::Geom_traits          Gt;

  typedef typename Cls::Point                Point;
  typedef typename Cls::Bare_point           Bare_point;
  typedef typename Cls::Weighted_point       Weighted_point;
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
  typedef typename Cls::size_type            size_type;

  CGAL_USE_TYPE(Gt);
  CGAL_USE_TYPE(Point);
  CGAL_USE_TYPE(Vertex);
  CGAL_USE_TYPE(Face);
  CGAL_USE_TYPE(Finite_vertices_iterator);
  CGAL_USE_TYPE(Finite_faces_iterator);
  CGAL_USE_TYPE(Finite_edges_iterator);
  CGAL_USE_TYPE(Vertex_circulator);
  CGAL_USE_TYPE(Face_circulator);
  CGAL_USE_TYPE(Edge_circulator);
  CGAL_USE_TYPE(Line_face_circulator);
  // activate verbose will print the number of vertices and hidden
  // vertices
  // in each triangulation tested by is_valid(verbose);
  //bool verbose = true ;
  bool verbose = false;

  // Build a few objects
  // p1,p3,p2,p9,p8 aligned in this order
  // p0,p5,p7 also

  Bare_point p0(5,6,1);
  Bare_point p1(1,9,1);
  Bare_point p2(6,14,1);
  Bare_point p3(4,12,1);
  Bare_point p4(3,29,1);
  Bare_point p5(6,7,1);
  Bare_point p6(6,39,1);
  Bare_point p7(8,9,1);
  Bare_point p8(10,18,1);
  Bare_point p9(75625,155625,10000); // intersection of p2,p8 and p6,p7
  Bare_point p10(10,50,2);
  Bare_point p11(6,15,2); // midpoint p1,p0
  Bare_point p12(6,16,2); // slightly above, in face
  Bare_point p13(10,11,1);
  Bare_point p14(10,40,1);
  Bare_point p15(60,-10,1);

  int px=1, py=1;
  int qx=-1, qy=2;


  Weighted_point wp0(p0,1);
  Weighted_point wp1(p1,20);
  Weighted_point wp2(p2,20);
  Weighted_point wp3(p3,1);
  //Weighted_point wp3(p3,1);
  Weighted_point wp4(p4,1);
  Weighted_point wp5(p5,8);
  Weighted_point wp6(p6,1);
  Weighted_point wp7(p7,1);
  Weighted_point wp8(p8,20);
  Weighted_point wp9(p9,1); // intersection of p2,p8 and p6,p7
  Weighted_point wp10(p10,2);
  Weighted_point wp11(p11,2); // midpoint p1,p0
  Weighted_point wp12(p12,2); // slightly above, in face
  Weighted_point wp13(p13,1);
  Weighted_point wp14(p14,1);
  Weighted_point wp15(p15,1);
  Weighted_point wp16(p2,1);
  Weighted_point wp17(p3,20);
  Weighted_point wp19(p9,0.5);
  Weighted_point wp29(p9,22);
  Weighted_point wp22(p12,300);

  {
    Weighted_point p15_bis(p15.x(), p15.y());
    assert(p15_bis == p15);
  }

  Cls T;

  std::cerr << wp1 << " " << wp1.x()  << std::endl;
  std::cerr << wp2 << std::endl;
  std::cerr << wp3 << std::endl;
  assert(T.power_test(wp1,wp2,wp3) == CGAL::ON_NEGATIVE_SIDE);
  assert(T.power_test(wp1,wp8,wp2) == CGAL::ON_POSITIVE_SIDE);
  assert(T.power_test(wp2,wp8,wp9) == CGAL::ON_NEGATIVE_SIDE);
  assert(T.power_test(wp1,wp9,wp3) == CGAL::ON_POSITIVE_SIDE);

  std::list<Weighted_point> lw; lw.push_back(wp0);
  lw.push_back(wp1); lw.push_back(wp2); lw.push_back(wp3);


  std::vector<Weighted_point> vw; vw.push_back(wp0);
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
  assert( T0_0.is_valid(verbose) );

  Cls T0_1;
  Vertex_handle v0_1_0 = T0_1.insert(wp0); assert( v0_1_0 != nullptr );
  assert( T0_1.dimension() == 0 );
  assert( T0_1.number_of_vertices() == 1 );
  assert( T0_1.is_valid(verbose) );

  Cls T0_2;
  T0_2.insert_first(wp0);
  assert( T0_2.dimension() == 0 );
  assert( T0_2.number_of_vertices() == 1 );
  assert( T0_2.is_valid(verbose) );


  /******** 1-dimensional triangulations ******/
  // T1_n denotes a 1-dimensional triangulation with n vertices
  // when there are several, we use T1_n_p
  std::cout << "    insertions 1-dim" << std::endl;

  Cls T1_5;
  Vertex_handle v1_5_1 = T1_5.insert(wp1);
  T1_5.is_valid(verbose);
  Vertex_handle v1_5_2 = T1_5.insert(wp2);
  T1_5.is_valid(verbose);
  Vertex_handle v1_5_3 = T1_5.insert(wp3);  //hidden vertex
  T1_5.is_valid(verbose);
  Vertex_handle v1_5_9 = T1_5.insert(wp9);
  T1_5.is_valid(verbose);
  Vertex_handle v1_5_8 = T1_5.insert(wp8); //hide wp9
  T1_5.is_valid(verbose);

  assert( T1_5.dimension() == 1 );
  assert( T1_5.number_of_vertices() == 3);
  assert( T1_5.is_valid(verbose) );
  Vertex_handle v1_5_16 =  T1_5.insert(wp16); T1_5.is_valid(verbose);
  Vertex_handle v1_5_17 =  T1_5.insert(wp17); T1_5.is_valid(verbose);

  // test insert_second()
  Cls T1_6 = T0_2;
  T1_6.insert_second( wp3);
  assert( T1_6.dimension() == 1 );
  assert( T1_6.number_of_vertices() == 2 );
  assert( T1_6.is_valid(verbose) );


  Cls T2_3;
  Vertex_handle v2_3_1 = T2_3.insert(wp1);
  T2_3.is_valid(verbose);
  Vertex_handle v2_3_3 = T2_3.insert(wp3);
  T2_3.is_valid(verbose);
  Vertex_handle v2_3_9 = T2_3.insert(wp9);
  T2_3.is_valid(verbose);
  Vertex_handle v2_3_2 = T2_3.insert(wp8);
  T2_3.is_valid(verbose);
  Vertex_handle v2_3_8 = T2_3.insert(wp2);
  T2_3.is_valid(verbose);
  assert( T2_3.dimension() == 1 );
  T2_3.is_valid(verbose);

    /******** 2-dimensional triangulations ******/
  std::cout << "    insertions 2-dim" << std::endl;
  Locate_type lt;
  Face_handle loc;
  int li;

  Vertex_handle v2_3_4 = T2_3.insert(wp4);
  assert( T2_3.dimension() == 2 );
  T2_3.is_valid(verbose);
  Vertex_handle v2_3_6 = T2_3.insert(wp6, T2_3.finite_faces_begin());
  T2_3.is_valid(verbose);
  Vertex_handle v2_3_0 = T2_3.insert(wp0, ++T2_3.finite_faces_begin());
  T2_3.is_valid(verbose);
  Vertex_handle v2_3_5 = T2_3.insert(wp5);
  T2_3.is_valid(verbose);
  Vertex_handle v2_3_7 = T2_3.insert(wp7);
  T2_3.is_valid(verbose);
  loc = T2_3.locate(wp10,lt,li);
  Vertex_handle v2_3_10 = T2_3.insert(wp10, lt, loc,li);
  assert( T2_3.dimension() == 2 );
  assert( T2_3.is_valid(verbose) );


  // test list iterator insert
  Cls T2_5;
  assert( T2_5.insert(lw.begin(), lw.end()) == 3);
  assert( T2_5.dimension() == 2 );
  assert( T2_5.number_of_vertices() == 3);
  assert( T2_5.is_valid(verbose) );

  // test list iterator insert
  Cls T2_6;
  assert( T2_6.insert(vw.begin(), vw.end()) == 3);
  assert( T2_6.dimension() == 2 );
  assert( T2_6.number_of_vertices() == 3);
  assert( T2_6.is_valid(verbose) );

  // test grid insert and make sure push_back exists
  Cls T2_7;
  int m, p;
  for (m=0; m<3; m++)
    for (p=0; p<3; p++)
      T2_7.push_back( Weighted_point(Bare_point(m*px+p*qx, m*py+p*qy), 1) );
  assert( T2_7.number_of_vertices() == static_cast<size_type>(m*p) );
  assert( T2_7.is_valid(verbose) );


 std::cout << "    constructors (2)" << std::endl;

  // test copy_constructor with non-empty 0-triangulation
  Cls T0_1_1( T0_1 );
  assert( T0_1_1.dimension() == 0 );
  assert( T0_1_1.number_of_vertices() == 1 );
  assert( T0_1_1.is_valid(verbose) );

  // test assignement
  Cls T0_1_2;
  T0_1_2 = T0_1;
  assert( T0_1_2.dimension() == 0 );
  assert( T0_1_2.number_of_vertices() == 1 );
  assert( T0_1_2.is_valid(verbose) );

  // test copy_constructor with non-empty 1-triangulation
  Cls T1_5_1( T1_5 );
  assert( T1_5_1.dimension() == T1_5.dimension() );
  assert( T1_5_1.number_of_vertices() == T1_5.number_of_vertices());
  assert( T1_5_1.number_of_hidden_vertices()==
          T1_5.number_of_hidden_vertices() );
  assert( T1_5_1.is_valid(verbose) );

  // Test assignment operator
  Cls T1_5_2 = T1_5;
  assert( T1_5_2.dimension() == T1_5.dimension());
  assert( T1_5_2.number_of_vertices() == T1_5.number_of_vertices());
  assert( T1_5_2.number_of_hidden_vertices()==
          T1_5.number_of_hidden_vertices() );
  assert( T1_5_2.is_valid(verbose) );

   // test copy_constructor with non-empty 2-triangulation
  Cls T2_3_1( T2_3 );
  assert( T2_3_1.dimension() == T2_3.dimension());
  assert( T2_3_1.number_of_vertices() == T2_3.number_of_vertices());
  assert( T2_3_1.number_of_hidden_vertices()==
          T2_3.number_of_hidden_vertices() );
  assert( T2_3_1.is_valid(verbose) );

  // test assignment operator
  Cls T2_3_4 = T2_3;
  assert( T2_3_4.dimension() == T2_3.dimension() );
  assert( T2_3_4.number_of_vertices() == T2_3.number_of_vertices());
  assert( T2_3_4.number_of_hidden_vertices()==
          T2_3.number_of_hidden_vertices() );
  assert( T2_3_4.is_valid(verbose) );


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
  T1_3_2.insert(wp1);
  T1_3_2.insert(wp2);
  T1_3_2.insert(wp9);
  T1_3_2.is_valid(verbose);
  loc = T1_3_2.locate(wp1,lt,li); assert( lt == Cls::VERTEX );
  assert( T1_3_2.xy_equal(loc->vertex(li)->point(), wp1) );
  loc = T1_3_2.locate(wp2,lt,li); assert( lt == Cls::VERTEX );
  assert( T1_3_2.xy_equal(loc->vertex(li)->point(), wp2) );
  loc = T1_3_2.locate(wp9,lt,li); assert( lt == Cls::VERTEX );
  assert( T1_3_2.xy_equal(loc->vertex(li)->point(), wp9) );
  loc = T1_3_2.locate(wp3,lt,li); assert( lt == Cls::EDGE );
  assert( (T1_3_2.xy_equal(loc->vertex(loc->ccw(li))->point(), wp1)
        && T1_3_2.xy_equal(loc->vertex(loc->cw(li))->point(), wp2))
       || (T1_3_2.xy_equal(loc->vertex(loc->ccw(li))->point(), wp2)
        && T1_3_2.xy_equal(loc->vertex(loc->cw(li))->point(), wp1)));
  loc = T1_3_2.locate(wp8,lt,li); assert( lt == Cls::OUTSIDE_CONVEX_HULL );
  loc = T1_3_2.locate(wp7,lt,li); assert( lt == Cls::OUTSIDE_AFFINE_HULL );
  loc = T1_3_2.locate(wp5,lt,li); assert( lt == Cls::OUTSIDE_AFFINE_HULL );
  loc = T1_3_2.locate(wp4,lt,li); assert( lt == Cls::OUTSIDE_AFFINE_HULL );
  loc = T1_3_2.locate(wp6,lt,li); assert( lt == Cls::OUTSIDE_AFFINE_HULL);


  // Check point location in 2-dimensional triangulations
  std::cout << "    point locations 2-dim" << std::endl;
  loc = T2_3.locate(wp0,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), wp0) );
  loc = T2_3.locate(wp1,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), wp1) );
  loc = T2_3.locate(wp2,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), wp2) );
  loc = T2_3.locate(wp4,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), wp4) );
  loc = T2_3.locate(wp5,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), wp5) );
  loc = T2_3.locate(wp6,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), wp6) );
  loc = T2_3.locate(wp7,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), wp7) );
  loc = T2_3.locate(wp8,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), wp8) );
  loc = T2_3.locate(wp10,lt,li); assert( lt == Cls::VERTEX );
  assert( T2_3.xy_equal(loc->vertex(li)->point(), wp10) );


  loc = T2_3.locate(wp3,lt,li); assert( lt == Cls::EDGE );
  loc = T2_3.locate(wp9,lt,li); assert( lt == Cls::EDGE );
  loc = T2_3.locate(wp11,lt,li); assert( lt == Cls::EDGE);
  assert( (T2_3.xy_equal(loc->vertex(loc->ccw(li))->point(), wp1)
        && T2_3.xy_equal(loc->vertex(loc->cw(li))->point(), wp0))
       || (T2_3.xy_equal(loc->vertex(loc->ccw(li))->point(), wp0)
        && T2_3.xy_equal(loc->vertex(loc->cw(li))->point(), wp1)));
  loc = T2_3.locate(wp12,lt,li); assert( lt == Cls::FACE );
  assert( T2_3.oriented_side(loc,wp12) == CGAL::ON_POSITIVE_SIDE );
  loc = T2_3.locate(wp13,lt,li,loc); assert( lt == Cls::OUTSIDE_CONVEX_HULL );
  li = loc->index(T2_3.infinite_vertex());
  assert( _test_is_to_the_left(T2_3,p13,loc,li) );
  loc = T2_3.locate(wp14,lt,li); assert( lt == Cls::OUTSIDE_CONVEX_HULL );
  li = loc->index(T2_3.infinite_vertex());
  assert( _test_is_to_the_left(T2_3,p14,loc,li) );
  loc = T2_3.locate(wp15,lt,li); assert( lt == Cls::OUTSIDE_CONVEX_HULL );
  li = loc->index(T2_3.infinite_vertex());
  assert( _test_is_to_the_left(T2_3,p15,loc,li) );


  /*********************************************/
  /****** FURTHER insertions test *******/
  std::cout << "    further insertions 2-dim" << std::endl;

  // insertion of hidden points - use copy T2_3_1 of  T2_3
  T2_3_1.insert(wp16);     //hidden on vertex
  T2_3_1.is_valid(verbose);
  T2_3_1.insert(wp19);    // hidden on edge
  T2_3_1.is_valid(verbose);

  loc = T2_3_1.locate(wp12,lt,li); assert( lt == Cls::FACE );
  assert( T2_3_1.oriented_side(loc,wp12) == CGAL::ON_POSITIVE_SIDE );
  assert( T2_3_1.power_test(loc,wp12) == CGAL::ON_NEGATIVE_SIDE);
  T2_3_1.insert(wp12); //hidden in face
  T2_3_1.is_valid(verbose);

  // insertion of hiding vertices (only wp22 is a hiding vertex)
  Vertex_handle V2_31_17 = T2_3_1.insert(wp17);
  T2_3_1.is_valid(verbose);
  Vertex_handle V2_31_29 = T2_3_1.insert(wp29);
  T2_3_1.is_valid(verbose);
  Vertex_handle V2_31_22  = T2_3_1.insert(wp22);
  T2_3_1.is_valid(verbose);


  //massive insertion
//  Cls tr;
//   int nn=10;
//    int sign = 1;
//    for(int i = 0; i < nn; i++) {
//      for(int j = 0; j< nn; j++) {
//        sign= -sign;
//        Bare_point p( i,j);
//        Weighted_point wp(p,500 + sign*500);
//        tr.insert(wp);
//        tr.is_valid();
//      }
//    }
//    tr.is_valid();

//  tr.clear();
//    Cls tr;
//    Weighted_point wp;
//    std::ifstream input("data"); CGAL::IO::set_ascii_mode(input);
//    int in = 0;
//    while(input){
//      in = in+1;
//      input >> wp;
//      std::cerr << in << " wpoint " << wp.point() << " " <<
//        wp.weight() <<std::endl;
//      tr.insert(wp);
//      tr.is_valid(true);
//    }

//    // remove all
//    std::cout <<  "  removal of all points" << std::endl;
//    while( tr.number_of_vertices() >0) {
//     tr.is_valid(true);
//     tr.remove(tr.finite_vertices_begin());
//   }

/*************************/
  /******* Iterators *******/
  std::cout << "    iterators" << std::endl;
  // In case of regular triangulation number_of_vertices() !=
  // of what the iterators can count
  // and this makes tests fail
  // _test_iterators(T0_0);
  // _test_iterators(T0_1);
   //_test_iterators(T1_5);
  // _test_iterators(T2_3);


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
  LFC fc= T2_3.line_walk(wp1,wp10);
  assert(fc != nullptr);
  assert(!fc.is_empty());
  LFC fc2=fc;
  assert(fc==fc2);
  fc++;
  fc--;
  ++fc;
  --fc;
  Bare_point pp(0,1,2);
  Weighted_point wpp(pp);
  loc = T2_3.locate(wpp,lt,li);
  fc= T2_3.line_walk(wpp,wp10,loc);
  fc2=fc;
  assert(fc==fc2);
  fc++;
  fc--;
  ++fc;
  --fc;
  // testing with dummy triangulations
  Cls T2_8;
  T2_8.insert(Weighted_point(Bare_point(0,0,1)));
  T2_8.insert(Weighted_point(Bare_point(1,0,1)));
  T2_8.insert(Weighted_point(Bare_point(0,1,1)));
  T2_8.insert(Weighted_point(Bare_point(1,1,1)));
  int n=0;
  do {fc2++ ; n = n+1;} while (fc2 != fc);
  assert(T2_8.number_of_vertices()>=2);
  assert(T2_8.is_valid(verbose));
  fc= T2_8.line_walk(Weighted_point(Bare_point(5,4,10)),Weighted_point(5,5));
  fc2=fc;
  n=0;
  assert(fc==fc2);
  do {fc2++ ; n = n+1;} while (fc2 != fc);
  assert(n==4);
  // the two point are vertices of the triangulation.
  Cls TT;
  TT.insert(Weighted_point(Bare_point(0,0)));
  TT.insert(Weighted_point(Bare_point(1,0)));
  TT.insert(Weighted_point(Bare_point(1,1)));
  TT.insert(Weighted_point(Bare_point(0,1)));
  assert(TT.dimension()==2);
  assert(TT.is_valid(verbose));
  assert(TT.number_of_vertices()==4);
  loc = TT.locate(Weighted_point(0,0));
  fc = TT.line_walk(Weighted_point(0,0),Weighted_point(1,1));
  fc2 = TT.line_walk(Weighted_point(0,0),Weighted_point(1,1),loc);
  if (fc != fc2)
    {
      TT.show_all();
      TT.show_face(loc);
      TT.show_face(fc);
      TT.show_face(fc2);
    }
  assert(fc==fc2);
  n=0;
  do {fc2++ ; n = n+1;} while (fc2 != fc);
  assert(n==4);

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
  loc = T2_3.locate(wp12,lt,li); // from section locate above
  Triangle t = T2_3.triangle(loc); assert( &t == &t );
  Segment  s = T2_3.segment(loc,0); assert( &s == &s );
  s = T2_3.segment(Edge(loc,1)); assert( &s == &s );
  s = T2_3.segment(T2_3.incident_edges(v2_3_6)); assert( &s == &s );
  s = T2_3.segment(T2_3.finite_edges_begin()); assert( &s == &s );


  std::cout << "    hidden vertices" << std::endl;
  {
      Cls rt;

      rt.insert (wp0);
      rt.insert (wp1);
      rt.insert (wp2);
      rt.insert (wp3);

      for (unsigned i = 4; i > 0; --i) {
          assert (rt.number_of_vertices() + rt.number_of_hidden_vertices() == i);
          assert (rt.number_of_vertices() > 0);
          rt.remove (rt.finite_vertices_begin ());
      }
      assert (rt.number_of_vertices() == 0);
      assert (rt.number_of_hidden_vertices() == 0);
  }

  /***************************/
  /***** get_conflicts ******/
  std:: cout << "    get conflicts " << std::endl;
  std::list<Face_handle> conflicts;
  std::list<Edge>  hole_bd;
  std::list<Vertex_handle>  hidden_vertices;
  std::back_insert_iterator<std::list<Face_handle> > c_inserter(conflicts);
  std::back_insert_iterator<std::list<Edge> > be_inserter(hole_bd);
  std::back_insert_iterator<std::list<Vertex_handle> >
    v_inserter(hidden_vertices);
  CGAL::Triple<std::back_insert_iterator<std::list<Face_handle> >,
     std::back_insert_iterator<std::list<Edge> >,
     std::back_insert_iterator<std::list<Vertex_handle> > >
    tit(c_inserter,be_inserter,v_inserter);
  std::pair<std::back_insert_iterator<std::list<Face_handle> >,
    std::back_insert_iterator<std::list<Edge> > >
    cbe_pit(c_inserter,be_inserter);
  std::pair<std::back_insert_iterator<std::list<Face_handle> >,
    std::back_insert_iterator<std::list<Vertex_handle> > >
    cv_pit(c_inserter,v_inserter);
  std::pair<std::back_insert_iterator<std::list<Edge> >,
    std::back_insert_iterator<std::list<Vertex_handle> > >
    bev_pit(be_inserter,v_inserter);
  //point that exists:
  tit = T2_3.get_conflicts_and_boundary_and_hidden_vertices
    (wp5,std::back_inserter(conflicts),
     std::back_inserter(hole_bd),
     std::back_inserter(hidden_vertices));
  assert(conflicts.empty() && hole_bd.empty() && hidden_vertices.empty());

  //on a point that exists with lower weight:
  tit = T2_3.get_conflicts_and_boundary_and_hidden_vertices
    (Weighted_point(p5,7),
     std::back_inserter(conflicts),
     std::back_inserter(hole_bd),
     std::back_inserter(hidden_vertices));
  assert(conflicts.empty() && hole_bd.empty() &&
         hidden_vertices.empty());

  //on a point that exists with higher weight:
  tit = T2_3.get_conflicts_and_boundary_and_hidden_vertices
    (Weighted_point(p5,9),
     std::back_inserter(conflicts),
     std::back_inserter(hole_bd),
     std::back_inserter(hidden_vertices));
  assert(!hidden_vertices.empty());
  assert(2*hidden_vertices.size() +  hole_bd.size() - conflicts.size()
         == 2);
  conflicts.clear();
  hole_bd.clear();
  hidden_vertices.clear();

  //on a point that exists with big weight:
  tit = T2_3.get_conflicts_and_boundary_and_hidden_vertices
    (Weighted_point(p5,150),
     std::back_inserter(conflicts),
     std::back_inserter(hole_bd),
     std::back_inserter(hidden_vertices));
  assert(2*hidden_vertices.size() +  hole_bd.size() - conflicts.size()
         == 2);
  conflicts.clear();
  hole_bd.clear();
  hidden_vertices.clear();

  //hidden vertices:
  tit = T2_3.get_conflicts_and_boundary_and_hidden_vertices
    (wp16,
     std::back_inserter(conflicts),
     std::back_inserter(hole_bd),
     std::back_inserter(hidden_vertices));
  tit = T2_3.get_conflicts_and_boundary_and_hidden_vertices
    (wp19,
     std::back_inserter(conflicts),
     std::back_inserter(hole_bd),
     std::back_inserter(hidden_vertices));
  loc = T2_3.locate(wp12,lt,li); assert( lt == Cls::FACE );
  tit = T2_3.get_conflicts_and_boundary_and_hidden_vertices
    (wp12,
     std::back_inserter(conflicts),
     std::back_inserter(hole_bd),
     std::back_inserter(hidden_vertices), loc);
  c_inserter = tit.first;
  be_inserter = tit.second;
  v_inserter = tit.third;
  assert(conflicts.empty() && hole_bd.empty() && hidden_vertices.empty());

  //non-hiding vertex:
  v_inserter = T2_3.get_hidden_vertices
    (wp17,std::back_inserter(hidden_vertices));
  cbe_pit = T2_3.get_conflicts_and_boundary
    (wp17,
     std::back_inserter(conflicts),
     std::back_inserter(hole_bd));
  assert(hidden_vertices.empty());
  assert(hole_bd.size() == conflicts.size() + 2);
  c_inserter = cbe_pit.first;
  be_inserter =cbe_pit.second;
  conflicts.clear();
  hole_bd.clear();

  //hiding vertex:
  cv_pit = T2_3.get_conflicts_and_hidden_vertices
    (wp22,
     std::back_inserter(conflicts),
     std::back_inserter(hidden_vertices));
  c_inserter = cv_pit.first;
  v_inserter = cv_pit.second;
  hidden_vertices.clear();

  bev_pit = T2_3.get_boundary_of_conflicts_and_hidden_vertices
    (wp22,
     std::back_inserter(hole_bd),
     std::back_inserter(hidden_vertices));
  be_inserter = bev_pit.first;
  v_inserter = bev_pit.second;
  assert(2*hidden_vertices.size() +  hole_bd.size() - conflicts.size()
         == 2);
  conflicts.clear();
  hole_bd.clear();
  hidden_vertices.clear();
  /********************/

  /***** Duality ******/
  std::cout << "    duality" << std::endl;
  _test_regular_duality(T1_5);
   _test_regular_duality(T2_3);

  /********************/
  /******** I/O *******/
   // INPUT-OUTPUT still to be rwertten
   // input output have not yet been overload
   // so they do not input output hidden vertices
  std::cout << "    output to a file" << std::endl;
  std::ofstream of1_5("T15.triangulation");
  CGAL::IO::set_ascii_mode(of1_5);
  of1_5 << T1_5; of1_5.close();
  std::ofstream of2_3("T23.triangulation");
  CGAL::IO::set_ascii_mode(of2_3);
  of2_3 << T2_3; of2_3.close();


//   std::cout << "    input from a file" << std::endl;

//   std::ifstream if1_5("T15.triangulation"); CGAL::IO::set_ascii_mode(if1_5);
//   Cls T1_5_copy; if1_5 >> T1_5_copy;
 //  assert( T1_5_copy.is_valid(verbose) &&
//           T1_5_copy.number_of_vertices() ==
//           T1_5.number_of_vertices() - T1_5.number_of_hidden_vertices());

//   std::ifstream if2_3("T23.triangulation"); CGAL::IO::set_ascii_mode(if2_3);
//   Cls T2_3_copy; if2_3 >> T2_3_copy;
  // assert( T2_3_copy.is_valid(verbose) &&
//           T2_3_copy.number_of_vertices() ==
//           T2_3.number_of_vertices() - T2_3.number_of_hidden_vertices());

  /**********************/
  /***** REMOVALS *******/
  std::cout << "    removals" << std::endl;

  // test remove_first()
  T0_1.remove_first(T0_1.finite_vertex());
  assert( T0_1.number_of_vertices() == 0 );

   // remove from 1-dimensional triangulations
  T1_5.is_valid(verbose);
  T1_5.remove(v1_5_16);
  T1_5.is_valid(verbose);
  T1_5.remove(v1_5_17);
  T1_5.is_valid(verbose);
  T1_5.remove(v1_5_1);
  T1_5.is_valid(verbose);
  T1_5.remove(v1_5_2);
  T1_5.is_valid(verbose);
  T1_5.remove(v1_5_3);
  T1_5.is_valid(verbose);
  T1_5.remove(v1_5_8);
  T1_5.is_valid(verbose);
  T1_5.remove(v1_5_9);
  T1_5.is_valid(verbose);
  assert( T1_5.number_of_vertices() == 0 );

  // remove from 2-dimensional triangulations
  T2_3.is_valid(verbose);
  T2_3.remove(v2_3_0);
  T2_3.is_valid(verbose);
  T2_3.remove(v2_3_1);
  T2_3.is_valid(verbose);
  T2_3.remove(v2_3_9);
  T2_3.is_valid(verbose);
  T2_3.remove(v2_3_8);
  T2_3.is_valid(verbose);
  T2_3.remove(v2_3_5);
  T2_3.is_valid(verbose);
  T2_3.remove(v2_3_3);
  T2_3.is_valid(verbose);
  T2_3.remove(v2_3_4);
  T2_3.is_valid(verbose);
  T2_3.remove(v2_3_2);
  T2_3.is_valid(verbose);
  T2_3.remove(v2_3_6);
  T2_3.is_valid(verbose);
  T2_3.remove(v2_3_7);
  T2_3.is_valid(verbose);
  T2_3.remove(v2_3_10);
  T2_3.is_valid(verbose);
  assert( T2_3.number_of_vertices() == 0 );

  // further remove
  std::cerr << "further removal " << std::endl;
  T2_3_1.is_valid(verbose);
  T2_3_1.remove(V2_31_22);
  T2_3_1.is_valid(verbose);
  T2_3_1.remove(V2_31_29);
  T2_3_1.is_valid(verbose);
  T2_3_1.remove(V2_31_17);
  T2_3_1.is_valid(verbose);

  typename Cls::size_type i;
  T2_5.clear();
  assert( T2_5.number_of_vertices() == 0 );

  for (i=T2_6.number_of_vertices() + T2_6.number_of_hidden_vertices(); i>0; i--){
    T2_6.remove(T2_6.finite_vertex());
  }
  assert( T2_6.number_of_vertices() == 0 );


  // test destructors and return
  std::cout << "    test destructors and return" << std::endl;
}


template <class Del>
void
_test_regular_duality( const Del &T )
{
  typedef typename Del::Geom_traits          Gt;
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
  Edge_circulator ec=T.incident_edges(T.finite_vertices_begin()), done(ec);
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

#endif //TEST_CLS_REGULAR_TRIANGULATION_C
