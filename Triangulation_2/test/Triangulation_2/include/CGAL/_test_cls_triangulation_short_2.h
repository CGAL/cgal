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
// file          : include/CGAL/_test_cls_triangulation_2.h
// revision      :
// revision_date :

// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)

// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <list>

#include <CGAL/_test_fct_is_infinite.h>
#include <CGAL/_test_triangulation_iterators.h>
#include <CGAL/_test_triangulation_circulators.h>
#include <CGAL/Testsuite/Triangulation_23/test_move_semantic.h>


template <class Triangul>
void
_test_cls_triangulation_short_2( const Triangul &)
{
  //typedef Triangul                     Cls;

  // We assume the traits class has been tested already
  // actually, any traits is good if it has been tested
  typedef typename Triangul::Geom_traits          Gt;

  typedef typename Triangul::Point                Point;
  typedef typename Triangul::Segment              Segment;
  typedef typename Triangul::Triangle             Triangle;

  typedef typename Triangul::Vertex               Vertex;
  typedef typename Triangul::Face                 Face;

  typedef typename Triangul::Vertex_handle        Vertex_handle;
  typedef typename Triangul::Face_handle          Face_handle;

  typedef std::pair<Face_handle,int>              Edge;

  typedef typename Triangul::Finite_vertices_iterator
                                                      Finite_vertices_iterator;
  typedef typename Triangul::Finite_faces_iterator    Finite_faces_iterator;
  typedef typename Triangul::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename Triangul::Point_iterator           Point_iterator;

  typedef typename Triangul::Vertex_circulator    Vertex_circulator;
  typedef typename Triangul::Face_circulator      Face_circulator;
  typedef typename Triangul::Edge_circulator      Edge_circulator;
  typedef typename Triangul::Line_face_circulator Line_face_circulator;

  typedef typename Triangul::Locate_type          Locate_type;

  CGAL_USE_TYPE(Gt);
  CGAL_USE_TYPE(Vertex);
  CGAL_USE_TYPE(Face);
  CGAL_USE_TYPE(Finite_vertices_iterator);
  CGAL_USE_TYPE(Finite_faces_iterator);
  CGAL_USE_TYPE(Finite_edges_iterator);
  CGAL_USE_TYPE(Vertex_circulator);
  CGAL_USE_TYPE(Face_circulator);
  CGAL_USE_TYPE(Edge_circulator);
  CGAL_USE_TYPE(Line_face_circulator);
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

  Locate_type lt;
  int li;
  Face_handle    f, loc;


  std::list<Point> l; l.push_back(p0);
  l.push_back(p1); l.push_back(p2); l.push_back(p3);

  std::vector<Point> v; v.push_back(p0);
  v.push_back(p1); v.push_back(p2); v.push_back(p3);



  /*****************************/
  /***** CONSTRUCTORS (1) ******/
  std::cout << "    constructors(1)" << std::endl;

  Triangul T1;
  assert( T1.dimension() == -1 );
  assert( T1.number_of_vertices() == 0 );

  Triangul T3(T1);
  Triangul T4 = T1;
  T3.swap(T1);


  /**************************/
  /******* INSERTIONS *******/

  // Tk denotes a k-dimensional triangulation
  // the handle returned when inserting pj into Tk_n is called vk_n_j
  // the asserts at the end of the insert() are to avoid compiler
  // warnings about unused variables
  // we use some of these variables for the remove (see the end)
  // so there is no need to put assert at the end of all of them

  /******* 0-dimensional triangulations ******/
  std::cout << "    insertions 0-dim" << std::endl;

  Triangul T0_0;
  assert( T0_0.dimension() == -1 );
  assert( T0_0.number_of_vertices() == 0 );
  assert( T0_0.is_valid() );

  Triangul T0_1;
  Vertex_handle v0_1_0 = T0_1.insert(p0); assert( v0_1_0 != nullptr );
  assert( T0_1.dimension() == 0 );
  assert( T0_1.number_of_vertices() == 1 );
  assert( T0_1.is_valid() );

  // test insert_first()
  Triangul T0_2;
  Vertex_handle v0_2_0 =   T0_2.insert_first(p0);
  assert( v0_2_0 != nullptr );
  assert( T0_2.dimension() == 0 );
  assert( T0_2.number_of_vertices() == 1 );
  assert( T0_2.is_valid() );

  /******** 1-dimensional triangulations ******/
  // T1_n denotes a 1-dimensional triangulation with n vertices
  // when there are several, we use T1_n_p
  std::cout << "    insertions 1-dim" << std::endl;

  Triangul T1_5;
  Vertex_handle v1_5_1 = T1_5.insert(p1);
  Vertex_handle v1_5_2 = T1_5.insert(p2);
  Vertex_handle v1_5_3 = T1_5.insert(p3);
  Vertex_handle v1_5_8 = T1_5.insert(p8);
  Vertex_handle v1_5_9 = T1_5.insert(p9);
  assert( T1_5.dimension() == 1 );
  assert( T1_5.number_of_vertices() == 5 );
  assert( T1_5.number_of_faces() == 0 );
  assert( T1_5.is_valid() );


  /******** 2-dimensional triangulations ******/
  std::cout << "    insertions 2-dim" << std::endl;

  Triangul T2_3;
  Vertex_handle v2_3_1 = T2_3.insert(p1);
  Vertex_handle v2_3_2 = T2_3.insert(p2);
  Vertex_handle v2_3_3 = T2_3.insert(p3);
  Vertex_handle v2_3_8 = T2_3.insert(p8);
  Vertex_handle v2_3_9 = T2_3.insert(p9);
  assert( T2_3.dimension() == 1 );
  Vertex_handle v2_3_4 = T2_3.insert(p4);
  assert( T2_3.dimension() == 2 );
  Vertex_handle v2_3_6 = T2_3.insert(p6, T2_3.finite_faces_begin());
  Vertex_handle v2_3_0 = T2_3.insert(p0, ++T2_3.finite_faces_begin());
  Vertex_handle v2_3_5 = T2_3.insert(p5);
  Vertex_handle v2_3_7 = T2_3.insert(p7);
  loc = T2_3.locate(p10,lt,li);
  Vertex_handle v2_3_10 = T2_3.insert(p10, lt, loc, li);
  assert( lt == Triangul::FACE );
  assert( T2_3.dimension() == 2 );
  assert( T2_3.number_of_vertices() == 11 );
  // assert( T2_3.number_of_faces() == 13 );
  assert( T2_3.is_valid() );

  // make sure inserting on a previous point does not insert it again
  assert( T2_3.insert(p10) == v2_3_10 );
  assert( T2_3.number_of_vertices() == 11 );

  // make sure push_back exists and does the same thing as insert
  assert( T2_3.push_back(p10) == v2_3_10 );
  assert( T2_3.number_of_vertices() == 11 );

  // test generic iterator insert
#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  Triangul T2_4; T2_4.insert( Point_iterator(T2_3.finite_vertices_begin()),
                         Point_iterator(T2_3.finite_vertices_end()) );
  assert( T2_3.dimension() == 2 );
  assert( T2_3.number_of_vertices() == 11 );
  // assert( T2_4.number_of_faces() == 13 );
  assert( T2_3.is_valid() );
#endif

  // test list iterator insert
  Triangul T2_5;
  assert( T2_5.insert(l.begin(), l.end()) == 4 );
  assert( T2_5.dimension() == 2 );
  assert( T2_5.number_of_vertices() == 4 );
  // assert( T2_5.number_of_faces() == 13 );
  assert( T2_5.is_valid() );

  // test list iterator insert
  Triangul T2_6;
  assert( T2_6.insert(v.begin(), v.end()) == 4 );
  assert( T2_6.dimension() == 2 );
  assert( T2_6.number_of_vertices() == 4 );
  // assert( T2_6.number_of_faces() == 13 );
  assert( T2_6.is_valid() );


  // test flip
     std::cout << "    test flip " << std::endl;
     Triangul T2_8;
     T2_8.insert(Point(0,0,1));
     T2_8.insert(Point(1,0,1));
     T2_8.insert(Point(1,1,1));
     T2_8.insert(Point(0,1,1));
     Face_handle ff = T2_8.locate(Point(1,1,2),lt,li);
     assert(lt == Triangul::EDGE);
     assert(!T2_8.is_infinite(ff));
     Face_handle f2 = ff->neighbor(li);
     assert(!T2_8.is_infinite(f2));
     int fli = ff->index(f2);
     T2_8.flip(ff,fli);
     assert( T2_8.is_valid() );


  /****************************/
  /***** CONSTRUCTORS (2) *****/
  std::cout << "    constructors (2)" << std::endl;

  // test copy_constructor with non-empty 0-triangulation
  Triangul T0_1_1( T0_1 );
  assert( T0_1_1.dimension() == 0 );
  assert( T0_1_1.number_of_vertices() == 1 );
  assert( T0_1_1.is_valid() );

  // test copy_constructor with non-empty 1-triangulation
  Triangul T1_5_1( T1_5 );
  assert( T1_5_1.dimension() == 1 );
  assert( T1_5_1.number_of_vertices() == 5 );
  assert( T1_5_1.is_valid() );

  // Test assignment operator
  Triangul T1_5_2 = T1_5;
  assert( T1_5_2.dimension() == 1 );
  assert( T1_5_2.number_of_vertices() == 5 );
  assert( T1_5_2.is_valid() );

  // test copy_constructor with non-empty 2-triangulation
  Triangul T2_8_1(T2_8);
  assert( T2_8_1.is_valid());

  // test assignment operator
  Triangul T2_3_4 = T2_3;
  assert( T2_3_4.dimension() == 2 );
  assert( T2_3_4.number_of_vertices() == 11 );
  assert( T2_3_4.is_valid() );

  /****************************/
  /******* MOVE SEMANTIC*******/

  std::cout << "    move constructors and move assignment" << std::endl;
  namespace test_tr_23 = CGAL::Testsuite::Triangulation_23;
  test_tr_23::test_move_semantic(T0_1);
  test_tr_23::test_move_semantic(T1_5);
  test_tr_23::test_move_semantic(T2_8);
  test_tr_23::test_move_semantic(T2_3);

  /*********************************************/
  /****** FINITE/INFINITE VERTICES/FACES *******/

  std::cout << "    finite/infinite vertices/faces" << std::endl;
  _test_fct_is_infinite( T0_0 );
  _test_fct_is_infinite( T0_1 );
  _test_fct_is_infinite( T1_5 );
  _test_fct_is_infinite( T2_3 );


  /*************************************/
  /******** POINT LOCATIONS ************/

  // Locate_type lt; // see above

  // Check point location in 0-dimensional triangulations
  // No need because of precondition (at least two vertices)

  // Check point location in 1-dimensional triangulations
  std::cout << "    point locations 1-dim" << std::endl;
  Triangul T1_3_2;
  T1_3_2.insert(p1);
  T1_3_2.insert(p2);
  T1_3_2.insert(p9);
  f = T1_3_2.locate(p1,lt,li); assert( lt == Triangul::VERTEX );
  assert( T1_3_2.xy_equal(f->vertex(li)->point(), p1) );
  f = T1_3_2.locate(p2,lt,li); assert( lt == Triangul::VERTEX );
  assert( T1_3_2.xy_equal(f->vertex(li)->point(), p2) );
  f = T1_3_2.locate(p9,lt,li); assert( lt == Triangul::VERTEX );
  assert( T1_3_2.xy_equal(f->vertex(li)->point(), p9) );
  f = T1_3_2.locate(p3,lt,li); assert( lt == Triangul::EDGE );
  assert( (T1_3_2.xy_equal(f->vertex(f->ccw(li))->point(), p1)
        && T1_3_2.xy_equal(f->vertex(f->cw(li))->point(), p2))
       || (T1_3_2.xy_equal(f->vertex(f->ccw(li))->point(), p2)
        && T1_3_2.xy_equal(f->vertex(f->cw(li))->point(), p1)));
  f = T1_3_2.locate(p8,lt,li); assert( lt == Triangul::OUTSIDE_CONVEX_HULL );
    assert(T1_3_2.is_infinite(f->vertex(li)));
  f = T1_3_2.locate(p0,lt,li); assert( lt == Triangul::OUTSIDE_AFFINE_HULL );
  f = T1_3_2.locate(p7,lt,li); assert( lt == Triangul::OUTSIDE_AFFINE_HULL );
  f = T1_3_2.locate(p5,lt,li); assert( lt == Triangul::OUTSIDE_AFFINE_HULL );
  f = T1_3_2.locate(p4,lt,li); assert( lt == Triangul::OUTSIDE_AFFINE_HULL );
  f = T1_3_2.locate(p6,lt,li); assert( lt == Triangul::OUTSIDE_AFFINE_HULL );


  // Check point location in 2-dimensional triangulations
  std::cout << "    point locations 2-dim" << std::endl;
  f = T2_3.locate(p0,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_3.xy_equal(f->vertex(li)->point(), p0) );
  f = T2_3.locate(p1,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_3.xy_equal(f->vertex(li)->point(), p1) );
  f = T2_3.locate(p2,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_3.xy_equal(f->vertex(li)->point(), p2) );
  f = T2_3.locate(p3,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_3.xy_equal(f->vertex(li)->point(), p3) );
  f = T2_3.locate(p4,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_3.xy_equal(f->vertex(li)->point(), p4) );
  f = T2_3.locate(p5,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_3.xy_equal(f->vertex(li)->point(), p5) );
  f = T2_3.locate(p6,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_3.xy_equal(f->vertex(li)->point(), p6) );
  f = T2_3.locate(p7,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_3.xy_equal(f->vertex(li)->point(), p7) );
  f = T2_3.locate(p8,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_3.xy_equal(f->vertex(li)->point(), p8) );
  f = T2_3.locate(p9,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_3.xy_equal(f->vertex(li)->point(), p9) );
  f = T2_3.locate(p10,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_3.xy_equal(f->vertex(li)->point(), p10) );
  f = T2_3.locate(p11,lt,li); assert( lt == Triangul::EDGE );
  assert( (T2_3.xy_equal(f->vertex(f->ccw(li))->point(), p1)
        && T2_3.xy_equal(f->vertex(f->cw(li))->point(), p0))
       || (T2_3.xy_equal(f->vertex(f->ccw(li))->point(), p0)
        && T2_3.xy_equal(f->vertex(f->cw(li))->point(), p1)));
  f = T2_3.locate(p12,lt,li); assert( lt == Triangul::FACE );
  assert( T2_3.oriented_side(f,p12) == CGAL::ON_POSITIVE_SIDE );
  f = T2_3.locate(p13,lt,li,f); assert( lt == Triangul::OUTSIDE_CONVEX_HULL );
  assert( T2_3.orientation(p13,
                           f->vertex(f->ccw(li))->point(),
                           f->vertex(f->cw(li))->point())
            == CGAL::COUNTERCLOCKWISE);
  f = T2_3.locate(p13,lt,li); assert( lt == Triangul::OUTSIDE_CONVEX_HULL );
  f = T2_3.locate(p14,lt,li); assert( lt == Triangul::OUTSIDE_CONVEX_HULL );
  assert( T2_3.orientation(p14,
                          f->vertex(f->ccw(li))->point(),
                          f->vertex(f->cw(li))->point())
            == CGAL::COUNTERCLOCKWISE);
  f = T2_3.locate(p15,lt,li); assert( lt == Triangul::OUTSIDE_CONVEX_HULL );
  assert( T2_3.orientation(p15,
                           f->vertex(f->ccw(li))->point(),
                           f->vertex(f->cw(li))->point())
            == CGAL::COUNTERCLOCKWISE);



  /*************************/
  /******* Iterators *******/
  std::cout << "    iterators" << std::endl;
  _test_iterators(T0_0);
  _test_iterators(T0_1);
  _test_iterators(T1_5);
  _test_iterators(T2_3);


  /***************************/
  /******* Circulators *******/
  std::cout << "    circulators" << std::endl;
  _test_circulators(T0_0);
  _test_circulators(T0_1);
  _test_circulators(T1_5);
  _test_circulators(T2_3);

  // Line_face_circulator
  std::cout << "    line face circulator  " << std::endl;
  typedef typename Triangul::Line_face_circulator LFC;
  // here == operator needed for Point!
  // testing with the grid triangulation
  LFC fc= T2_3.line_walk(p1,p10);
  assert(fc != nullptr);
  assert(!fc.is_empty());
  LFC fc2=fc;
  assert(fc==fc2);
  fc++;
  fc--;
  ++fc;
  --fc;
  Point pp(0,1,2); //Point pp(0,0.5);
  f = T2_3.locate(pp,lt,li);
  //assert(lt==Triangul::FACE); no longer true, I've change the triangulation;
  fc= T2_3.line_walk(pp,p10,f);
  fc2=fc;
  assert(fc==fc2);
  fc++;
  fc--;
  ++fc;
  --fc;

  // testing with dummy triangulations
  int n=0;
  do {fc2++ ; n = n+1;} while (fc2 != fc);
  assert(T2_8.number_of_vertices()>=2);
  assert(T2_8.is_valid());
  fc= T2_8.line_walk(Point(1,4,10),Point(20,5,10));
  fc2=fc;
  n=0;
  assert(fc==fc2);
  do {fc2++ ; n = n+1;} while (fc2 != fc);
  assert(n==4);
  // the two point are vertices of the triangulation.
  Triangul TT;
  TT.insert(Point(0,0)); TT.insert(Point(1,0));TT.insert(Point(1,1));
  TT.insert(Point(0,1));
  assert(TT.dimension()==2);
  assert(TT.is_valid());
  assert(TT.number_of_vertices()==4);
  f = TT.locate(Point(0,0));
  fc = TT.line_walk(Point(0,0),Point(1,1));
  fc2 = TT.line_walk(Point(0,0),Point(1,1),f);
  assert(fc==fc2);
  n=0;
  do {fc2++ ; n = n+1;} while (fc2 != fc);
  assert(n==4 || n==3);

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
  f = T2_3.locate(p12,lt,li); // from section locate above
  Triangle t = T2_3.triangle(f); assert( &t == &t );
  Segment  s = T2_3.segment(f,0); assert( &s == &s );
  s = T2_3.segment(Edge(f,1)); assert( &s == &s );
  s = T2_3.segment(T2_3.incident_edges(v2_3_6)); assert( &s == &s );
  s = T2_3.segment(T2_3.finite_edges_begin()); assert( &s == &s );


  /********************/
  /******** I/O *******/
  std::string T15fname, T23fname;
  std::stringstream ss15;
  ss15 << "T15" << &T1_5 << ".triangulation";
  ss15 >> T15fname;
  std::stringstream ss23;
  ss23 << "T23" << &T2_3 << ".triangulation";
  ss23 >> T23fname;

  std::cout << "    output to a file" << std::endl;
  std::ofstream of1_5(T15fname.c_str());
  CGAL::IO::set_ascii_mode(of1_5);
  of1_5 << T1_5; of1_5.close();
  std::ofstream of2_3(T23fname.c_str());
  CGAL::IO::set_ascii_mode(of2_3);
  of2_3 << T2_3; of2_3.close();


  std::cout << "    input from a file" << std::endl;

  std::ifstream if1_5(T15fname.c_str()); CGAL::IO::set_ascii_mode(if1_5);
  Triangul T1_5_copy; if1_5 >> T1_5_copy;
  assert( T1_5_copy.is_valid() &&
          T1_5_copy.number_of_vertices() == T1_5.number_of_vertices() );

  std::ifstream if2_3(T23fname.c_str()); CGAL::IO::set_ascii_mode(if2_3);
  Triangul T2_3_copy; if2_3 >> T2_3_copy;
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
  assert( T1_5.is_valid());
  T1_5.remove(v1_5_2);
  T1_5.remove(v1_5_3);
  T1_5.remove(v1_5_8);
  T1_5.remove(v1_5_9);
  assert( T1_5.number_of_vertices() == 0 );

  // remove from 2-dimensional triangulations
  T2_3.remove(v2_3_10); assert( T2_3.is_valid());
  T2_3.remove(v2_3_9);
  T2_3.remove(v2_3_8);
  T2_3.remove(v2_3_7);
  T2_3.remove(v2_3_6);
  T2_3.remove(v2_3_5);
  T2_3.remove(v2_3_4);
  T2_3.remove(v2_3_3); assert(T2_3.is_valid());
  T2_3.remove(v2_3_0);
  assert(T2_3.is_valid());
  T2_3.remove(v2_3_2); assert(T2_3.is_valid());
  T2_3.remove(v2_3_1); assert(T2_3.is_valid());
  assert( T2_3.number_of_vertices() == 0 );

  typename Triangul::size_type i;
  for (i=T2_4.number_of_vertices(); i>0; i--)
    T2_4.remove(T2_4.finite_vertex());
  assert( T2_4.number_of_vertices() == 0 );

  T2_5.clear();
  assert( T2_5.number_of_vertices() == 0 );


  // test destructors and return
  std::cout << "    test destructors and return" << std::endl;
}
