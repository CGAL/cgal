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

#include <utility>
#include <list>

#include <CGAL/_test_fct_is_infinite.h>
#include <CGAL/_test_triangulation_iterators.h>
#include <CGAL/_test_triangulation_circulators.h>
#include <CGAL/_test_line_face_circulator.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/use.h>

template <class Triangul>
void
_test_cls_triangulation_2( const Triangul & )
{
  //typedef Triangulation                      Cls;

  // We assume the traits class has been tested already
  // actually, any traits is good if it has been tested
  typedef typename Triangul::Geom_traits          Gt;

  typedef typename Triangul::size_type            size_type;
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
  typedef typename Triangul::Finite_faces_iterator     Finite_faces_iterator;
  typedef typename Triangul::Finite_edges_iterator     Finite_edges_iterator;
  typedef typename Triangul::Point_iterator            Point_iterator;

  typedef typename Triangul::Vertex_circulator    Vertex_circulator;
  typedef typename Triangul::Face_circulator      Face_circulator;
  typedef typename Triangul::Edge_circulator      Edge_circulator;
  typedef typename Triangul::Line_face_circulator Line_face_circulator;

  typedef typename Triangul::Locate_type          Locate_type;
  
  CGAL_USE_TYPE(Gt);
  CGAL_USE_TYPE(Vertex);
  CGAL_USE_TYPE(Face);
  CGAL_USE_TYPE(Finite_vertices_iterator);
  CGAL_USE_TYPE(Finite_edges_iterator);
  CGAL_USE_TYPE(Vertex_circulator);
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

  int px=1, py=1;
  int qx=-1, qy=2;
  Locate_type lt;
  int li;
  Face_handle    f, ff, loc;
  Vertex_handle v0, v1, v2, v3;


  std::list<Point> l; l.push_back(p0);
  l.push_back(p1); l.push_back(p2); l.push_back(p3);
  l.push_back(p4); l.push_back(p5); l.push_back(p6);
  l.push_back(p7); l.push_back(p8); l.push_back(p9);
   
  std::vector<Point> v; v.push_back(p0);
  v.push_back(p1); v.push_back(p2); v.push_back(p3);
  v.push_back(p4); v.push_back(p5); v.push_back(p6);
  v.push_back(p7); v.push_back(p8); v.push_back(p9);
  
  
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
  assert( T0_0.number_of_faces() == 0);
  assert( T0_0.is_valid() );

  Triangul T0_1; 
  Vertex_handle v0_1_0 = T0_1.insert(p0); assert( v0_1_0 != NULL );
  assert( T0_1.dimension() == 0 );
  assert( T0_1.number_of_vertices() == 1 );
  assert( T0_1.number_of_faces() == 0);
  assert( T0_1.is_valid() );

  // test insert_first()
  Triangul T0_2; 
  Vertex_handle v0_2_0 =   T0_2.insert_first(p0);
  assert( v0_2_0 != NULL );
  assert( T0_2.dimension() == 0 );
  assert( T0_2.number_of_vertices() == 1 );
  assert( T0_2.number_of_faces() == 0);
  assert( T0_2.is_valid() );
  
  /******** 1-dimensional triangulations ******/
  // T1_n denotes a 1-dimensional triangulation with n vertices
  // when there are several, we use T1_n_p
  std::cout << "    insertions 1-dim" << std::endl;
  
  Triangul T1_2;
  Vertex_handle v1_2_1 = T1_2.insert(p1);
  Vertex_handle v1_2_2 = T1_2.insert(p2);
  assert( T1_2.dimension() == 1 );
  assert( T1_2.number_of_vertices() == 2 );
  assert( T1_2.number_of_faces() == 0 );
  assert( T1_2.is_valid() );
  
  // p1,p3,p2  [endpoints first]
  Triangul T1_3_0;
  Vertex_handle v1_3_0_1 = T1_3_0.insert(p1); assert( v1_3_0_1 != NULL );
  Vertex_handle v1_3_0_3 = T1_3_0.insert(p3); assert( v1_3_0_3 != NULL );
  Vertex_handle v1_3_0_2 = T1_3_0.insert(p2); assert( v1_3_0_2 != NULL );
  assert( T1_3_0.dimension() == 1 );
  assert( T1_3_0.number_of_vertices() == 3 );
  assert( T1_3_0.number_of_faces() == 0 );
  assert( T1_3_0.is_valid() );
  
  // p1,p2,p3  [middle point first]
  Triangul T1_3_1;
  Vertex_handle v1_3_1_1 = T1_3_1.insert(p1); assert( v1_3_1_1 != NULL );
  Vertex_handle v1_3_1_3 = T1_3_1.insert(p2); assert( v1_3_1_3 != NULL );
  Vertex_handle v1_3_1_2 = T1_3_1.insert(p3); assert( v1_3_1_2 != NULL );
  assert( T1_3_1.dimension() == 1 );
  assert( T1_3_1.number_of_vertices() == 3 );
  assert( T1_3_1.number_of_faces() == 0 );
  assert( T1_3_1.is_valid() );

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

  // test insert_second()
  Triangul T1_6 = T0_2; 
  Vertex_handle v1_6_2 = T1_6.insert_second(p3); assert( v1_6_2 != NULL );
  assert( T1_6.dimension() == 1 );
  assert( T1_6.number_of_vertices() == 2 );
  assert( T1_6.is_valid() ); 
  
  /******** 2-dimensional triangulations ******/ 
  std::cout << "    insertions 2-dim" << std::endl;
  
  Triangul T2_1;
  Vertex_handle v2_1_0 = T2_1.insert(p0);
  Vertex_handle v2_1_1 = T2_1.insert(p1);
  Vertex_handle v2_1_2 = T2_1.insert(p2);
  Vertex_handle v2_1_3 = T2_1.insert(p3); // on the edge p1,p2, on the convex hull
  Vertex_handle v2_1_4 = T2_1.insert(p4); // outside, with two visible collineaar edges
  Vertex_handle v2_1_5 = T2_1.insert(p5); 
  Vertex_handle v2_1_6 = T2_1.insert(p6);  // outside, collinear with p2,p5
  Vertex_handle v2_1_7 = T2_1.insert(p7);  // outside with two visible collinear edges
                      		           // but also collinear with and extending p0,p5 
  Vertex_handle v2_1_8 = T2_1.insert(p8); 
  Vertex_handle v2_1_9 = T2_1.insert(p9);    // inside, on the edge p6,p7
  Vertex_handle v2_1_10 = T2_1.insert(p10);  // inside the face p2,p4,p6
  assert( T2_1.dimension() == 2 );
  assert( T2_1.number_of_vertices() == 11 );
#ifndef CGAL_NO_DEPRECATED_CODE
  assert( T2_1.number_of_faces() == 2 * 12 - 4 
                                 - T2_1.infinite_vertex()->degree() );
#endif

  assert( T2_1.number_of_faces() == 2 * 12 - 4 
                                 - T2_1.degree(T2_1.infinite_vertex()) );
  
  // test is_valid for 2-triangulations
  assert( T2_1.is_valid() );

  // we now test the other insert functions
  // more vicious, we insert all the points on a single line first
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
  assert( T2_3.is_valid() );
  
  // make sure inserting on a previous point does not insert it again
  assert( T2_3.insert(p10) == v2_3_10 );
  assert( T2_3.number_of_vertices() == 11 );

  // make sure push_back exists and does the same thing as insert
  assert( T2_3.push_back(p10) == v2_3_10 );
  assert( T2_3.number_of_vertices() == 11 );

  // test generic iterator insert
#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  Triangul T2_4; T2_4.insert( Point_iterator(T2_1.finite_vertices_begin()),
                         Point_iterator(T2_1.finite_vertices_end()) );
  assert( T2_4.dimension() == 2 );
  assert( T2_4.number_of_vertices() == 11 );
  assert( T2_4.is_valid() );
#endif

  // test list iterator insert
  Triangul T2_5;
  assert( T2_5.insert(l.begin(), l.end()) == 10 );
  assert( T2_5.dimension() == 2 );
  assert( T2_5.number_of_vertices() == 10 );
  assert( T2_5.is_valid() );

  // test list iterator insert
  Triangul T2_6;
  assert( T2_6.insert(v.begin(), v.end()) == 10 );
  assert( T2_6.dimension() == 2 );
  assert( T2_6.number_of_vertices() == 10 );
  assert( T2_6.is_valid() );
  
  // test grid insert
  Triangul T2_7;
  int m, p;
  for (m=0; m<3; m++)
    for (p=0; p<3; p++)
      T2_7.insert( Point(m*px+p*qx, m*py+p*qy, 1) );
  assert( T2_7.number_of_vertices() == static_cast<size_type> (m*p) );
  assert( T2_7.is_valid() );

  // test flip
     std::cout << "    test flip " << std::endl;
     Triangul T2_8;
     T2_8.insert(Point(0,0,1));
     T2_8.insert(Point(1,0,1));
     T2_8.insert(Point(1,1,1));
     T2_8.insert(Point(0,1,1));
     ff = T2_8.locate(Point(1,1,2),lt,li);
     assert(lt == Triangul::EDGE);
     assert(!T2_8.is_infinite(ff));
     Face_handle f2 = ff->neighbor(li);
     assert(!T2_8.is_infinite(f2));
     int fli = ff->index(f2);
     T2_8.flip(ff,fli);
     assert( T2_8.is_valid() );
     
  //make_hole star_hole
  std::list<Edge> hole;
  T2_3.make_hole(v2_3_10, hole);
  T2_3.delete_vertex(v2_3_10);
  v2_3_10 = T2_3.star_hole(p10, hole.begin(), hole.end());
  assert( T2_3.is_valid());

  /**************************/
  /******* MOVE *********/
  std::cout << "    displacements" << std::endl;

  std::cout << "    degenerate cases: " << std::endl;
  
  Triangul TM_0, TM_1;
  Vertex_handle tmv1 = TM_0.insert(Point(0,0));
  Vertex_handle tmv2 = TM_0.insert(Point(1,0));
  Vertex_handle tmv3 = TM_0.insert(Point(2,0));
  Vertex_handle tmv4 = TM_0.insert(Point(1,1));
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv4, Point(2, 1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv4, Point(2, -1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv4, Point(3, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  TM_0.move_if_no_collision(tmv3, Point(1, 1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv3, Point(-1, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  TM_0.move_if_no_collision(tmv2, Point(-1, 0, 2));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  TM_0.move_if_no_collision(tmv2, Point(-1, 0, 4));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  TM_0.move_if_no_collision(tmv2, Point(-1, 0, 2));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  TM_0.move_if_no_collision(tmv2, Point(-1, 1, 2));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv1, Point(0, 2));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv1, Point(0, 1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv1, Point(0, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);
  assert(TM_0.move_if_no_collision(tmv1, Point(3, 0)) != tmv1);

  TM_0.move_if_no_collision(tmv1, Point(0, 1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv4, Point(1, 2));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);   

  TM_0.move_if_no_collision(tmv4, Point(3, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv1, Point(2, 3));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv4, Point(1, 2));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  std::cout << "    non-degenerate cases: " << std::endl;
  // non-degenerate cases
  std::list<Point> points;
  for(int count=0; count<50; count++) {
    points.push_back(Point(rand()%30000, rand()%30000));
  }
  TM_1.insert(points.begin(), points.end());
	Vertex_handle vTM_1;
  for(int i=0; i<5; i++) {
    for(typename Triangul::Finite_vertices_iterator 
         fvi = TM_1.finite_vertices_begin();
         fvi != TM_1.finite_vertices_end(); fvi++) {
      Point p = Point(rand()%30000, rand()%30000);
      vTM_1 = TM_1.move_if_no_collision(fvi, p);
      assert(TM_1.is_valid());
    }
  }

  // A simple test to see if move returns the good vertex
  // when there is a collision
  assert(TM_1.move(TM_1.finite_vertices_begin(), vTM_1->point()) == vTM_1);

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
  Triangul T2_1_1( T2_1 );
  assert( T2_1_1.dimension() == 2 );
  assert( T2_1_1.number_of_vertices() == 11 );
  assert( T2_1_1.is_valid() );
  
  // test assignment operator
  Triangul T2_1_4 = T2_1;
  assert( T2_1_4.dimension() == 2 );
  assert( T2_1_4.number_of_vertices() == 11 );
  assert( T2_1_4.is_valid() );
  

  /*********************************************/
  /****** FINITE/INFINITE VERTICES/FACES *******/

  std::cout << "    finite/infinite vertices/faces" << std::endl;
  _test_fct_is_infinite( T0_0 );
  _test_fct_is_infinite( T0_1 );
  _test_fct_is_infinite( T1_2 );
  _test_fct_is_infinite( T1_5 );
  _test_fct_is_infinite( T2_1 );
  _test_fct_is_infinite( T2_3 );
  _test_fct_is_infinite( T2_4 );
  _test_fct_is_infinite( T2_5 );
  _test_fct_is_infinite( T2_6 );

  std::cout << "    is_edge, is_face, includes_edge" << std::endl;
  assert(  T1_5.is_edge(v1_5_1,v1_5_3));
  assert( !T1_5.is_edge(v1_5_1,v1_5_2));
  assert(  T1_5.is_edge(v1_5_1,v1_5_3,f,li));
  assert(! T1_5.is_edge(v1_5_1,v1_5_2,f,li));
  assert(  T2_1.is_edge(v2_1_1,v2_1_3));
  assert( !T2_1.is_edge(v2_1_1,v2_1_2));
  assert(  T2_1.is_edge(v2_1_1,v2_1_3,f,li));
  assert(! T2_1.is_edge(v2_1_1,v2_1_2,f,li));
  f = T2_1.finite_faces_begin();
  v0 = f->vertex(0);
  v1 = f->vertex(1);
  v2 = f->vertex(2);
  v3 = T2_1.mirror_vertex(f,0);
  assert(T2_1.is_face(v0,v1,v2));
  assert(T2_1.is_face(v0,v2,v1));
  assert(T2_1.is_face(v1,v2,v3,ff) && ff == f->neighbor(0));
  const Edge e = T2_1.mirror_edge(Edge(f, 0));
  assert(e.first == ff);
  assert(e.second == ff->index(v3));
  assert(! T2_1.is_face(v0,v3,v1));
  assert(T1_5.includes_edge(v1_5_1,v1_5_2,v0,f,li));
  assert(v0 == v1_5_3);
  assert(T2_1.includes_edge(v2_1_1,v2_1_2,v0,f,li));
  assert(v0 == v2_1_3);

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
  f = T2_1.locate(p0,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_1.xy_equal(f->vertex(li)->point(), p0) );
  f = T2_1.locate(p1,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_1.xy_equal(f->vertex(li)->point(), p1) );
  f = T2_1.locate(p2,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_1.xy_equal(f->vertex(li)->point(), p2) );
  f = T2_1.locate(p3,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_1.xy_equal(f->vertex(li)->point(), p3) );
  f = T2_1.locate(p4,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_1.xy_equal(f->vertex(li)->point(), p4) );
  f = T2_1.locate(p5,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_1.xy_equal(f->vertex(li)->point(), p5) );
  f = T2_1.locate(p6,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_1.xy_equal(f->vertex(li)->point(), p6) );
  f = T2_1.locate(p7,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_1.xy_equal(f->vertex(li)->point(), p7) );
  f = T2_1.locate(p8,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_1.xy_equal(f->vertex(li)->point(), p8) );
  f = T2_1.locate(p9,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_1.xy_equal(f->vertex(li)->point(), p9) );
  f = T2_1.locate(p10,lt,li); assert( lt == Triangul::VERTEX );
  assert( T2_1.xy_equal(f->vertex(li)->point(), p10) );
  f = T2_1.locate(p11,lt,li); assert( lt == Triangul::EDGE );
  assert( (T2_1.xy_equal(f->vertex(f->ccw(li))->point(), p1)
        && T2_1.xy_equal(f->vertex(f->cw(li))->point(), p0))
       || (T2_1.xy_equal(f->vertex(f->ccw(li))->point(), p0)
        && T2_1.xy_equal(f->vertex(f->cw(li))->point(), p1)));
  f = T2_1.locate(p12,lt,li); assert( lt == Triangul::FACE );
  assert( T2_1.oriented_side(f,p12) == CGAL::ON_POSITIVE_SIDE );
  f = T2_1.locate(p13,lt,li,f); assert( lt == Triangul::OUTSIDE_CONVEX_HULL );
  assert( T2_1.orientation(p13,
			   f->vertex(f->ccw(li))->point(),
			   f->vertex(f->cw(li))->point())
	    == CGAL::COUNTERCLOCKWISE);
  f = T2_1.locate(p14,lt,li); assert( lt == Triangul::OUTSIDE_CONVEX_HULL );
  assert( T2_1.orientation(p14,
			  f->vertex(f->ccw(li))->point(),
			  f->vertex(f->cw(li))->point())
	    == CGAL::COUNTERCLOCKWISE);
  f = T2_1.locate(p15,lt,li); assert( lt == Triangul::OUTSIDE_CONVEX_HULL );
  assert( T2_1.orientation(p15,
			   f->vertex(f->ccw(li))->point(),
			   f->vertex(f->cw(li))->point())
	    == CGAL::COUNTERCLOCKWISE);

  // test grid locate
  for (m=0; m<1; m++)
    for (p=0; p<1; p++)
      {
	Point q= Point(m*px+p*qx, m*py+p*qy, 1);
       	f = T2_7.locate(q,lt,li); assert( lt == Triangul::VERTEX );
  	assert( T2_7.xy_equal(f->vertex(li)->point(), q) );
      }
  for (m=0; m<1; m++)
    for (p=0; p<1; p++)
      {
	Point q= Point(2*m*px+(2*p+1)*qx, 2*m*py+(2*p+1)*qy, 2);
	Point r= Point(m*px+p*qx, m*py+p*qy, 1);
	Point s= Point(m*px+(p+1)*qx, m*py+(p+1)*qy, 1);
       	f = T2_7.locate(q,lt,li); assert( lt == Triangul::EDGE );
        assert( (T2_7.xy_equal(f->vertex(f->ccw(li))->point(), r)
              && T2_7.xy_equal(f->vertex(f->cw(li))->point(), s))
             || (T2_7.xy_equal(f->vertex(f->ccw(li))->point(), s)
              && T2_7.xy_equal(f->vertex(f->cw(li))->point(), r)));

      }
  for (m=0; m<1; m++)
    for (p=0; p<1; p++)
      {
	Point q= Point((50*m+1)*px+(50*p+1)*qx, (50*m+1)*py+(50*p+1)*qy, 50);
       	f = T2_7.locate(q,lt,li); assert( lt == Triangul::FACE );
	assert( T2_7.oriented_side(f,q) == CGAL::ON_POSITIVE_SIDE );
      }

  /*************************/
  /******* Iterators *******/
  std::cout << "    iterators" << std::endl;
  _test_iterators(T0_0);
  _test_iterators(T0_1);
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

  //test iterators as arguments
  Finite_faces_iterator fit = T2_7.finite_faces_begin();
  assert(!T2_7.is_infinite(fit));
  while(!T2_7.is_infinite(fit->neighbor(0)) ) ++fit;
 

  /***************************/
  /******* Circulators *******/
  std::cout << "    circulators" << std::endl;
  _test_circulators(T0_0);
  _test_circulators(T0_1);
  _test_circulators(T1_2);
  _test_circulators(T1_3_0);
  _test_circulators(T1_3_1);
  _test_circulators(T1_5);
  _test_circulators(T1_6);
  _test_circulators(T2_1);
  _test_circulators(T2_3);
  _test_circulators(T2_5);
  _test_circulators(T2_6);
  _test_circulators(T2_7);
  
  // Line_face_circulator
  std::cout << "    line face circulator  " << std::endl;
  _test_line_face_circulator(Triangul());
  typedef typename Triangul::Line_face_circulator LFC;
  // here == operator needed for Point!
  // testing with the grid triangulation
  LFC fc= T2_7.line_walk(p1,p10);
  assert(fc!=NULL);
  assert(!fc.is_empty());
  LFC fc2=fc;
  assert(fc==fc2);
  fc++;
  fc--;
  ++fc;
  --fc;
  Point pp(0,1,2); //Point pp(0,0.5);
  f = T2_7.locate(pp,lt,li);
  assert(lt==Triangul::FACE);
  fc= T2_7.line_walk(pp,p10,f);
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
  //fc= T2_8.line_walk(Point(0.5,0.4),Point(5,5));
  fc= T2_8.line_walk(Point(5,4,10),Point(5,5));
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
  f = TT.locate(Point(0,0),lt,li);
  assert(lt == Triangul::VERTEX);
  Face_circulator fcc = TT.incident_faces(f->vertex(li),f);
  while( TT.is_infinite(fcc)) ++fcc;
  fc = TT.line_walk(Point(0,0),Point(1,1));
  fc2 = TT.line_walk(Point(0,0),Point(1,1),fcc);
  assert(fc==fc2);
  //fc.print();
  n=0;
  do { 
    //fc2.print(); 
    fc2++ ; 
    n = n+1;} while (fc2 != fc);
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
  f = T2_1.locate(p12,lt,li); // from section locate above
  Triangle t = T2_1.triangle(f); assert( &t == &t );
  Segment  s = T2_1.segment(f,0); assert( &s == &s );
  s = T2_1.segment(Edge(f,1)); assert( &s == &s );
#ifndef CGAL_NO_DEPRECATED_CODE
  s = T2_1.segment(v2_1_6->incident_edges()); assert( &s == &s );
#endif
  s = T2_1.segment(T2_1.incident_edges(v2_1_6)); assert( &s == &s );
  s = T2_1.segment(T2_1.finite_edges_begin()); assert( &s == &s );

  // finite/infinite vertex
  // T2_1.set_infinite_vertex(T2_1.infinite_vertex());
  // assert( T2_1.is_valid() );

  /********************/
  /******** I/O *******/
  std::cout << "    output to a file" << std::endl;
  std::ofstream of0_0("T00.triangulation", std::ios::out);
  CGAL::set_ascii_mode(of0_0); 
  of0_0 << T0_0; of0_0.close();
  std::ofstream of0_1("T01.triangulation");
  CGAL::set_ascii_mode(of0_1); 
  of0_1 << T0_1; of0_1.close();
  std::ofstream of1_2("T12.triangulation");
  CGAL::set_ascii_mode(of1_2); 
  of1_2 << T1_2; of1_2.close();
  std::ofstream of1_5("T15.triangulation");
  CGAL::set_ascii_mode(of1_5); 
  of1_5 << T1_5; of1_5.close();
  std::ofstream of1_6("T16.triangulation");
  CGAL::set_ascii_mode(of1_6); 
  of1_6 << T1_6; of1_6.close();
  std::ofstream of2_1("T21.triangulation");
  CGAL::set_ascii_mode(of2_1); 
  of2_1 << T2_1; of2_1.close();
  std::ofstream of2_3("T23.triangulation");
  CGAL::set_ascii_mode(of2_3); 
  of2_3 << T2_3; of2_3.close();
  std::ofstream of2_5("T25.triangulation");
  CGAL::set_ascii_mode(of2_5); 
  of2_5 << T2_5; of2_5.close();
  std::ofstream of2_6("T26.triangulation");
  CGAL::set_ascii_mode(of2_6); 
  of2_6 << T2_6; of2_6.close();

  std::cout << "    input from a file" << std::endl;
  std::ifstream if0_0("T00.triangulation"); CGAL::set_ascii_mode(if0_0);
  Triangul T0_0_copy;   if0_0 >> T0_0_copy;
  assert( T0_0_copy.is_valid() &&
	  T0_0_copy.number_of_vertices() == T0_0.number_of_vertices() );
  std::ifstream if0_1("T01.triangulation"); CGAL::set_ascii_mode(if0_1);
  Triangul T0_1_copy; if0_1 >> T0_1_copy;
  assert( T0_1_copy.is_valid() &&
	  T0_1_copy.number_of_vertices() == T0_1.number_of_vertices() );
  std::ifstream if1_2("T12.triangulation"); CGAL::set_ascii_mode(if1_2); 
  Triangul T1_2_copy; if1_2 >> T1_2_copy;
  assert( T1_2_copy.is_valid() &&
	  T1_2_copy.number_of_vertices() == T1_2.number_of_vertices() );
  std::ifstream if1_5("T15.triangulation"); CGAL::set_ascii_mode(if1_5); 
  Triangul T1_5_copy; if1_5 >> T1_5_copy;
  assert( T1_5_copy.is_valid() &&
	  T1_5_copy.number_of_vertices() == T1_5.number_of_vertices() );
  std::ifstream if1_6("T16.triangulation"); CGAL::set_ascii_mode(if1_6);
  Triangul T1_6_copy; if1_6 >> T1_6_copy;
  assert( T1_6_copy.is_valid() &&
	  T1_6_copy.number_of_vertices() == T1_6.number_of_vertices() );
  std::ifstream if2_1("T21.triangulation"); CGAL::set_ascii_mode(if2_1);
  Triangul T2_1_copy; if2_1 >> T2_1_copy;
  assert( T2_1_copy.is_valid() &&
	  T2_1_copy.number_of_vertices() == T2_1.number_of_vertices() );
  std::ifstream if2_3("T23.triangulation"); CGAL::set_ascii_mode(if2_3);
  Triangul T2_3_copy; if2_3 >> T2_3_copy;
  assert( T2_3_copy.is_valid() &&
	  T2_3_copy.number_of_vertices() == T2_3.number_of_vertices() );
  std::ifstream if2_5("T25.triangulation"); CGAL::set_ascii_mode(if2_5); 
  Triangul T2_5_copy; if2_5 >> T2_5_copy;
  assert( T2_5_copy.is_valid() &&
	  T2_5_copy.number_of_vertices() == T2_5.number_of_vertices() );
  std::ifstream if2_6("T26.triangulation"); CGAL::set_ascii_mode(if2_6);
  Triangul T2_6_copy; if2_6 >> T2_6_copy;
  assert( T2_6_copy.is_valid() &&
	  T2_6_copy.number_of_vertices() == T2_6.number_of_vertices() );

  

  /**********************/
  /***** REMOVALS *******/ 
  std::cout << "    removals" << std::endl;

  // test remove_first()
  T0_1.remove_first(T0_1.finite_vertex());
  assert( T0_1.number_of_vertices() == 0 );

  // test remove_second()
  T1_6.remove_second(T1_6.finite_vertex());
  assert( T1_6.is_valid());
  assert( T1_6.number_of_vertices() == 1 );

  // remove from 1-dimensional triangulations
  T1_2.remove(v1_2_1);
  assert( T1_2.is_valid());
  T1_2.remove(v1_2_2);
  assert( T1_2.number_of_vertices() == 0 );

  T1_5.remove(v1_5_1);
  assert( T1_5.is_valid());
  T1_5.remove(v1_5_2);
  T1_5.remove(v1_5_3);
  T1_5.remove(v1_5_8);
  T1_5.remove(v1_5_9);
  assert( T1_5.number_of_vertices() == 0 );

 // remove from 2-dimensional triangulations
  T2_1.remove(v2_1_10); assert( T2_1.is_valid());
  T2_1.remove(v2_1_9);
  T2_1.remove(v2_1_8);
  T2_1.remove(v2_1_7);
  T2_1.remove(v2_1_6);
  T2_1.remove(v2_1_5);
  T2_1.remove(v2_1_4);
  T2_1.remove(v2_1_3); assert(T2_1.is_valid());
  T2_1.remove(v2_1_0);
  assert(T2_1.is_valid());
  T2_1.remove(v2_1_2); assert(T2_1.is_valid());
  T2_1.remove(v2_1_1); assert(T2_1.is_valid());
  assert( T2_1.number_of_vertices() == 0 );

   T2_3.remove(v2_3_0); assert(T2_3.is_valid());
  T2_3.remove(v2_3_1); assert(T2_3.is_valid());
  T2_3.remove(v2_3_9); assert(T2_3.is_valid());
  T2_3.remove(v2_3_8); assert(T2_3.is_valid()); 
  T2_3.remove(v2_3_5); assert(T2_3.is_valid()); 
  T2_3.remove(v2_3_3); assert(T2_3.is_valid());
  T2_3.remove(v2_3_4); assert(T2_3.is_valid());
  T2_3.remove(v2_3_2); assert(T2_3.is_valid());
  T2_3.remove(v2_3_6); assert(T2_3.is_valid());
  T2_3.remove(v2_3_7); assert(T2_3.is_valid()); 
  T2_3.remove(v2_3_10); assert(T2_3.is_valid());
  assert( T2_3.number_of_vertices() == 0 );

  typename Triangul::size_type i;
  for (i=T2_4.number_of_vertices(); i>0; i--)
    T2_4.remove(T2_4.finite_vertex());
  assert( T2_4.number_of_vertices() == 0 );
  
  T2_5.clear();
  assert( T2_5.number_of_vertices() == 0 );

  for (i=T2_6.number_of_vertices(); i>0; i--) {
    assert(T2_6.is_valid());
    T2_6.remove(T2_6.finite_vertex());
  }
  assert( T2_6.number_of_vertices() == 0 );
  
  for (i=T2_7.number_of_vertices(); i>0; i--)
    T2_7.remove(T2_7.finite_vertex());
  assert( T2_7.number_of_vertices() == 0 );

  //test access to tds
  assert (T2_1.tds().is_valid());
  assert (T1_5.tds().is_valid());

  // test destructors and return
  std::cout << "    test destructors and return" << std::endl;
}
