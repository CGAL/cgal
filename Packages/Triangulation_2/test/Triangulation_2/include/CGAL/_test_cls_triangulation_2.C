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
#include <vector>

#include <CGAL/_test_cls_triangulation_vertex.C>
#include <CGAL/_test_cls_triangulation_face.C>
#include <CGAL/_test_fct_is_infinite.C>
#include <CGAL/_test_triangulation_iterators.C>
#include <CGAL/_test_triangulation_circulators.C>



template < class VIt >
class V2p_adaptor : public VIt {
public:
  typedef typename VIt::Vertex::Point Point;
  V2p_adaptor(const VIt &vit) : VIt(vit) {}
  V2p_adaptor(VIt &vit) : VIt(vit) {}
  const Point& operator*() const { return (VIt::operator*()).point(); }
};


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
_test_cls_triangulation_2( const Triangulation & )
{
  typedef Triangulation                      Cls;

  // We assume the traits class has been tested already
  // actually, any traits is good if it has been tested
  typedef typename Cls::Geom_traits          Gt;

  typedef typename Cls::Point                Point;
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
  
  /***********************/
  /***** SUBCLASSES ******/
  _test_cls_triangulation_vertex(Vertex());
  _test_cls_triangulation_face(Face());
  
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
  // assert( T1.geom_traits() == T2.geom_traits() );
  // assert( T3.geom_traits() == T4.geom_traits() );

  /**************************/
  /******* INSERTIONS *******/
  
  // Tk denotes a k-dimensional triangulation
  // the handle returned when inserting pj into Tk_n is called vk_n_j
  // the asserts at the end of the insert() are to avoid compiler
  // warnings about unused variables
  // we use some of these variables for the remove (see the end)
  // so there is no need to put assert at the end of all of them

  /******* 0-dimensional triangulations ******/
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

  // test insert_first()
  Cls T0_2; 
  // T0_2.insert_first(Vertex(p0).handle());
  // this  statement cause a segmentation fault on Linux
  // when the whole procedure is leaved
  Vertex_handle v0_2_0(new Vertex(p0));
  T0_2.insert_first(v0_2_0);

  assert( T0_2.dimension() == 0 );
  assert( T0_2.number_of_vertices() == 1 );
  assert( T0_2.is_valid() );
  
  /******** 1-dimensional triangulations ******/
  // T1_n denotes a 1-dimensional triangulation with n vertices
  // when there are several, we use T1_n_p
  cout << "    insertions 1-dim" << endl;
  
  Cls T1_2;
  Vertex_handle v1_2_1 = T1_2.insert(p1);
  Vertex_handle v1_2_2 = T1_2.insert(p2);
  assert( T1_2.dimension() == 1 );
  assert( T1_2.number_of_vertices() == 2 );
  // assert( T1_2.number_of_faces() == 0 );
  assert( T1_2.is_valid() );
  
  // p1,p2,p3  [endpoints first]
  Cls T1_3_0;
  Vertex_handle v1_3_0_1 = T1_3_0.insert(p1); assert( !v1_3_0_1.is_null() );
  Vertex_handle v1_3_0_3 = T1_3_0.insert(p3); assert( !v1_3_0_3.is_null() );
  Vertex_handle v1_3_0_2 = T1_3_0.insert(p2); assert( !v1_3_0_2.is_null() );
  assert( T1_3_0.dimension() == 1 );
  assert( T1_3_0.number_of_vertices() == 3 );
  // assert( T1_3_0.number_of_faces() == 0 );
  assert( T1_3_0.is_valid() );
  
  // p1,p2,p3  [middle point first]
  Cls T1_3_1;
  Vertex_handle v1_3_1_1 = T1_3_1.insert(p1); assert( !v1_3_1_1.is_null() );
  Vertex_handle v1_3_1_3 = T1_3_1.insert(p3); assert( !v1_3_1_3.is_null() );
  Vertex_handle v1_3_1_2 = T1_3_1.insert(p2); assert( !v1_3_1_2.is_null() );
  assert( T1_3_1.dimension() == 1 );
  assert( T1_3_1.number_of_vertices() == 3 );
  // assert( T1_3_1.number_of_faces() == 0 );
  // assert( T1_3_0 == T1_3_1 ); // operator== is not defined!
  assert( T1_3_1.is_valid() );

  Cls T1_5;
  Vertex_handle v1_5_1 = T1_5.insert(p1);
  Vertex_handle v1_5_2 = T1_5.insert(p2);
  Vertex_handle v1_5_3 = T1_5.insert(p3);
  Vertex_handle v1_5_8 = T1_5.insert(p8);
  Vertex_handle v1_5_9 = T1_5.insert(p9);
  assert( T1_5.dimension() == 1 );
  assert( T1_5.number_of_vertices() == 5 );
  // assert( T1_5.number_of_faces() == 0 );
  assert( T1_5.is_valid() );

  // test insert_second()
  Cls T1_6 = T0_2; 
  //T1_6.insert_second(Vertex(p3).handle());
  // the following statement cause a segmentation fault on Linux
  // when the whole procedure is leaved
  T1_6.insert_second( Vertex_handle(new Vertex(p3)));
  assert( T1_6.dimension() == 1 );
  assert( T1_6.number_of_vertices() == 2 );
  assert( T1_6.is_valid() ); 
  
  /******** 2-dimensional triangulations ******/ 
  cout << "    insertions 2-dim" << endl;
  
  Cls T2_1;
  Vertex_handle v2_1_0 = T2_1.insert(p0);
  Vertex_handle v2_1_1 = T2_1.insert(p1);
  Vertex_handle v2_1_2 = T2_1.insert(p2);
  Vertex_handle v2_1_3 = T2_1.insert(p3);  // on the edge p1,p2, on the convex hull
  Vertex_handle v2_1_4 = T2_1.insert(p4);  // outside, with two visible collineaar edges
  Vertex_handle v2_1_5 = T2_1.insert(p5); 
  Vertex_handle v2_1_6 = T2_1.insert(p6);  // outside, collinear with p2,p5
  Vertex_handle v2_1_7 = T2_1.insert(p7);  // outside with two visible collinear edges
                      		           // but also collinear with and extending p0,p5 
  Vertex_handle v2_1_8 = T2_1.insert(p8); 
  Vertex_handle v2_1_9 = T2_1.insert(p9);    // inside, on the edge p6,p7
  Vertex_handle v2_1_10 = T2_1.insert(p10);  // inside the face p2,p4,p6
  assert( T2_1.dimension() == 2 );
  assert( T2_1.number_of_vertices() == 11 );
  
  // test is_valid for 2-triangulations
  assert( T2_1.is_valid() );

  // we now test the other insert functions
  // more vicious, we insert all the points on a single line first
  Cls T2_3;
  Locate_type lt;
  Vertex_handle v2_3_1 = T2_3.insert(p1);
  Vertex_handle v2_3_2 = T2_3.insert(p2);
  Vertex_handle v2_3_3 = T2_3.insert(p3, lt);
  assert( lt == Cls::EDGE );
  Vertex_handle v2_3_8 = T2_3.insert(p8, lt);
  assert( lt == Cls::COLLINEAR_OUTSIDE );
  Vertex_handle v2_3_9 = T2_3.insert(p9);
  assert( T2_3.dimension() == 1 );
  Vertex_handle v2_3_4 = T2_3.insert(p4);
  assert( T2_3.dimension() == 2 );
  Vertex_handle v2_3_6 = T2_3.insert(p6, T2_3.faces_begin());
  Vertex_handle v2_3_0 = T2_3.insert(p0, lt, ++T2_3.faces_begin());
  assert( lt == Cls::OUTSIDE );
  Vertex_handle v2_3_5 = T2_3.insert(p5);
  Vertex_handle v2_3_7 = T2_3.insert(p7);
  Vertex_handle v2_3_10 = T2_3.insert(p10, lt, ++(++(T2_3.faces_begin())));
  assert( lt == Cls::FACE );
  assert( T2_3.dimension() == 2 );
  assert( T2_3.number_of_vertices() == 11 );
  // assert( T2_3.number_of_faces() == 13 );
  assert( T2_3.is_valid() );
  
  // make sure inserting on a previous point does not insert it again
  assert( T2_3.insert(p10, lt) == v2_3_10 );
  assert( lt == Cls::VERTEX );
  assert( T2_3.number_of_vertices() == 11 );

  // make sure push_back exists and does the same thing as insert
  assert( T2_3.push_back(p10) == v2_3_10 );
  assert( T2_3.number_of_vertices() == 11 );

  // test generic iterator insert
#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  typedef V2p_adaptor<Vertex_iterator> Point_iterator;
  Cls T2_4; T2_4.insert( Point_iterator(T2_1.vertices_begin()),
                         Point_iterator(T2_1.vertices_end()) );
  assert( T2_4.dimension() == 2 );
  assert( T2_4.number_of_vertices() == 11 );
  // assert( T2_4.number_of_faces() == 13 );
  assert( T2_4.is_valid() );
#endif

  // test list iterator insert
  Cls T2_5;
  assert( T2_5.insert(l.begin(), l.end()) == 10 );
  assert( T2_5.dimension() == 2 );
  assert( T2_5.number_of_vertices() == 10 );
  // assert( T2_5.number_of_faces() == 13 );
  assert( T2_5.is_valid() );

  // test list iterator insert
  Cls T2_6;
  assert( T2_6.insert(v.begin(), v.end()) == 10 );
  assert( T2_6.dimension() == 2 );
  assert( T2_6.number_of_vertices() == 10 );
  // assert( T2_6.number_of_faces() == 13 );
  assert( T2_6.is_valid() );
  
  // test grid insert
  Cls T2_7;
  int m, p;
  for (m=0; m<20; m++)
    for (p=0; p<20; p++)
      T2_7.insert( Point(m*px+p*qx, m*py+p*qy, 1) );
  assert( T2_7.number_of_vertices() == m*p );
  assert( T2_7.is_valid() );

  // test flip
     cout << "    test flip " << endl;
     Cls T2_8=T2_7;
     assert( T2_8.is_valid() );
     Face_handle ff = T2_8.locate(Point(5*px+5*qx,5*py+5*qy));
     assert(!T2_8.is_infinite(ff));
     Face_handle f2 = ff->neighbor(0);
     assert(!T2_8.is_infinite(f2));
     // T2_8.flip(ff,0);
     // Ok precondition violation  (non convex quadrilater)   
     assert( T2_8.is_valid() );
     // try with a valid face (and neighbor)
     T2_8.clear();
     T2_8.insert(Point(0,0,1));
     T2_8.insert(Point(1,0,1));
     T2_8.insert(Point(1,1,1));
     T2_8.insert(Point(0,1,1));
     ff = T2_8.locate(Point(1,1,2));
     assert(!T2_8.is_infinite(ff));
     f2 = ff->neighbor(0);
     assert(!T2_8.is_infinite(f2));
     T2_8.flip(ff,0);
     assert( T2_8.is_valid() );
     

  /****************************/
  /***** CONSTRUCTORS (2) *****/
  cout << "    constructors (2)" << endl;

  // test copy_constructor with non-empty 0-triangulation
  Cls T0_1_1( T0_1 );
  assert( T0_1_1.dimension() == 0 );
  assert( T0_1_1.number_of_vertices() == 1 );
  assert( T0_1_1.is_valid() );

  // test constructor that takes a vertex handle
  // this must(!) be the vertex at infinity
  Cls T0_1_2( T0_1.infinite_vertex() );
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
  Cls T1_5_3( T1_5_2.infinite_vertex() );
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
  Cls T2_1_2( T2_1.infinite_vertex() );
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
#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  _test_fct_is_infinite( T2_4 );
#endif
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
  f = T2_1.locate(p11,lt,li); assert( lt == Cls::EDGE );
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

  // test grid locate
  for (m=0; m<20; m++)
    for (p=0; p<20; p++)
      {
	Point q= Point(m*px+p*qx, m*py+p*qy, 1);
       	f = T2_7.locate(q,lt,li); assert( lt == Cls::VERTEX );
  	assert( T2_7.geom_traits().compare(f->vertex(li)->point(), q) );
      }
  for (m=0; m<20; m+=19)
    for (p=0; p<19; p++)
      {
	Point q= Point(2*m*px+(2*p+1)*qx, 2*m*py+(2*p+1)*qy, 2);
	Point r= Point(m*px+p*qx, m*py+p*qy, 1);
	Point s= Point(m*px+(p+1)*qx, m*py+(p+1)*qy, 1);
       	f = T2_7.locate(q,lt,li); assert( lt == Cls::EDGE );
        assert( (T2_7.geom_traits().compare(f->vertex(f->ccw(li))->point(), r)
              && T2_7.geom_traits().compare(f->vertex(f->cw(li))->point(), s))
             || (T2_7.geom_traits().compare(f->vertex(f->ccw(li))->point(), s)
              && T2_7.geom_traits().compare(f->vertex(f->cw(li))->point(), r)));

      }
  for (m=0; m<19; m++)
    for (p=0; p<19; p++)
      {
	Point q= Point((50*m+1)*px+(50*p+1)*qx, (50*m+1)*py+(50*p+1)*qy, 50);
       	f = T2_7.locate(q,lt,li); assert( lt == Cls::FACE );
	assert( T2_7.oriented_side(f,q) == CGAL::ON_POSITIVE_SIDE );
      }

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
  Point pp(0,1,2); //Point pp(0,0.5);
  f = T2_7.locate(pp,lt,li);
  assert(lt==Cls::FACE);
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
  Cls TT;
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

  // finite/infinite vertex
  // T2_1.set_finite_vertex(v2_1_1);
  // assert( T2_1.is_valid() );
  // T2_1.set_finite_vertex(v2_1_6);
  // assert( T2_1.is_valid() );
  // T2_1.set_infinite_vertex(T2_1.infinite_vertex());
  // assert( T2_1.is_valid() );

  /********************/
  /******** I/O *******/
  cout << "    output to a file" << endl;
  ofstream of0_0("T00.triangulation", ios::out);
  of0_0 << T0_1; of0_0.close();
  ofstream of0_1("T01.triangulation");
  of0_1 << T0_1; of0_1.close();
  ofstream of1_2("T12.triangulation");
  of1_2 << T1_2; of1_2.close();
  ofstream of1_5("T15.triangulation");
  of1_5 << T1_5; of1_5.close();
  ofstream of1_6("T16.triangulation");
  of1_6 << T1_6; of1_6.close();
  ofstream of2_1("T21.triangulation");
  of2_1 << T2_1; of2_1.close();
  ofstream of2_3("T23.triangulation");
  of2_3 << T2_3; of2_3.close();
  ofstream of2_5("T25.triangulation");
  of2_5 << T2_5; of2_5.close();
  ofstream of2_6("T26.triangulation");
  of2_6 << T2_6; of2_6.close();

  cout << "    input from a file" << endl;
  Cls T0_0_copy; ifstream if0_0("T00.triangulation");
  if0_0 >> T0_0_copy;
  // assert( T0_0_copy.number_of_vertices() == T0_0.number_of_vertices() );
  Cls T0_1_copy; ifstream if0_1("T01.triangulation");
  if0_1 >> T0_1_copy;
  // assert( T0_1_copy.number_of_vertices() == T0_1.number_of_vertices() );
  Cls T1_2_copy; ifstream if1_2("T12.triangulation");
  if1_2 >> T1_2_copy;
  // assert( T1_2_copy.number_of_vertices() == T1_2.number_of_vertices() );
  Cls T1_5_copy; ifstream if1_5("T15.triangulation");
  if1_5 >> T1_5_copy;
  // assert( T1_5_copy.number_of_vertices() == T1_5.number_of_vertices() );
  Cls T1_6_copy; ifstream if1_6("T16.triangulation");
  if1_6 >> T1_6_copy;
  // assert( T1_6_copy.number_of_vertices() == T1_6.number_of_vertices() );
  Cls T2_1_copy; ifstream if2_1("T21.triangulation");
  if2_1 >> T2_1_copy;
  // assert( T2_1_copy.number_of_vertices() == T2_1.number_of_vertices() );
  Cls T2_3_copy; ifstream if2_3("T23.triangulation");
  if2_3 >> T2_3_copy;
  // assert( T2_3_copy.number_of_vertices() == T2_3.number_of_vertices() );
  Cls T2_5_copy; ifstream if2_5("T25.triangulation");
  if2_5 >> T2_5_copy;
  // assert( T2_5_copy.number_of_vertices() == T2_5.number_of_vertices() );
  Cls T2_6_copy; ifstream if2_6("T26.triangulation");
  if2_6 >> T2_6_copy;
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
#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  for (i=T2_4.number_of_vertices(); i>0; i--)
    T2_4.remove(T2_4.finite_vertex());
  assert( T2_4.number_of_vertices() == 0 );
#endif
  
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
