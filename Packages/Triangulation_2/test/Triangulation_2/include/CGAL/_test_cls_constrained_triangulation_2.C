// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL
// release
// of the Computational Geometry Algorithms Library (CGAL). It is
// not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        : 
// file          : 
// revision      : 
// revision_date : 
// author(s)     : Francois Rebufat
// (Francois.Rebufat@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <list.h>
#include <CGAL/_test_cls_triangulation_2.C>

template <class Triangulation>
void 
CGAL__test_cls_constrained_triangulation(const Triangulation &)
{
  typedef Triangulation                       Cls;

  // We assume the traits class has been tested already
  // actually, any traits is good if it has been tested
  typedef typename Cls::Geom_traits          Gt;


  typedef typename Cls::Point                Point;
  typedef typename Cls::Segment              Segment;
  typedef typename Cls::Triangle             Triangle;




  typedef typename Cls::Vertex_handle         Vertex_handle;
  typedef typename Cls::Face_handle           Face_handle;
  typedef pair<Face_handle,int>               Edge;

  typedef typename Cls::Locate_type           Locate_type;

  typedef pair<Point,Point>                   Constraint ;
  typedef list<Constraint>                    list_constraints;
  typedef typename list_constraints::iterator list_iterator;


  // Constructors
  cout << "    constructors    " << endl;
  list_constraints l;

  // Empty triangulation (0-dimensional)
  cout << "    0-Dim " <<endl;
  Cls T0_1;
  assert( T0_1.dimension() == 0 );
  assert( T0_1.number_of_vertices() == 0 );
  assert( T0_1.is_valid() );
  
  l.push_back(Constraint(Point(0,0),Point(0,0)));
  Cls T0_2(l);
  assert( T0_2.dimension() == 0 );
  assert( T0_2.number_of_vertices() == 1 );
  assert( T0_2.is_valid() );
  //l.clear();
  l.erase(l.begin(),l.end());

  // Build dummy triangulations, 1-dimensional

  cout << "    1-Dim "<<endl;
  int m;
  for (m=0; m<20; m++)
      l.push_back(Constraint(Point(3*m, 2*m),Point(3*m,2*m) ));
  Cls T1_1(l);
  assert( T1_1.dimension() == 1 );
  assert( T1_1.number_of_vertices() == 20 );
  assert( T1_1.is_valid() );
  
  //l.clear();
  l.erase(l.begin(),l.end());
  for (m=0; m<19; m++)
      l.push_back(Constraint(Point(3*m, 2*m),Point(3*(m+1),2*(m+1)) ));
  Cls T1_2(l);
  assert( T1_2.dimension() == 1 );
  assert( T1_2.number_of_vertices() == 20);
  assert( T1_2.is_valid() );
  //l.clear();
  l.erase(l.begin(),l.end());

  // Build triangulations, 2-dimensional

  cout << "    2-Dim "<< endl;
  Point lp[5]= {Point(0,0),Point(1,0),Point(0,1),Point(-1,0),Point(0,-1)};
  for (m=1;m<5;m++)
    l.push_back(Constraint(lp[0],lp[m]));
  Cls T2_1(l);
  assert( T2_1.dimension() == 2 );
  assert( T2_1.number_of_vertices() == 5);
  assert( T2_1.is_valid() );
  
  //l.clear();
  l.erase(l.begin(),l.end());

  Point lpt[20] = {
  Point(0,0), Point(1,0), Point(2,0), Point(3,0),Point(4,0),
  Point(4,1), Point(3,1), Point(2,1), Point(1,1),Point(0,1),
  Point(0,2), Point(1,2), Point(2,2), Point(3,2),Point(4,2),
  Point(4,3), Point(3,3), Point(2,3), Point(1,3),Point(0,3)
  };
  for (m=0;m<19;m++) 
    l.push_back(Constraint(lpt[m],lpt[m+1]));
  Cls T2_2(l);
  assert( T2_2.dimension() == 2 );
  assert( T2_2.number_of_vertices() == 20);
  assert( T2_2.is_valid() );

 
  // Build triangulation with iterator
   cout << "    With input iterator" << endl;
   list_iterator first=l.begin();
   list_iterator last=l.end();
   Cls T2_3(first,last);
   assert( T2_3.is_valid() );
   assert(T2_3.number_of_vertices()==T2_2.number_of_vertices());

   // test assignement operator
    Cls Taux = T2_2;
    assert( Taux.dimension() == 2 );
    assert( Taux.number_of_vertices() == 20);
    assert( Taux.is_valid() );

   // Points locations
    // 1-dimensional
   cout << "    point locations 1-dim" << endl;
   Locate_type lt; 
   int            li;
   Face_handle    f;
   f = T1_1.locate(Point(0,0),lt,li); assert( lt == Cls::VERTEX );
   assert( T1_1.geom_traits().compare(f->vertex(li)->point(), Point(0,0)) );
   f = T1_1.locate(Point(9,6),lt,li); assert( lt == Cls::VERTEX );
   assert( T1_1.geom_traits().compare(f->vertex(li)->point(), Point(9,6)) );
   f = T1_1.locate(Point(1.5,1),lt,li); assert( lt == Cls::EDGE );
  assert( (T1_1.geom_traits().compare(f->vertex(f->ccw(li))->point(), Point(0,0))
        && T1_1.geom_traits().compare(f->vertex(f->cw(li))->point(), Point(3,2)))
       || (T1_1.geom_traits().compare(f->vertex(f->ccw(li))->point(), Point(3,2))
        && T1_1.geom_traits().compare(f->vertex(f->cw(li))->point(), Point(0,0))));
  f = T1_1.locate(Point(-3,-2),lt,li); assert( lt == Cls::COLLINEAR_OUTSIDE );
  assert( T1_1.geom_traits().compare(f->vertex(li)->point(), Point(0,0)) );
  f = T1_1.locate(Point(1,1),lt,li); assert( lt == Cls::OUTSIDE );
  li = f->index(T1_1.infinite_vertex());
  assert( CGAL__test_is_to_the_left(T1_1,Point(1,1),f,li) );
  f = T1_1.locate(Point(0,100),lt,li); assert( lt == Cls::OUTSIDE );
  li = f->index(T1_1.infinite_vertex());
  assert( CGAL__test_is_to_the_left(T1_1,Point(0,100),f,li) );

   // 2-dimensional
   cout << "    point locations 2-dim" << endl;
   f = T2_2.locate(Point(0,0),lt,li); assert( lt == Cls::VERTEX );
   assert( T2_2.geom_traits().compare(f->vertex(li)->point(), Point(0,0)) );
   f = T2_2.locate(Point(3,2),lt,li); assert( lt == Cls::VERTEX );
   assert( T2_2.geom_traits().compare(f->vertex(li)->point(), Point(3,2)) );

   f = T2_2.locate(lpt[0],lt,li); assert( lt == Cls::VERTEX );
   assert( T2_2.geom_traits().compare(f->vertex(li)->point(), lpt[0]) );
   f = T2_2.locate(lpt[2],lt,li); assert( lt == Cls::VERTEX );
   assert( T2_2.geom_traits().compare(f->vertex(li)->point(), lpt[2]) );
   f = T2_2.locate(lpt[4],lt,li); assert( lt == Cls::VERTEX );
   assert( T2_2.geom_traits().compare(f->vertex(li)->point(), lpt[4]) );
   f = T2_2.locate(lpt[6],lt,li); assert( lt == Cls::VERTEX );
   assert( T2_2.geom_traits().compare(f->vertex(li)->point(), lpt[6]) );
   f = T2_2.locate(lpt[8],lt,li); assert( lt == Cls::VERTEX );
   assert( T2_2.geom_traits().compare(f->vertex(li)->point(), lpt[8]) );
   f = T2_2.locate(lpt[15],lt,li); assert( lt == Cls::VERTEX );
   assert( T2_2.geom_traits().compare(f->vertex(li)->point(), lpt[15]) );
   f = T2_2.locate(lpt[19],lt,li); assert( lt == Cls::VERTEX );
   assert( T2_2.geom_traits().compare(f->vertex(li)->point(), lpt[19]) );

   f = T2_2.locate(Point(0.5,0),lt,li); assert( lt == Cls::EDGE );
   assert( (T2_2.geom_traits().compare(f->vertex(f->ccw(li))->point(),lpt[1])
        && T2_2.geom_traits().compare(f->vertex(f->cw(li))->point(),lpt[0]))
       || (T2_2.geom_traits().compare(f->vertex(f->ccw(li))->point(), lpt[0])
        && T2_2.geom_traits().compare(f->vertex(f->cw(li))->point(),lpt[1])));
   f = T2_2.locate(Point(10,10),lt,li); assert( lt == Cls::OUTSIDE );
   li = f->index(T2_2.infinite_vertex());
   assert( CGAL__test_is_to_the_left(T2_2,Point(10,10),f,li) );
   f = T2_2.locate(Point(-1,3),lt,li); assert( lt == Cls::OUTSIDE );
   li = f->index(T2_2.infinite_vertex());
   assert( CGAL__test_is_to_the_left(T2_2,Point(-1,3),f,li) );



   /*************************/
  /******* Iterators *******/
   cout << "    iterators" << endl;
   CGAL__test_iterators(T1_1);
   CGAL__test_iterators(T1_2);
   CGAL__test_iterators(T2_1);
   CGAL__test_iterators(T2_2);

   /***************************/
  /******* Circulators *******/
   cout << "    circulators" << endl;
   CGAL__test_circulators(T1_1);
   CGAL__test_circulators(T1_2);
   CGAL__test_circulators(T2_1);
   CGAL__test_circulators(T2_2);


// Line_face_circulator
  cout << "    line face circulator  " << endl;
  typedef typename Cls::Line_face_circulator LFC;

  LFC fc= T2_2.line_walk(Point(-1,-1),Point(10,10));
  assert(fc.ptr()!=NULL);
  assert(!fc.is_empty());
  LFC fc2=fc;
  assert(fc==fc2);
  fc++;
  fc--;
  ++fc;
  --fc;
  Point pp(0.5,0.2);
  f = T2_2.locate(pp,lt,li);
  assert(lt==Cls::FACE);
  fc= T2_2.line_walk(pp,Point(10,10),f);
  fc2=fc;
  assert(fc==fc2);
  fc++;
  fc--;
  ++fc;
  --fc;
  int n=0;
  do {fc2++ ; n = n+1;} while (fc2 != fc);
  // testing with dummy triangulations
  assert(T2_1.number_of_vertices()==5);
  assert(T2_1.is_valid());
  fc= T2_1.line_walk(Point(-5,-5),Point(5,5));
  fc2=fc;
  n=0;
  assert(fc==fc2);
  do {fc2++ ; n = n+1;} while (fc2 != fc);
  assert(n==4);
  fc= T2_1.line_walk(Point(0,-1),Point(1,0));
  fc2=fc;
  n=0;
  assert(fc==fc2);
  do {fc2++ ; n = n+1;} while (fc2 != fc);
  assert(n==3);
  fc= T2_1.line_walk(Point(-1,0),Point(0,1));
  assert(fc.ptr()==NULL);
  assert(fc.is_empty());

   /*****************************/
  /******** Miscellaneaous *****/
  cout << "    misc." << endl;
  assert( T0_1.ccw(0) == 1 );
  assert( T0_1.ccw(1) == 2 );
  assert( T0_1.ccw(2) == 0 );
  assert( T0_1.cw(0) == 2 );
  assert( T0_1.cw(1) == 0 );
  assert( T0_1.cw(2) == 1 );


  f = T2_1.locate(Point(0.2,0.5),lt,li); 
  Triangle t = T2_1.triangle(f); assert( &t == &t );
  Segment  s = T2_1.segment(f,li); assert( &s == &s );
  s = T2_1.segment(Edge(f,li)); assert( &s == &s );
  s = T2_1.segment(T2_1.finite_vertex()->incident_edges()); assert( &s == &s );
  s = T2_1.segment(T2_1.edges_begin()); assert( &s == &s );


   /********************/
  /******** I/O *******/
  cout << "    output to a file" << endl;

  ofstream of0_1("T01.triangulation", ios::out);
  of0_1 << T0_1; of0_1.close();

  ofstream of0_2("T02.triangulation");
  of0_2 << T0_2; of0_2.close();

  ofstream of1_1("T11.triangulation");
  of1_1 << T1_1; of1_1.close();

  ofstream of1_2("T12.triangulation");
  of1_2 << T1_2; of1_2.close();

  ofstream of2_1("T21.triangulation");
  of2_1 << T2_1; of2_1.close();

  ofstream of2_2("T22.triangulation");
  of2_2 << T2_2; of2_2.close();

  cout << "    input from a file" << endl;
  Cls T0_1_copy; ifstream if0_1("T01.triangulation");
  if0_1 >> T0_1_copy;

  Cls T0_2_copy; ifstream if0_2("T02.triangulation");
  if0_2 >> T0_2_copy;

  Cls T1_1_copy; ifstream if1_1("T11.triangulation");
  if1_1 >> T1_1_copy;

  Cls T1_2_copy; ifstream if1_2("T12.triangulation");
  if1_2 >> T1_2_copy;

  Cls T2_1_copy; ifstream if2_1("T21.triangulation");
  if2_1 >> T2_1_copy;

  Cls T2_2_copy; ifstream if2_2("T22.triangulation");
  if2_2 >> T2_2_copy;

   // Actually need overloading functions insert and remove.
   // cout << " Insert and remove "<< endl;
//   // Insert and remove vertex
//    Point p1(0.5,0.5);
//    T4.insert(p1);
//    cout << "IV1"<<endl;
//    assert( T4.is_valid() );
//    Locate_type lt;
//    int li ;
//    Face_handle f=T4.locate(p1,lt,li);
//    Vertex_handle v=f->vertex(li);
//    T4.remove(v);
//    cout << "IV2"<<endl;
//    assert( T4.is_valid() );
//    assert(T5.number_of_vertices()==T4.number_of_vertices());

//   // Insert on edge
//    Point p2(1.5,0);
//    T4.insert(p2);
//   cout << "IV3"<<endl;
//    assert( T4.is_valid() );

//    f=T4.locate(p2,lt,li);
//    assert(f->is_constrained(li));
    
//    v=f->vertex(li);
  
//    T4.remove(v);
//    cout << "IV4"<<endl;
//    assert( T4.is_valid() );
//    assert(T5.number_of_vertices()==T4.number_of_vertices());
  
//   // test set_constraint
//    p1=Point(0.1,0.1);
//    f=T4.locate(p1);
//    f->set_constraint(1,false);
//    assert(f->is_constrained(1)==false);
//    f->set_constraint(1,true);
//    assert(f->is_constrained(1));

   
}
