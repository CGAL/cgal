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
// file          : include/CGAL/_test_cls_tds_2.C
// revision      : 
// revision_date : 
// author(s)     : Mariette Yvinec (Mariette.Yvinec@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <cassert>
#include <fstream>

#include <CGAL/_test_cls_tds_vertex.C>
#include <CGAL/_test_cls_tds_face.C>



template <class Tds, class Gt>
void
_test_cls_tds_2( const Tds &, const Gt &)
{
  
  typedef typename Tds::Vertex            Vertex;
  typedef typename Tds::Face              Face;
  typedef typename Tds::Edge              Edge;

  typedef typename Tds::Vertex_iterator   Vertex_iterator;
  typedef typename Tds::Face_iterator     Face_iterator;
  typedef typename Tds::Edge_iterator     Edge_iterator;
  
  typedef typename Tds::Vertex_circulator Vertex_circulator;
  typedef typename Tds::Face_circulator   Face_circulator;
  typedef typename Tds::Edge_circulator   Edge_circulator;

  typedef typename Gt::Point              Point;

  // Test subclasses
  CGAL::_test_cls_tds_vertex( Vertex(), Gt() );
  CGAL::_test_cls_tds_face( Face(), Gt() );

  // Test constructors
  cout << "    constructors" << endl;
  Tds tds0;
  Tds tds1;
  Tds tds2(tds1); 
  Tds tds3 = tds1;
  Tds tds4 ;
  tds4.swap(tds1);
  tds1.is_valid();
  tds1.clear();
  tds1.is_valid();
  // (), = and  swap to be tested later again with non trivial tds 

//   // Setting functions - access functions
//   cout << "    setting functions  access functions" << endl;
//   assert(tds1.dimension() == -1);
//   assert(tds1.number_of_vertices() == 0 );
  
//   tds1.set_number_of_vertices(1);
//   assert( tds1.number_of_vertices() == 1 );
//   tds1.set_dimension(0);
//   assert(tds1.dimension() == 0);

//   Vertex* vt1 = new Vertex; 
//   Face* f1 = new Face(vt1,NULL,NULL);
//   Vertex *vt2 = new Vertex; 
//   Face *f2  = new Face(vt2,NULL,NULL);
//   vt2->set_face(f2);
//   tds1.set_infinite_vertex(vt2);
 
//   // Finite and infinite vertices and faces
//   cout << "    finite/infinite faces and vertices" << endl;
//   assert( !tds1.is_infinite(vt1) );
//   assert( tds1.is_infinite(vt2) );
//   assert( !tds1.is_infinite(f1) );
//   assert( tds1.is_infinite(f2) );
//   assert( tds1.is_infinite(f2,1) );

//   //assert( tds1.infinite_face() == f2 );
//   assert( tds1.infinite_vertex() == vt2 );
//   //assert( tds1.finite_vertex() == vt1 );

   // make tds1 valid in order to allow clear() to work
//   tds1.set_number_of_vertices(0);
//   tds1.set_dimension(-1);
  //  tds1.set_infinite_vertex(NULL);
  //  assert( tds1.is_valid() );

 
  // misc.
  cout << "    miscellaneous" << endl;
  assert( tds1.ccw(0) == 1 );
  assert( tds1.ccw(1) == 2 );
  assert( tds1.ccw(2) == 0 );
  assert( tds1.cw(0) == 2 );
  assert( tds1.cw(1) == 0 );
  assert( tds1.cw(2) == 1 );
 
  // test insert, remove  and  flip
  // tds1 , tds2 0 dim 
  // tds3 1 dim
  // tds4 2dim

  cout << "    insert and flip" << endl;
  Vertex* w1 = tds1.insert_first();
  assert(tds1.dimension()== -1); 
  assert(tds1.number_of_vertices() == 1);
  assert(tds1.is_valid() );

  tds2.insert_first();
  Vertex* v2 = tds2.insert_second();
  assert(tds2.dimension()== 0); 
  assert(tds2.number_of_vertices() == 2);
  assert(tds2.is_valid() );

  Vertex* w3 = tds3.insert_first();
  tds3.insert_second();
  // test insert_dim_up, remove _dim_down  from dimsension 0 to 1
  Vertex * v3 = tds3.insert_dim_up(w3,false);
  assert(tds3.is_valid());
  assert(tds3.number_of_vertices() == 3);
  assert(tds3.dimension() == 1);
  tds3.remove_dim_down(v3);
  assert(tds3.is_valid());
  v3 = tds3.insert_dim_up(w3,true);
  assert(tds3.is_valid());
  // test insert_in_egde dim==1
  tds3.insert_in_edge(v3->face(), 2);
  assert(tds3.dimension()== 1);
  assert(tds3.number_of_vertices() == 4);
  assert(tds3.is_valid() );

 
  Vertex* w4 = tds4.insert_first();
  Vertex* v4_1 = tds4.insert_second();
  Vertex* v4_2 = tds4.insert_dim_up(w4,true);
   // from now on, coordinates have to be introduced for
  // the ierators in is_valid() to work
  Point p1(0,0);
  Point p2(2,0);
  Point p3(1,2);
  Point p3bis(1,-2);
  Point p4(1,1);
  Point p5(1,0);
  v4_1->set_point(p1);
  v4_2->set_point(p2);
  //test insert_dim_up, remove _dim_down  from dimsension 1 to 2
  Vertex* v4_3 = tds4.insert_dim_up(w4,false);
  v4_3->set_point(p3bis);
  assert(tds4.dimension()== 2);
  assert(tds4.number_of_vertices() == 4);
  assert(tds4.is_valid() );
  tds4.remove_dim_down(v4_3);
  assert(tds4.is_valid() );
  v4_3 = tds4.insert_dim_up(w4,true);
  v4_3->set_point(p3);
  assert(tds4.is_valid() );
  // test insert-in-face, insert_in_egde dim==2
  // Find the face  v4_1 v4_2 v4_3 for insertion
  Face_circulator fc= v4_1->incident_faces();
  while( ! (fc->has_vertex(v4_2) && fc->has_vertex(v4_3)) ) fc++;
  Vertex* v4_4 = tds4.insert_in_face(&( *fc));
  v4_4->set_point(p4);
  assert(tds4.is_valid() );
  // Find the edge v4_1v4_2 for insertion
  fc= v4_1->incident_faces();
  int ic;
  while(! (fc->has_vertex(v4_2, ic ) && ic == fc->ccw(fc->index(v4_1)))) 
    fc++;
  Vertex* v4_5 = tds4.insert_in_edge(&(*fc), ic);
  v4_5->set_point(p5);
  assert(tds4.is_valid() );
  assert(tds4.dimension()== 2);
  assert(tds4.number_of_vertices() == 6);
  assert(tds4.is_valid() );
  //flip
  tds4.flip(v4_1->face(),1);
  assert(tds4.is_valid() );

  //remove_degree_3 remove_1D
  Vertex* u4 = tds4.insert_in_face(v4_1->face());
  tds4.remove_degree_3(u4);
  assert(tds4.is_valid() );
  
  Vertex* u3 = tds3.insert_in_edge(v3->face(),2);
  tds3.remove_1D(u3);
  assert(tds3.is_valid() );
  
  //remove_second, remove first
  tds2.remove_second(v2);
  assert(tds2.is_valid() && tds2.number_of_vertices()==1);
  v2 = tds2.insert_second();
  tds1.remove_first(w1);
  assert(tds1.is_valid()&& tds1.number_of_vertices()==0); 
  w1 = tds1.insert_first();

  //access
  cout << "    test access" << endl;
  assert(tds0.dimension() == -1     && tds0.number_of_vertices() == 0 &&
	 tds0.number_of_faces()== 0 && tds0.number_of_edges()    == 0 &&
	 tds0.number_of_full_dim_faces() == 0);
  assert(tds1.dimension() == -1     && tds1.number_of_vertices() == 1 &&
	 tds1.number_of_faces()== 0 && tds1.number_of_edges()    == 0 &&
	 tds1.number_of_full_dim_faces() == 0);
  assert(tds2.dimension() == 0      && tds2.number_of_vertices() == 2 &&
	 tds2.number_of_faces()== 0 && tds2.number_of_edges()    == 0 &&
	 tds2.number_of_full_dim_faces() == 2);
  assert(tds3.dimension() == 1      && tds3.number_of_vertices() == 4 &&
	 tds3.number_of_faces()== 0 && tds3.number_of_edges()    == 4 &&
	 tds3.number_of_full_dim_faces() == 4);
  assert(tds4.dimension() == 2      && tds4.number_of_vertices() == 6 &&
	 tds4.number_of_faces()== 8 && tds4.number_of_edges()    == 12 &&
	 tds4.number_of_full_dim_faces() == 8);

 
   //clear(), swap() and copy_constructor
  cout << "    clear, swap, copy_constructor and assign " << endl;
  Tds tds1b(tds1);
  assert(tds1b.is_valid());
  Tds tds1c = tds1b;
  assert(tds1c.is_valid());
  Tds tds1d;
  tds1d.swap(tds1c);
  assert(tds1d.is_valid() & tds1d.number_of_vertices()==1);
  tds1d.clear();

  Tds tds2b(tds2);
  assert(tds2b.is_valid());
  Tds tds2c = tds2b;
  assert(tds2c.is_valid());
  Tds tds2d;
  tds2d.swap(tds2c);
  assert(tds2d.is_valid() & tds2d.number_of_vertices()==2);
  tds2d.clear();

  Tds tds3b(tds3);
  assert(tds3b.is_valid());
  Tds tds3c = tds3b;
  assert(tds3c.is_valid());
  Tds tds3d;
  tds3d.swap(tds3c);
  assert(tds3d.is_valid() & tds3d.number_of_vertices()==4);
  tds3d.clear();

  Tds tds4b(tds4);
  assert(tds4b.is_valid());
  Tds tds4c = tds4b;
  assert(tds4c.is_valid());
  Tds tds4d;
  tds4d.swap(tds4c);
  assert(tds4d.is_valid() && tds4d.number_of_vertices()==6);
  tds4d.clear();
	     

   //iterators are tested by is_valid()
  //test circulators and v->degree()
  cout << "    circulators" <<endl;
  _test_tds_circulators(tds0);
  _test_tds_circulators(tds1);
  _test_tds_circulators(tds2);
  _test_tds_circulators(tds3);
  _test_tds_circulators(tds4);

  //TODO
  //test input, output
  cout << "    output to a file" << endl;
  ofstream of0("file_tds0");
  CGAL::set_ascii_mode(of0); 
  of0 << tds0 ; 
  of0.close();
  ofstream of1("file_tds1");
  CGAL::set_ascii_mode(of1); 
  of1 << tds1 ; 
  of1.close();
  ofstream of2("file_tds2");
  CGAL::set_ascii_mode(of2); 
  of2 << tds2 ; 
  of2.close();
  ofstream of3("file_tds3");
  CGAL::set_ascii_mode(of3); 
  of3 << tds3 ; 
  of3.close();
  ofstream of4("file_tds4");
  CGAL::set_ascii_mode(of4); 
  of4 << tds4 ; 
  of4.close();

  cout << "    input from a file" << endl;
  ifstream if0("file_tds0"); CGAL::set_ascii_mode(if0);
  Tds tds0e; if0 >> tds0e ;   assert( tds0e.is_valid());
  ifstream if1("file_tds1"); CGAL::set_ascii_mode(if1);
  Tds tds1e; if1 >> tds1e;  assert( tds1e.is_valid());
  ifstream if2("file_tds2"); CGAL::set_ascii_mode(if2);
  Tds tds2e; if2 >> tds2e ;   assert( tds2e.is_valid());
  ifstream if3("file_tds3"); CGAL::set_ascii_mode(if3);
  Tds tds3e; if3 >> tds3e ;   assert( tds3e.is_valid());
  ifstream if4("file_tds4"); CGAL::set_ascii_mode(if4);
  Tds tds4e; if4 >> tds4e ;   assert( tds4e.is_valid());

  // test destructor and return
  cout << "    destructors and return" << endl;
}

template< class Tds>
void
_test_tds_circulators( const Tds&  tds)
{
  typedef typename Tds::Vertex_iterator   Vertex_iterator;
  typedef typename Tds::Face_iterator     Face_iterator;
  typedef typename Tds::Edge_iterator     Edge_iterator;
  
  typedef typename Tds::Vertex_circulator Vertex_circulator;
  typedef typename Tds::Face_circulator   Face_circulator;
  typedef typename Tds::Edge_circulator   Edge_circulator;

  int countf = 0;
  int counte = 0;
  int countv = 0;
  int countvv = 0;

  for( Vertex_iterator vit = tds.vertices_begin();
	               vit != tds.vertices_end(); vit++) {
    
    Face_circulator fc = vit->incident_faces(), fdone(fc);
    if (! fc.is_empty()) {
      do {
	countf +=1;
      } while (++fc != fdone);
    }

    Edge_circulator ec = vit->incident_edges(), edone(ec);
    if (! ec.is_empty()) {
      do {
	counte +=1;
      } while (++ec != edone);
    }

    countvv = 0;
    Vertex_circulator vc = vit->incident_vertices(), vdone(vc);
    if (! vc.is_empty()) {
      do {
	countv +=1;
	countvv +=1;
      } while (++vc != vdone);
    }

    assert( vit->degree() == countvv);
  }	       
					     
  assert( countf == 3 * tds.number_of_faces());
  assert( counte == 2 * tds.number_of_edges());
  assert( countv == counte);
	  
}
