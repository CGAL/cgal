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
// file          : test/Triangulation/include/CGAL/_test_cls_tds_2.C
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
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

  typedef typename Gt::Point_2              Point;

  // Test subclasses
  CGAL::_test_cls_tds_vertex( Vertex(), Gt() );
  CGAL::_test_cls_tds_face( Face(), Gt() );

  // Test constructors
  std::cout << "    constructors" << std::endl;
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


 
  // misc.
  std::cout << "    miscellaneous" << std::endl;
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

  std::cout << "    insert and flip" << std::endl;
  Vertex* w1 = tds1.insert_first();
  assert(tds1.dimension()== -1); 
  assert(tds1.number_of_vertices() == 1);
  assert(tds1.is_valid() );

  Vertex* w2 = tds2.insert_first();
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

  //  make_hole, star_hole
  u4 = tds4.insert_in_face(v4_1->face());
  typename Tds::List_edges hole;
  tds4.make_hole(u4, hole);
  u4 = tds4.star_hole(hole);
  tds4.remove_degree_3(u4);

 //access
  std::cout << "    test access" << std::endl;
  assert(tds0.dimension() <= -1     && tds0.number_of_vertices() == 0 &&
	 tds0.number_of_faces()== 0 && tds0.number_of_edges()    == 0 &&
	 tds0.number_of_full_dim_faces() == 0);
  assert(tds1.dimension() == -1     && tds1.number_of_vertices() == 1 &&
	 tds1.number_of_faces()== 0 && tds1.number_of_edges()    == 0 &&
	 tds1.number_of_full_dim_faces() >= 0);
  assert(tds2.dimension() == 0      && tds2.number_of_vertices() == 2 &&
	 tds2.number_of_faces()== 0 && tds2.number_of_edges()    == 0 &&
	 tds2.number_of_full_dim_faces() == 2);
  assert(tds3.dimension() == 1      && tds3.number_of_vertices() == 4 &&
	 tds3.number_of_faces()== 0 && tds3.number_of_edges()    == 4 &&
	 tds3.number_of_full_dim_faces() == 4);
  assert(tds4.dimension() == 2      && tds4.number_of_vertices() == 6 &&
	 tds4.number_of_faces()== 8 && tds4.number_of_edges()    == 12 &&
	 tds4.number_of_full_dim_faces() == 8);

 
   //clear(), swap() and copy_constructor and copy
  std::cout << "    clear, swap, assign and copy " << std::endl;
  Tds tds1b(tds1);
  assert(tds1b.is_valid());
  Tds tds1c = tds1b;
  assert(tds1c.is_valid());
  Tds tds1d;
  tds1d.swap(tds1c);
  assert(tds1d.is_valid() & tds1d.number_of_vertices()==1);
  tds1d.clear();
  Tds tds1e;
  Vertex *w1e= tds1e.copy_tds(tds1,w1);
  assert(tds1e.is_valid() && tds1e.is_vertex(w1e));

  Tds tds2b(tds2);
  assert(tds2b.is_valid());
  Tds tds2c = tds2b;
  assert(tds2c.is_valid());
  Tds tds2d;
  tds2d.swap(tds2c);
  assert(tds2d.is_valid() & tds2d.number_of_vertices()==2);
  tds2d.clear();
  Tds tds2e;
  Vertex *w2e= tds2e.copy_tds(tds2,w2);
  assert(tds2e.is_valid() && tds2e.is_vertex(w2e));

  Tds tds3b(tds3);
  assert(tds3b.is_valid());
  Tds tds3c = tds3b;
  assert(tds3c.is_valid());
  Tds tds3d;
  tds3d.swap(tds3c);
  assert(tds3d.is_valid() & tds3d.number_of_vertices()==4);
  tds3d.clear();
  Tds tds3e;
  Vertex *w3e= tds3e.copy_tds(tds3,w3);
  assert(tds3e.is_valid()&& tds3e.is_vertex(w3e));


  Tds tds4b(tds4);
  assert(tds4b.is_valid());
  Tds tds4c = tds4b;
  assert(tds4c.is_valid());
  Tds tds4d;
  tds4d.swap(tds4c);
  assert(tds4d.is_valid() && tds4d.number_of_vertices()==6);
  tds4d.clear();
  Tds tds4e;
  Vertex *w4e= tds4e.copy_tds(tds4,w4);
  assert(tds4e.is_valid() && tds4e.is_vertex(w4e));

  //test circulators and v->degree()
  std::cout << "    circulators" <<std::endl;
  _test_tds_circulators(tds0);
  _test_tds_circulators(tds1);
  _test_tds_circulators(tds2);
  _test_tds_circulators(tds3);
  _test_tds_circulators(tds4);

  _test_tds_iterators(tds0);
  _test_tds_iterators(tds1);
  _test_tds_iterators(tds2);
  _test_tds_iterators(tds3);
  _test_tds_iterators(tds4);

  //is_vertex
  std::cout << "    is_vertex" << std::endl;
  assert (tds4.is_vertex(v4_5));
  assert (tds3.is_vertex(v3));
  assert (tds2.is_vertex(v2));
  


  //test input, output

  // first set_point to vertices which have not been set
  //for clean input output of points
  w1->set_point(p1);
  w2->set_point(p1);
  v2->set_point(p1);
  w3->set_point(p1);
  v3->set_point(p1);
  w4->set_point(p1); 
  Vertex_iterator vit = tds3.vertices_begin();
  for ( ; vit != tds3.vertices_end(); vit++) {
    vit->set_point(p1);
  }
  std::cout << "    output to a file" << std::endl;
  std::ofstream of0("file_tds0");
  CGAL::set_ascii_mode(of0); 
  of0 << tds0 ; 
  of0.close();
  std::ofstream of1("file_tds1");
  CGAL::set_ascii_mode(of1); 
  of1 << tds1 ; 
  of1.close();
  std::ofstream of2("file_tds2");
  CGAL::set_ascii_mode(of2); 
  of2 << tds2 ; 
  of2.close();
  std::ofstream of3("file_tds3");
  CGAL::set_ascii_mode(of3); 
  of3 << tds3 ; 
  of3.close();
  std::ofstream  of4("file_tds4");
  CGAL::set_ascii_mode(of4); 
  of4 << tds4 ; 
  of4.close();

  std::cout << "    input from a file" << std::endl;
  std::ifstream if0("file_tds0"); CGAL::set_ascii_mode(if0);
  Tds tds0f; if0 >> tds0f ;   assert( tds0f.is_valid());
  std::ifstream if1("file_tds1"); CGAL::set_ascii_mode(if1);
  Tds tds1f; if1 >> tds1f;  
  assert( tds1f.is_valid());
  std::ifstream if2("file_tds2"); 
  CGAL::set_ascii_mode(if2);
  Tds tds2f; if2 >> tds2f ;   
  assert( tds2f.is_valid());
  std::ifstream if3("file_tds3"); 
  CGAL::set_ascii_mode(if3);
  Tds tds3f; if3 >> tds3f ;   
  assert( tds3f.is_valid());
  std::ifstream if4("file_tds4"); 
  CGAL::set_ascii_mode(if4);
  Tds tds4f; if4 >> tds4f ;   
  assert( tds4f.is_valid());

  // test destructor and return
  std::cout << "    destructors and return" << std::endl;
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

template< class Tds>
void
_test_tds_iterators( const Tds&  tds)
{
  typedef typename Tds::Edge     Edge;
  typedef typename Tds::Vertex   Vertex;
  typedef typename Tds::Face     Face;

  typedef typename Tds::Vertex_iterator   Vertex_iterator;
  typedef typename Tds::Face_iterator     Face_iterator;
  typedef typename Tds::Edge_iterator     Edge_iterator;
  typedef typename Tds::Iterator_base     Iterator_base;

  int nv, ne, nf, nb;
  nv = ne = nf = nb = 0;

  Vertex v;
  Face f;
  Edge e;

  for (Iterator_base itbp = tds.iterator_base_begin();
       itbp != tds.iterator_base_end();
       ++itbp) {
   f = *itbp;
   nb += 1;
  }
  assert(nb == tds.number_of_full_dim_faces());
  for (Iterator_base itbm = tds.iterator_base_end();
       itbm != tds.iterator_base_begin();
       --itbm) {
    nb -= 1;
  }
  assert(nb == 0);

  for (Vertex_iterator vitp = tds.vertices_begin();
       vitp != tds.vertices_end();
       ++vitp) {
   v = *vitp;
   nv += 1;
  }
  assert(nv == tds.number_of_vertices());
  for (Vertex_iterator vitm = tds.vertices_end();
       vitm != tds.vertices_begin();
       --vitm) {
    nv -= 1;
  }
  assert(nv == 0);

  
  for (Face_iterator fitp = tds.faces_begin();
       fitp != tds.faces_end();
       ++fitp) {
    f = *fitp;
    nf += 1;
  }
  assert(nf == tds.number_of_faces());
  for (Face_iterator fitm = tds.faces_end();
       fitm != tds.faces_begin();
       --fitm) {
    nf -= 1;
  }
  assert(nf == 0);

  for (Edge_iterator eitp = tds.edges_begin();
       eitp != tds.edges_end();
       ++eitp) {
    e = *eitp;
    ne += 1;
  }
  assert(ne == tds.number_of_edges());
  for (Edge_iterator eitm = tds.edges_end();
       eitm != tds.edges_begin();
       --eitm) {
    ne -= 1;
  }
  assert(ne == 0);

  return;
}
