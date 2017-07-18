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
// file          : test/Triangulation/include/CGAL/_test_cls_constrained...
// source        : $URL$
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Francois Rebufat, Mariette Yvinec
//
// coordinator   : INRIA Sophia-Antipolis Mariette.Yvinec@sophia.inria.fr
// ============================================================================

#include <list>
#include <CGAL/_test_cls_triangulation_short_2.h>

template <class Triang>
void 
_test_cls_constrained_triangulation(const Triang &)
{
  // typedef Triangulation                       Cls;
  typedef typename Triang::Geom_traits          Gt;

  typedef typename Triang::Point                Point;
  typedef typename Triang::Segment              Segment;
  typedef typename Triang::Triangle             Triangle;

  typedef typename Triang::Vertex_handle        Vertex_handle;
  typedef typename Triang::Face_handle          Face_handle;
  typedef std::pair<Face_handle,int>            Edge;

  typedef typename Triang::Locate_type          Locate_type;
  typedef typename Triang::All_faces_iterator   All_faces_iterator;

  typedef typename Triang::Constrained_edges_iterator Constrained_edges_iterator;

  typedef std::pair<Point,Point>                Constraint ;
  typedef std::list<Edge>                       List_edges;
  typedef std::list<Constraint>                 list_constraints;
  typedef typename list_constraints::iterator   list_iterator;

  CGAL_USE_TYPE(Gt);
  CGAL_USE_TYPE(Segment);
  CGAL_USE_TYPE(Triangle);

  /***********************/
  /***** SUBCLASSE ******/
   _test_cls_triangulation_short_2( Triang() );

  // Constructors
  std::cout << "constructors    " << std::endl;
  list_constraints l;

  // Empty triangulation (0-dimensional)
  std::cout << "    0-Dim " <<std::endl;
  Triang T0_1;
  assert( T0_1.dimension() == -1 );
  assert( T0_1.number_of_vertices() == 0 );
  assert( T0_1.is_valid() );
  
  l.push_back(Constraint(Point(0,0),Point(0,0)));
  Triang T0_2(l);
  assert( T0_2.dimension() == 0 );
  assert( T0_2.number_of_vertices() == 1 );
  assert( T0_2.is_valid() );
  l.erase(l.begin(),l.end());

  // Build triangulations, 1-dimensional
  std::cout << "    1-Dim "<<std::endl;
  int m;
  for (m=0; m<4; m++)
    l.push_back(Constraint(Point(3*m, 2*m),Point(3*m,2*m) ));
  Triang T1_1(l);
  assert( T1_1.dimension() == 1 );
  assert( T1_1.number_of_vertices() == 4 );
  assert( T1_1.is_valid() );
  
  l.erase(l.begin(),l.end());
  for (m=0; m<4; m++)
    l.push_back(Constraint(Point(3*m, 2*m),Point(3*(m+1),2*(m+1)) ));
  Triang T1_2(l);
  assert( T1_2.dimension() == 1 );
  assert( T1_2.number_of_vertices() == 5);
  assert( T1_2.is_valid() );
  l.erase(l.begin(),l.end());

  // Build triangulations, 2-dimensional
  std::cout << "    2-Dim "<< std::endl;
  Point lp[5]= {Point(0,0),Point(1,0),Point(0,1),Point(-1,0),Point(0,-1)};
  for (m=1;m<5;m++)
    l.push_back(Constraint(lp[0],lp[m]));
  Triang T2_1(l);
  assert( T2_1.dimension() == 2 );
  assert( T2_1.number_of_vertices() == 5);
  assert( T2_1.is_valid() );
 
  l.erase(l.begin(),l.end());

  Point lpt[20] = {
  Point(0,0), Point(1,0), Point(2,0), Point(3,0),Point(4,0),
  Point(4,1), Point(3,1), Point(2,1), Point(1,1),Point(0,1),
  Point(0,2), Point(1,2), Point(2,2), Point(3,2),Point(4,2),
  Point(4,3), Point(3,3), Point(2,3), Point(1,3),Point(0,3)
  };
  for (m=0;m<19;m++) 
    l.push_back(Constraint(lpt[m],lpt[m+1]));
  Triang T2_2(l);
  assert( T2_2.dimension() == 2 );
  assert( T2_2.number_of_vertices() == 20);
  assert( T2_2.is_valid() );

 
  // Build triangulation with iterator
   std::cout << "    with input iterator" << std::endl;
   list_iterator first=l.begin();
   list_iterator last=l.end();
   Triang T2_3(first,last);
   assert( T2_3.is_valid() );
   assert(T2_3.number_of_vertices()==T2_2.number_of_vertices());

   // There are too many constraint here to make difference
   // between constrained and Constrained Delaunay
   Triang T2_5;
   for (int j=0; j < 20; j++) T2_5.insert(lpt[j]);
   T2_5.insert( Point(1,0.5), Point(2.5, 3.5));
   T2_5.is_valid();


   // test assignement operator
    Triang Taux = T2_2;
    assert( Taux.dimension() == 2 );
    assert( Taux.number_of_vertices() == 20);
    assert( Taux.is_valid() );


   /********************/
  /******** I/O *******/
  std::cout << "output to a file" << std::endl;

  std::ofstream of0_1("T01.triangulation", std::ios::out);
  CGAL::set_ascii_mode(of0_1);
   of0_1 << T0_1; of0_1.close();

  std::ofstream of0_2("T02.triangulation");
  CGAL::set_ascii_mode(of0_2);
  of0_2 << T0_2; of0_2.close();

  std::ofstream of1_1("T11.triangulation");
  CGAL::set_ascii_mode(of1_1);
  of1_1 << T1_1; of1_1.close();

  std::ofstream of1_2("T12.triangulation");
  CGAL::set_ascii_mode(of1_2);
   of1_2 << T1_2; of1_2.close();
  
  std::ofstream of2_1("T21.triangulation");
  CGAL::set_ascii_mode(of2_1);
  of2_1 << T2_1; of2_1.close();

  std::ofstream of2_2("T22.triangulation");
  CGAL::set_ascii_mode(of2_2);
  of2_2 << T2_2; of2_2.close();

  std::cout << "input from a file" << std::endl;
  std::ifstream if0_1("T01.triangulation"); CGAL::set_ascii_mode(if0_1);
  Triang T0_1_copy;   if0_1 >> T0_1_copy;

  std::ifstream if0_2("T02.triangulation"); CGAL::set_ascii_mode(if0_2);
  Triang T0_2_copy;  if0_2 >> T0_2_copy;

  std::ifstream if1_1("T11.triangulation"); CGAL::set_ascii_mode(if1_1);
  Triang T1_1_copy; if1_1 >> T1_1_copy;

  std::ifstream if1_2("T12.triangulation"); CGAL::set_ascii_mode(if1_2);
   Triang T1_2_copy; if1_2 >> T1_2_copy;

  std::ifstream if2_1("T21.triangulation"); CGAL::set_ascii_mode(if2_1);
  Triang T2_1_copy; if2_1 >> T2_1_copy;

  std::ifstream if2_2("T22.triangulation"); CGAL::set_ascii_mode(if2_2);
  Triang T2_2_copy; if2_2 >> T2_2_copy;

  // test copy of constrained Triangulation
   Triang T2_4(T2_2);
  std::ofstream of2_2_bis("T22.triangulation");
  CGAL::set_ascii_mode(of2_2_bis);
  of2_2_bis << T2_4; of2_2_bis.close();
  All_faces_iterator fit2 = T2_2.all_faces_begin();
  All_faces_iterator fit2_bis = T2_4.all_faces_begin();
  for( ; fit2 != T2_2.all_faces_end(); ++fit2, ++fit2_bis) {
    for(int i=0; i<3 ; i++)  
    assert( fit2->is_constrained(i) ==  fit2_bis->is_constrained(i) );
  }
  
  
  
  // remove_constraint and remove _1 dim
  std::cout << "remove_constrained_edge and remove 1-dim" << std::endl;
  Face_handle fh;
  int ih;
  Vertex_handle vha, vhb;
  Locate_type lt; 
  int li;
  fh  =  T1_2.locate(Point(0,0),lt,li); assert( lt == Triang::VERTEX );
  vha = fh->vertex(li);
  fh  =  T1_2.locate(Point(3,2),lt,li); assert( lt == Triang::VERTEX );
  vhb =  fh->vertex(li);
  assert(T1_2.is_edge(vha,vhb, fh, ih));
  assert(fh->is_constrained(ih));
  T1_2.remove_constrained_edge(fh,ih);
  assert(T1_2.is_valid());
  T1_2.insert(Point(0,0),Point(3,2));
  fh  =  T1_2.locate(Point(3,2),lt,li); assert( lt == Triang::VERTEX );
  vhb =  fh->vertex(li);
  assert(T1_2.are_there_incident_constraints(vhb));
  T1_2.remove_incident_constraints(vhb);
  T1_2.remove(vhb);
  fh  =  T1_2.locate(Point(6,4),lt,li); assert( lt == Triang::VERTEX );
  vha = fh->vertex(li);
  List_edges edges;
  assert(T1_2.are_there_incident_constraints(vha,
					     std::back_inserter(edges)));
  List_edges ic_edges;
  std::back_insert_iterator<List_edges> inserter(ic_edges);
  inserter  = T1_2.incident_constraints(vha, inserter);
  assert(ic_edges.size() == 1); 
  T1_2.remove_incident_constraints(vha);
  inserter = T1_2.incident_constraints(vha, inserter );
  assert(ic_edges.size() == 1); 
  T1_2.remove(vha);
  assert(T1_2.is_valid());
  

   // remove_constraint and remove 2 dim
  std::cout << "remove_constrained_edge and remove 2-dim " << std::endl;
  m=2;
  fh = T2_2.locate(lpt[m], lt,li); assert( lt == Triang::VERTEX );
  vha = fh->vertex(li);
  fh  =  T2_2.locate(lpt[m+1],lt,li); assert( lt == Triang::VERTEX );
  vhb =  fh->vertex(li);
  bool check = T2_2.is_edge(vha,vhb, fh, ih);
  assert(check);
  assert(fh->is_constrained(ih));
  T2_2.remove_constrained_edge(fh,ih);
  T2_2.insert(lpt[m], lpt[m+1]);
  assert(T2_2.is_valid());
  fh  =  T2_2.locate(lpt[m+1],lt,li); assert( lt == Triang::VERTEX );
  vhb =  fh->vertex(li);  
  assert(T2_2.are_there_incident_constraints(vhb));
  T2_2.remove_incident_constraints(vhb);
  T2_2.remove(vhb);

  m=15;
  fh  =  T2_2.locate(lpt[m],lt,li); assert( lt == Triang::VERTEX );
  vha = fh->vertex(li);
  edges.clear();
  ic_edges.clear();
  assert(T2_2.are_there_incident_constraints(vha));
  assert(T2_2.are_there_incident_constraints(vha, 
					     std::back_inserter(edges)));
  ic_edges.clear();
  inserter = std::back_insert_iterator<List_edges>(ic_edges);
  inserter  = T2_2.incident_constraints(vha,inserter);
  assert(ic_edges.size() == 2); 
  T2_2.remove_incident_constraints(vha);
  inserter = T2_2.incident_constraints(vha, inserter);
  assert(ic_edges.size() == 2);
  T2_2.remove(vha);
  assert(T2_2.is_valid());

  // test Constained_edges_iterator
  Triang T2_6;

  T2_6.insert(Point(0,0));
  T2_6.insert(Point(10,0));
  T2_6.insert(Point(10,10));
  assert(T2_6.constrained_edges_begin() == T2_6.constrained_edges_end());

  T2_6.insert_constraint(Point(1,0.1), Point(2,0.2));
  assert(std::distance(T2_6.constrained_edges_begin(),  T2_6.constrained_edges_end()) == 1);
}
