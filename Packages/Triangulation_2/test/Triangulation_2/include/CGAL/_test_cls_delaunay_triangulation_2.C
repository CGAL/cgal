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
// source        : $RCSfile$
// file          : include/CGAL/_test_cls_delaunay_triangulation_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann,Mariette Yvinec
//
// coordinator   : INRIA Sophia-Antipolis <Mariette Yvinec@sophia.inria.fr>
// ============================================================================

#include <iostream>
#include <iterator>
//#include <vector>

#include <CGAL/_test_cls_triangulation_short_2.C>

template <class Del>
void
_test_cls_delaunay_triangulation_2( const Del & )
{
  //typedef Del  Delaunay;
  typedef typename Del::Point                Point;
  typedef typename Del::Locate_type          Locate_type;
  typedef typename Del::Vertex_handle        Vertex_handle;
  typedef typename Del::Face_handle          Face_handle;
  typedef typename Del::Edge                 Edge;

  /***********************/
  /***** SUBCLASSES ******/
   _test_cls_triangulation_short_2( Del() );

  // Constructors
  std::cout << "    constructors(3)" << std::endl;

  // Build dummy delaunay triangulations, 1- and 2-dimensional
  Del T1;
  int m,p;
  for (m=0; m<20; m++)
      T1.insert( Point(3*m, 2*m) );
  assert( T1.is_valid() );
   
  Del T2;
  for (m=0; m<20; m++)
      for (p=0; p<20; p++)
	//	  T2.insert( Point(3*m+p, m-2*p) );
	T2.insert(Point(m,p));
  assert( T2.is_valid() );
  
  Del T3;
  // All these points are on a circle of radius 325
  Point pt[28] = {
      Point(36,323), Point(80,315), Point(91,-312),
      Point(125,300), Point(165,280), Point(195,-260), Point(204,253),
      Point(36,-323), Point(80,-315), Point(91,-312),
      Point(125,-300), Point(165,-280), Point(195,-260), Point(204,-253),
      Point(-36,-323), Point(-80,-315), Point(-91,-312),
      Point(-125,-300), Point(-165,-280), Point(-195,-260), Point(-204,-253),
      Point(-36,323), Point(-80,315), Point(-91,312),
      Point(-125,300), Point(-165,280), Point(-195,260), Point(-204,253)
  };
  for (m=0; m<28; m++) 
      T3.insert( Point(pt[m]) );
  assert( T3.is_valid() );
 

  // test nearest_vertex
  Vertex_handle vnn;
  Face_handle cible;
  int i; 
  Locate_type lt;
  cible = T2.locate(Point(0,0,1),lt,i);
  assert( lt == Del::VERTEX);
  vnn = T2.nearest_vertex(Point(1,1,3));
  assert(vnn == cible->vertex(i));
  vnn = T2.nearest_vertex(Point(0,1,3));
  assert(vnn == cible->vertex(i));
  vnn = T2.nearest_vertex(Point(0,0,1));
  assert(vnn == cible->vertex(i));
  vnn = T2.nearest_vertex(Point(-1,-1,2));
  assert(vnn == cible->vertex(i));

  // test get_conflicts
  std:: cout << "    get conflicts" << std::endl;
  std::list<Face_handle> conflicts;
  std::list<Edge>  hole_bd;
  assert(T2.get_conflicts(Point(1,1,2), std::back_inserter(conflicts)));
  conflicts.clear();	 
  assert(T2.get_conflicts_and_boundary(Point(1,1,2), 
				       std::back_inserter(conflicts),
				       std::back_inserter(hole_bd)));
  assert(hole_bd.size() == conflicts.size() + 2);
  conflicts.clear();
  hole_bd.clear();
  assert(T2.get_conflicts(Point(0,1,2), std::back_inserter(conflicts)));
  assert(T2.get_boundary_of_conflicts(Point(0,1,2), 
				      std::back_inserter(hole_bd)));
  assert(hole_bd.size() == conflicts.size() + 2);
  conflicts.clear();
  hole_bd.clear();
  assert(!T2.get_conflicts(Point(0,0,1), std::back_inserter(conflicts)));
  conflicts.clear();
  assert(T2.get_conflicts(Point(-1,-1,1), std::back_inserter(conflicts)));
  unsigned int ns = conflicts.size();
  conflicts.clear();
  assert(T2.find_conflicts(Point(-1,-1,1), conflicts));
  assert(conflicts.size() == ns);

  // test insertion through get_conflicts + star_hole
  conflicts.clear();
  hole_bd.clear();
  T2.get_conflicts_and_boundary(Point(1,1,2), 
				std::back_inserter(conflicts),
				std::back_inserter(hole_bd));
  T2.star_hole (Point(1,1,2), hole_bd.begin(), hole_bd.end(),
		      conflicts.begin(), conflicts.end() );

  

  /********************/
  /***** Duality ******/
  std::cout << "    duality" << std::endl;
   _test_delaunay_duality(T1);
   _test_delaunay_duality(T2);
   _test_delaunay_duality(T3);
}


template <class Del>
void
_test_delaunay_duality( const Del &T )
{
  typedef typename Del::Geom_traits          Gt;
  typedef typename Del::Finite_faces_iterator        Face_iterator;
  typedef typename Del::Finite_edges_iterator        Edge_iterator;
  typedef typename Del::Edge_circulator              Edge_circulator;

  // Test dual(face iterator)
  Face_iterator fit;
  for (fit = T.finite_faces_begin(); fit !=  T.finite_faces_end(); ++fit)
    {
      assert( T.side_of_oriented_circle(fit, T.dual(fit)) == 
	      CGAL::ON_POSITIVE_SIDE );
    }
  
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
