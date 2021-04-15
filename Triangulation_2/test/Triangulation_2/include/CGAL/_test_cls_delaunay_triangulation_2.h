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
// source        : $URL$
// file          : include/CGAL/_test_cls_delaunay_triangulation_2.C
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Herve Bronnimann,Mariette Yvinec
//
// coordinator   : INRIA Sophia-Antipolis <Mariette Yvinec@sophia.inria.fr>
// ============================================================================

#include <CGAL/_test_cls_triangulation_short_2.h>

#include <cassert>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <cstdlib>

using std::rand;

template <class Del>
void
_test_delaunay_duality( const Del &T );

template <class Del>
void
_test_cls_delaunay_triangulation_2( const Del & )
{
  static_assert(std::is_nothrow_move_constructible<Del>::value,
                "move cstr is missing");
  static_assert(std::is_nothrow_move_assignable<Del>::value,
                "move assignment is missing");

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
        //          T2.insert( Point(3*m+p, m-2*p) );
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

  // test get_conflicts and star_hole
  std:: cout << "    get conflicts and star_hole" << std::endl;
  std::list<Face_handle> conflicts;
  std::list<Edge>  hole_bd;
  std::back_insert_iterator<std::list<Face_handle> > c_inserter(conflicts);
  std::back_insert_iterator<std::list<Edge> > be_inserter(hole_bd);
  std::pair<std::back_insert_iterator<std::list<Face_handle> >,
            std::back_insert_iterator<std::list<Edge> > >
    pit(c_inserter,be_inserter);
  c_inserter = T2.get_conflicts(Point(1,1,2), std::back_inserter(conflicts));
  conflicts.clear();
  pit = T2.get_conflicts_and_boundary(Point(1,1,2),
                                      std::back_inserter(conflicts),
                                      std::back_inserter(hole_bd));
  c_inserter = pit.first;
  be_inserter = pit.second;
  assert(hole_bd.size() == conflicts.size() + 2);
  conflicts.clear();
  hole_bd.clear();
  T2.get_conflicts(Point(0,1,2), std::back_inserter(conflicts));
  T2.get_boundary_of_conflicts(Point(0,1,2),
                               std::back_inserter(hole_bd));
  assert(hole_bd.size() == conflicts.size() + 2);
  conflicts.clear();
  hole_bd.clear();
  T2.get_conflicts(Point(0,0,1), std::back_inserter(conflicts));
  assert(conflicts.empty());
  conflicts.clear();
  T2.get_conflicts(Point(-1,-1,1), std::back_inserter(conflicts));
  std::size_t ns = conflicts.size();
  conflicts.clear();
  T2.find_conflicts(Point(-1,-1,1), conflicts);
  assert(conflicts.size() == ns);

  // test insertion through get_conflicts + star_hole
  conflicts.clear();
  hole_bd.clear();
  Point query(1,1,2);
  T2.get_conflicts_and_boundary(query,
                                std::back_inserter(conflicts),
                                std::back_inserter(hole_bd));

  // check the sanity of the boundary (faces are not in conflict && edges are ccw ordered)
  typename std::list<Edge>::iterator curr = hole_bd.begin(), last = --(hole_bd.end());
  Vertex_handle prev_vh = last->first->vertex(T2.ccw(last->second));
  do
  {
    assert(curr->first->vertex(T2.cw(curr->second)) == prev_vh);
    assert(!T2.test_conflict(query, curr->first));
    prev_vh = curr->first->vertex(T2.ccw(curr->second));
    ++curr;
  }
  while(curr != hole_bd.end());

  T2.star_hole (Point(1,1,2), hole_bd.begin(), hole_bd.end(),
                      conflicts.begin(), conflicts.end() );
  assert(T2.is_valid());


  // check get_conflict for a large enough point set (to use the non-recursive function)
  double x, y;
  std::vector<Point> layer_pts;

  std::ifstream in("data/layers.xy");
  assert(in);
  while(in >> x >> y)
    layer_pts.push_back(Point(x, y));

  Del T2b(layer_pts.begin(), layer_pts.end());

  conflicts.clear();
  hole_bd.clear();

  query = Point(12.25, 0.031250);
  T2b.get_conflicts_and_boundary(query,
                                 std::back_inserter(conflicts),
                                 std::back_inserter(hole_bd));

  // check the sanity of the boundary (faces are not in conflict && edges are ccw ordered)
  curr = hole_bd.begin(), last = --(hole_bd.end());
  prev_vh = last->first->vertex(T2b.ccw(last->second));
  do
  {
    assert(curr->first->vertex(T2b.cw(curr->second)) == prev_vh);
    assert(!T2b.test_conflict(query, curr->first));
    prev_vh = curr->first->vertex(T2b.ccw(curr->second));
    ++curr;
  }
  while(curr != hole_bd.end());

  /********************/
  /***** Duality ******/
  std::cout << "    duality" << std::endl;
   _test_delaunay_duality(T1);
   _test_delaunay_duality(T2);
   _test_delaunay_duality(T3);


  /**********************/
  /******* MOVE *********/
  std::cout << "    displacements" << std::endl;

  std::cout << "    degenerate cases: " << std::endl;

  Del TM_0, TM_1;
  Vertex_handle tmv1 = TM_0.insert(Point(0,0));
  Vertex_handle tmv2 = TM_0.insert(Point(1,0));

  TM_0.move_if_no_collision(tmv1, Point(2, 1));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  TM_0.move_if_no_collision(tmv1, Point(0, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

  Vertex_handle tmv3 = TM_0.insert(Point(2,1));
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv3, Point(2, 2));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 2);

  TM_0.move_if_no_collision(tmv3, Point(2, 0));
  assert(TM_0.tds().is_valid());
  assert(TM_0.is_valid());
  assert(TM_0.dimension() == 1);

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
  Vertex_handle mtmv1 = TM_0.move_if_no_collision(tmv1, Point(3, 0));
  assert(mtmv1 != tmv1);

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
    for(typename Del::Finite_vertices_iterator
         fvi = TM_1.finite_vertices_begin();
         fvi != TM_1.finite_vertices_end(); fvi++) {
      Point p = Point(rand()%30000, rand()%30000);
      vTM_1 = TM_1.move_if_no_collision(fvi, p);
      assert(TM_1.is_valid());
    }
  }

  // A simple test to see if move return the good vertex
  // when there is a collision
  Vertex_handle mvTM_1 = TM_1.move(TM_1.finite_vertices_begin(), vTM_1->point());
  assert(mvTM_1 == vTM_1);
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
