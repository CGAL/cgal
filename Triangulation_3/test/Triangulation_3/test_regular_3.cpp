// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Monique Teillaud (Monique.Teillaud@sophia.inria.fr)
//               : Mariette Yvinec (Mariette.Yvinec@sophia.inria.fr)

#include <CGAL/Regular_triangulation_3.h>
#include <iostream>
#include <cassert>
#include <list>

#include <CGAL/_test_types.h>
#include <CGAL/_test_cls_regular_3.h>
#include <CGAL/_test_cls_parallel_triangulation_3.h>

bool del=true;

typedef CGAL::Exact_predicates_inexact_constructions_kernel traits;


// Explicit instantiation of the whole class :
template class CGAL::Regular_triangulation_3<traits>;

template<typename RT>
void test_RT()
{
  typedef RT                 Cls;

  _test_cls_regular_3( Cls() );
  typedef typename RT::Bare_point                    Point;
  typedef typename RT::Weighted_point                Weighted_point;

  typedef typename Cls::Vertex_handle                Vertex_handle;
  typedef typename Cls::Cell_handle                  Cell_handle;
  typedef typename Cls::Facet                        Facet;
  typedef typename Cls::Edge                         Edge;

  typedef std::list<Weighted_point>                  list_point;
  typedef typename Cls::Finite_cells_iterator        Finite_cells_iterator;

  // temporary version

  int n, m;
  int count = 0;

  // For dimension 0, we need to check that the point of highest weight is the
  // one that finally ends up in the vertex.
  std::cout << " test dimension 0 " << std::endl;
  Cls T0;
  T0.insert(Weighted_point( Point (0,0,0), 0) );
  T0.insert(Weighted_point( Point (0,0,0), 1) );
  T0.insert(Weighted_point( Point (0,0,0), -1) );
  assert(T0.dimension() == 0);
  assert(T0.number_of_vertices() == 1);
  assert(T0.finite_vertices_begin()->point().weight() == 1);

  std::cout << " test dimension 1 " << std::endl;
  Cls T1;
  std::cout << " number of inserted points : " ;
  Weighted_point p[5];
  for ( m=0; m<5; m++) {
    if ( (m%2)== 0 )
      p[m] = Weighted_point( Point( 2*m,0,0 ), 2 );
    else
      p[m] = Weighted_point( Point( -2*m+1,0,0 ), 2 );
    T1.insert( p[m] );
    count++;
    if (count <10)
      std::cout << count << '\b' ;
    else
      if (count < 100)
        std::cout << count << '\b' << '\b' ;
      else
        std::cout << count << '\b' << '\b' << '\b' ;
    std::cout.flush();
  }
  assert( T1.is_valid() );
  std::cout << std::endl << " number of vertices : "
      << T1.number_of_vertices() << std::endl;

  std::cout << " number of inserted points : " ;
  Weighted_point q[5];
  for ( m=0; m<5; m++) {
    if ( (m%2)== 0 )
      q[m] = Weighted_point( Point( 2*m+1,0,0 ), 5 );
    else
      q[m] = Weighted_point( Point( -2*m+1,0,0 ), 5 );
    T1.insert( q[m] );
    count++;
    if (count <10)
      std::cout << count << '\b' ;
    else
      if (count < 100)
  std::cout << count << '\b' << '\b' ;
      else
  std::cout << count << '\b' << '\b' << '\b' ;
    std::cout.flush();
  }
  assert( T1.is_valid() );
  std::cout << std::endl << " number of vertices : "
      << T1.number_of_vertices() << std::endl;

  std::cout << " number of inserted points : " ;
  Weighted_point r[10];
  for ( m=0; m<10; m++) {
    if ( (m%2)== 0 )
      r[m] = Weighted_point( Point( m,0,0 ), 1 );
    else
      r[m] = Weighted_point( Point( -m,0,0 ), 1 );
    T1.insert( r[m] );
    count++;
    if (count <10)
      std::cout << count << '\b' ;
    else
      if (count < 100)
  std::cout << count << '\b' << '\b' ;
      else
  std::cout << count << '\b' << '\b' << '\b' ;
    std::cout.flush();
  }
  assert( T1.is_valid() );
  std::cout << std::endl << " number of vertices : "
      << T1.number_of_vertices() << std::endl;
  assert( T1.dimension()==1 );

  // The following is distilled from a bug report by Wulue Zhao
  // (zhao.88@osu.edu), a student of Tamal Dey.
  Point pt0(0,0,0);
  Point pt1( 1,0,0), pt2(2,0,0),  pt3(3,0,0);
  Point pt4(-1,0,0), pt5(-2,0,0), pt6(-3,0,0);

  Weighted_point wp0(pt0,10.0);
  Weighted_point wp1(pt1,0.0),  wp2(pt2,0.0),  wp3(pt3,0.0);
  Weighted_point wp4(pt4,0.0),  wp5(pt5,0.0),  wp6(pt6,0.0);

  Cls T11;

  T11.insert(wp0);
  T11.insert(wp1);
  T11.insert(wp2);
  T11.insert(wp3);
  T11.insert(wp4);
  T11.insert(wp5);
  T11.insert(wp6);

  assert(T11.is_valid());

  // And another distilled bug report from the same guy.
 {
  Point p1(-0.07, 0.04, 0.04);
  Point p2(0.09, 0.04, 0.04);
  Point p3(0.09, -0.05, 0.04);
  Point p4(0.05, -0.05, 0.04);
  Point p5(0.05, 0.0, 0.04);
  Point p6(-0.07, 0.0, 0.04);
  Point p7(-0.07, 0.04, -0.04);
  Point p8(0.09, 0.04, -0.04);
  Point p9(0.09, -0.05, -0.04);
  Point p10(0.05, -0.05, -0.04);
  Point p11(0.05, 0.0, -0.04);
  Point p12(-0.07, 0.0, -0.04);

  Weighted_point wp1(p1,0);
  Weighted_point wp2(p2,0);
  Weighted_point wp3(p3,0);
  Weighted_point wp4(p4,0);
  Weighted_point wp5(p5,0);
  Weighted_point wp6(p6,0);
  Weighted_point wp7(p7,0);
  Weighted_point wp8(p8,0);
  Weighted_point wp9(p9,0);
  Weighted_point wp10(p10,0);
  Weighted_point wp11(p11,0);
  Weighted_point wp12(p12,0);
  Weighted_point wp13(p3,0.3); // wp13 has the same coordinates with wp3

  Cls T111;

  T111.insert(wp1);
  T111.insert(wp2);
  T111.insert(wp3);
  T111.insert(wp13); // it doesnot work inserting wp13 here
  T111.insert(wp4);
  T111.insert(wp5);
  T111.insert(wp6);
  T111.insert(wp7);
  T111.insert(wp8);
  T111.insert(wp9);
  T111.insert(wp10);
  T111.insert(wp11);
  T111.insert(wp12);

  assert(T111.is_valid());
 }

  std::cout << " test dimension 2 " << std::endl;
  std::cout << " number of inserted points : " ;
  Cls T2;

  count = 0 ;
  int px=1, py=1;
  int qx=-1, qy=2;
  Weighted_point s[400];
  for (m=0; m<10; m++)
    for (n=0; n<10; n++) {
      s[m+20*n] = Weighted_point( Point(m*px+n*qx, m*py+n*qy, 0), 1 );
      T2.insert( s[m+20*n] );
      count++;
      if (count <10)
  std::cout << count << '\b' ;
      else
  if (count < 100)
    std::cout << count << '\b' << '\b' ;
  else
    std::cout << count << '\b' << '\b' << '\b' ;
      std::cout.flush();
    }
  for (m=10; m<20; m++)
    for (n=0; n<10; n++) {
      s[m+20*n] = Weighted_point( Point(m*px+n*qx, m*py+n*qy, 0), -1 );
      T2.insert( s[m+20*n] );
      count++;
      if (count <10)
  std::cout << count << '\b' ;
      else
  if (count < 100)
    std::cout << count << '\b' << '\b' ;
  else
    std::cout << count << '\b' << '\b' << '\b' ;
      std::cout.flush();
    }
  for (m=0; m<10; m++)
    for (n=10; n<20; n++) {
      s[m+20*n] = Weighted_point( Point(m*px+n*qx, m*py+n*qy, 0), -2 );
      T2.insert( s[m+20*n] );
      count++;
      if (count <10)
  std::cout << count << '\b' ;
      else
  if (count < 100)
    std::cout << count << '\b' << '\b' ;
  else
    std::cout << count << '\b' << '\b' << '\b' ;
      std::cout.flush();
    }
  for (m=10; m<20; m++)
    for (n=10; n<20; n++) {
      s[m+20*n] = Weighted_point( Point(m*px+n*qx, m*py+n*qy, 0), 5 );
      T2.insert( s[m+20*n] );
      count++;
      if (count <10)
  std::cout << count << '\b' ;
      else
  if (count < 100)
    std::cout << count << '\b' << '\b' ;
  else
    std::cout << count << '\b' << '\b' << '\b' ;
      std::cout.flush();
    }

  std::cout << std::endl << " number of vertices : "
      << T2.number_of_vertices() << std::endl;
  assert( T2.dimension()==2 );
  assert( T2.is_valid() );

 // dimension 3
  std::cout << " test dimension 3" << std::endl;
  Cls T;

  list_point lp;
  int a, b, d;
  for (a=0;a!=10;a++)
    //    for (b=0;b!=10;b++)
    for (b=0;b!=5;b++)
      //      for (d=0;d!=10;d++)
      for (d=0;d!=5;d++)
  lp.push_back(Weighted_point( Point(a*b-d*a + (a-b)*10 +a ,
             a-b+d +5*b,
             a*a-d*d+b),
             a*b-a*d) );
  typename list_point::iterator it;
  count = 0 ;
  std::cout << " number of inserted points : " ;
  for (it=lp.begin(); it!=lp.end(); ++it){
    count++;
    T.insert(*it);
    if (count <10)
      std::cout << count << '\b' ;
    else
      if (count < 100)
        std::cout << count << '\b' << '\b' ;
      else
        if (count < 1000)
          std::cout << count << '\b' << '\b' << '\b' ;
        else
    std::cout << count << std::endl;
    std::cout.flush();
  }
  std::cout << std::endl;

  std::cout << " number of vertices : "
      << T.number_of_vertices() << std::endl;
  assert(T.is_valid());
  assert(T.dimension()==3);

  T.clear();
  std::cout << " test iterator range insert" << std::endl;
  T.insert (lp.begin(), lp.end());

  std::cout << " number of vertices : "
      << T.number_of_vertices() << std::endl;
  assert(T.is_valid());
  assert(T.dimension()==3);


    //test nearest_power_vertex
  std::cout << " test nearest_power_vertex " << std::endl;
  Point pp1(0.0, 0.0, 0.0);
  Point pp2(1.0, 0.0, 0.0);
  Point pp3(0.0, 1.0, 0.0);
  Point pp4(0.0, 0.0, 1.0);
  Point pp5(1.0, 1.0, 0.0);
  Point pp6(0.0, 1.0, 1.0);
  Point pp7(1.0, 0.0, 1.0);
  Point pp8(1.0, 1.0, 1.0);
  Point pp9(0.5, 0.5, 0.5);

  Weighted_point wpp1(pp1, 1.0);
  Weighted_point wpp2(pp2, 2.0);
  Weighted_point wpp3(pp3, 1.0);
  Weighted_point wpp4(pp4, 4.0);
  Weighted_point wpp5(pp5, 1.0);
  Weighted_point wpp6(pp6, 1.0);
  Weighted_point wpp7(pp7, 1.0);
  Weighted_point wpp8(pp8, 8.0);
  Weighted_point wpp9(pp9, -8.0);

  Cls T3;

  T3.insert(wpp1);
  Vertex_handle v2 = T3.insert(wpp2);
  assert( T3.nearest_power_vertex(Point(0.5,0.5,0.5)) == v2);

  T3.insert(wpp3);
  Vertex_handle v4 = T3.insert(wpp4);
  assert( T3.nearest_power_vertex(Point(0.5,0.5,0.5)) == v4);

  T3.insert(wpp5);
  T3.insert(wpp6);
  T3.insert(wpp7);

  Vertex_handle v8 = T3.insert(wpp8);
  Point query(0.5,0.5,0.5);
  assert(T3.nearest_power_vertex(query) == v8);
  assert(T3.nearest_power_vertex_in_cell(query ,v8->cell()) == v8);

  Vertex_handle v9 = T3.insert(wpp9);
  assert(v9 == Vertex_handle()); // hidden point

  // test dual
  std::cout << " test dual member functions" << std::endl;
  Finite_cells_iterator fcit = T3.finite_cells_begin();
  for( ; fcit != T3.finite_cells_end(); ++fcit) {
    Point cc = T3.dual(fcit);
    Vertex_handle ncc = T3.nearest_power_vertex(cc);
    assert(fcit->has_vertex(ncc));
  }

  // test Gabriel
  std::cout << " test is_Gabriel " << std::endl;
  Point q0(0.,0.,0.);
  Point q1(2.,0.,0.);
  Point q2(0.,2.,0.);
  Point q3(0.,0.,2.);

  Weighted_point wq0(q0,0.);
  Weighted_point wq1(q1,0.);
  Weighted_point wq2(q2,0.);
  Weighted_point wq3(q3,0.);
  Weighted_point wq01(q0,2.);

  Cls T4;
  Vertex_handle v0 = T4.insert(wq0);
  Vertex_handle v1 = T4.insert(wq1);
  v2 = T4.insert(wq2);
  Vertex_handle v3 = T4.insert(wq3);
  Cell_handle c;
  int i,j,k,l;
  assert(T4.is_facet(v0,v1,v2,c,j,k,l));
  i = 6 - (j+k+l);
  Facet f = std::make_pair(c,i);
  assert(T4.is_Gabriel(c,i));
  assert(T4.is_Gabriel(f));
  assert(T4.is_facet(v1,v2,v3,c,j,k,l));
  i = 6 - (j+k+l);
  assert(!T4.is_Gabriel(c,i));
  assert(T4.is_edge(v0,v1,c,i,j));
  assert(T4.is_Gabriel(c,i,j));
  Edge e = make_triple(c,i,j);
  assert(T4.is_Gabriel(e));
  assert(T4.is_edge(v2,v3,c,i,j));
  assert(T4.is_Gabriel(c,i,j));

  Vertex_handle v01 = T4.insert(wq01);
  (void) v01; // kill warning
  assert(T4.is_edge(v2,v3,c,i,j));
  assert(!T4.is_Gabriel(c,i,j));

  Weighted_point wwq0(q0,0.);
  Weighted_point wwq1(q1,0.);
  Weighted_point wwq2(q2,0.);
  Weighted_point wwq3(q3,5.);
  Cls T5;
  v0 = T5.insert(wwq0);
  v1 = T5.insert(wwq1);
  v2 = T5.insert(wwq2);
  v3 = T5.insert(wwq3);
  assert(T5.nearest_power_vertex(v3->point().point()) == v3);
  assert(T5.nearest_power_vertex(v0->point().point()) == v3);
  assert(T5.is_Gabriel(v3));
  assert(!T5.is_Gabriel(v0));
}

int main()
{
  test_RT<CGAL::Regular_triangulation_3<traits> >();

#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Spatial_lock_grid_3<
    CGAL::Tag_priority_blocking>                      Lock_ds;
  typedef CGAL::Triangulation_data_structure_3<
    CGAL::Regular_triangulation_vertex_base_3<traits>,
    CGAL::Regular_triangulation_cell_base_3<traits>,
    CGAL::Parallel_tag >                                    Tds_parallel;
  typedef CGAL::Regular_triangulation_3<
    traits, Tds_parallel, Lock_ds>                    RT_parallel;

  // The following test won't do things in parallel since it doesn't provide a lock data structure
  test_RT<RT_parallel>();

  // This test performs parallel operations
  _test_cls_parallel_triangulation_3( RT_parallel() );
#endif

  std::cout << " quit " << std::endl;
  return 0;
}

