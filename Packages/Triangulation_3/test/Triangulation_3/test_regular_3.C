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
// file          : test_regular_3.C
// revision      : 
// revision_date : 
// author(s)     : Monique Teillaud (Monique.Teillaud@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include <iostream>
#include <cassert>
#include <list>

#include <CGAL/_test_types.h>
#include <CGAL/_test_cls_regular_3.C>

bool del=true;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K> traits;

// Explicit instantiation of the whole class :
template class CGAL::Regular_triangulation_3<traits>;

int main()
{
  std::cout << " with CGAL::Regular_triangulation_euclidean_traits_3: " << std::endl;

  typedef CGAL::Regular_triangulation_3<traits>                 Cls;

  //  _test_cls_regular_3( Cls() );
  typedef traits::Bare_point Point;
  typedef traits::Weighted_point Weighted_point;

  typedef std::list<Weighted_point>                        list_point;

  // temporary version

  int n,m;
  int count = 0 ;

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
  list_point::iterator it;
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

  return 0;
}
