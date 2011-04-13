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

#include <CGAL/basic.h>
#include <iostream>
#include <cassert>
#include <list>
#include <vector>
#include <CGAL/triple.h>

#include <CGAL/_test_types.h>

#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include <CGAL/Regular_triangulation_3.h>

#include <CGAL/_test_cls_regular_3.C>

bool del=true;

int main()
{
  std::cout << " with CGAL::Regular_triangulation_euclidean_traits_3: " << std::endl;

  typedef CGAL::Regular_triangulation_euclidean_traits_3<Test_rep_cartesian> traits;

  typedef CGAL::Regular_triangulation_3<traits>                 Cls;

  //  _test_cls_regular_3( Cls() );
  typedef traits::Bare_point Point;
  typedef traits::Weighted_point Weighted_point;

  typedef std::list<Weighted_point>                        list_point;

  // temporary version

  int n,m;
  int count = 0 ;

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
  for (it=lp.begin(); it!=lp.end();it++){
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

  std::cout << " number of vertices : " 
	    << T.number_of_vertices() << std::endl;
  assert(T.is_valid());
  assert(T.dimension()==3);

  return 0;
}
