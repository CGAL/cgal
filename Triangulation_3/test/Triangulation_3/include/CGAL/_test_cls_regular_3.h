// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Monique Teillaud (Monique.Teillaud@sophia.inria.fr)

// This stuff is not used (obsolete) at the moment.
#if 0
#include <iostream>
#include <fstream>
#include <list>

template <class Triangulation>
void
_test_cls_regular_3(const Triangulation &)
{
  typedef Triangulation                      Cls;
  typedef typename Triangulation::Traits Tr;
  typedef typename Tr::Traits_base Tb;

  typedef typename Tb::Point Bare_point;
  // We assume the traits class has been tested already
  // actually, any traits is good if it has been tested
  //  typedef  typename Cls::Bare_point Point;
  typedef  typename Cls::Weighted_point Weighted_point;

  typedef std::list<Weighted_point>                        list_point;

  // temporary version

  int n,m;
  int count = 0 ;

  std::cout << " test dimension 1 " << std::endl;
  Cls T1;
  std::cout << " number of inserted points : " ;
  for ( m=0; m<5; m++) {
    if ( (m%2)== 0 ) 
      T1.insert( Weighted_point( Point( 2*m,0,0 ), 2 ) );
    else T1.insert( Weighted_point( Point( -2*m+1,0,0 ), 2 ) );
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
  for ( m=0; m<5; m++) {
    if ( (m%2)== 0 ) 
      T1.insert( Weighted_point( Point( 2*m+1,0,0 ), 5 ) );
    else T1.insert( Weighted_point( Point( -2*m+1,0,0 ), 5 ) );
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
  for ( m=0; m<10; m++) {
    if ( (m%2)== 0 ) 
      T1.insert( Weighted_point( Point( m,0,0 ), 1 ) );
    else T1.insert( Weighted_point( Point( -m,0,0 ), 1 ) );
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
  for (m=0; m<10; m++)
    for (n=0; n<10; n++) {
      T2.insert( Weighted_point( Point(m*px+n*qx, m*py+n*qy, 0), 1 ) );
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
      T2.insert( Weighted_point( Point(m*px+n*qx, m*py+n*qy, 0), -1 ) );
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
      T2.insert( Weighted_point( Point(m*px+n*qx, m*py+n*qy, 0), -2 ) );
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
      T2.insert( Weighted_point( Point(m*px+n*qx, m*py+n*qy, 0), 5 ) );
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
    for (b=0;b!=10;b++)
      for (d=0;d!=10;d++)
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
}

#endif
