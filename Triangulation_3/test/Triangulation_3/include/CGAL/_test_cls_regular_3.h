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

// This stuff is not used (obsolete) at the moment.

#include <iostream>
#include <fstream>
#include <list>
#include <type_traits>
#include <CGAL/use.h>
#include <CGAL/Testsuite/Triangulation_23/test_move_semantic.h>

template <class Triangulation>
void
_test_cls_regular_3(const Triangulation &)
{
  typedef Triangulation                       Cls;

  static_assert(std::is_nothrow_move_constructible<Cls>::value,
                "move cstr is missing");
  static_assert(std::is_nothrow_move_assignable<Cls>::value,
                "move assignment is missing");

  typedef typename Triangulation::Geom_traits Gt;
  CGAL_USE_TYPE(Gt);

  // We assume the traits class has been tested already
  // actually, any traits is good if it has been tested

  typedef typename Cls::Bare_point           Bare_point;
  typedef typename Cls::Weighted_point       Weighted_point;

  typedef std::list<Weighted_point>           list_point;

  // temporary version

  int n,m;
  int count = 0 ;

  std::cout << " test dimension 1 " << std::endl;
  Cls T1;
  std::cout << " number of inserted points : " ;
  for ( m=0; m<5; m++) {
    if ( (m%2)== 0 )
      T1.insert( Weighted_point( Bare_point( 2*m,0,0 ), 2 ) );
    else T1.insert( Weighted_point( Bare_point( -2*m+1,0,0 ), 2 ) );
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
      T1.insert( Weighted_point( Bare_point( 2*m+1,0,0 ), 5 ) );
    else T1.insert( Weighted_point( Bare_point( -2*m+1,0,0 ), 5 ) );
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
      T1.insert( Weighted_point( Bare_point( m,0,0 ), 1 ) );
    else T1.insert( Weighted_point( Bare_point( -m,0,0 ), 1 ) );
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
      T2.insert( Weighted_point( Bare_point(m*px+n*qx, m*py+n*qy, 0), 1 ) );
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
      T2.insert( Weighted_point( Bare_point(m*px+n*qx, m*py+n*qy, 0), -1 ) );
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
      T2.insert( Weighted_point( Bare_point(m*px+n*qx, m*py+n*qy, 0), -2 ) );
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
      T2.insert( Weighted_point( Bare_point(m*px+n*qx, m*py+n*qy, 0), 5 ) );
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
        lp.push_back(Weighted_point( Bare_point(a*b-d*a + (a-b)*10 +a ,
                                           a-b+d +5*b,
                                           a*a-d*d+b),
                                     a*b-a*d) );
  typename list_point::iterator it;
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

  namespace test_tr_23 = CGAL::Testsuite::Triangulation_23;
  test_tr_23::test_move_semantic(T);
}
