// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra
 

#ifndef CGAL__TEST_FCT_DIRECTION_2_H
#define CGAL__TEST_FCT_DIRECTION_2_H

bool
ccw(int i, int j, int k)
// i between j and k
{
  return   ((j < i)&&(i < k))  // j < i < k
         ||((k < j)&&(j < i))  // k < j < i
         ||((i < k)&&(k < j)); // i < k < j
}

template <class R>
bool
_test_fct_direction_2(const R& )
{
 std::cout << "Testing functions Direction_2" ;

 typedef typename  R::RT    RT;

 RT n0 =  0;
 RT n1 =  1;
 RT n2 =  2;
 RT n3 =  3;
 RT n4 =  4;
 RT n5 =  5;
 RT n6 =  6;
 RT n8 =  8;

 CGAL::Direction_2<R> d[8];
 d[0] = CGAL::Direction_2<R>( n5, n0 );   // ( 5, 0)
 d[1] = CGAL::Direction_2<R>( n8, n3 );   // ( 8, 3)
 d[2] = CGAL::Direction_2<R>( n0, n4 );   // ( 0, 4)
 d[3] = CGAL::Direction_2<R>(-n4, n3 );   // (-4, 3)
 d[4] = CGAL::Direction_2<R>(-n2, n0 );   // (-2, 0)
 d[5] = CGAL::Direction_2<R>(-n4,-n6 );   // (-4,-6)
 d[6] = CGAL::Direction_2<R>( n0,-n1 );   // ( 0,-1)
 d[7] = CGAL::Direction_2<R>( n4,-n5 );   // ( 4,-5)

 std::cout << '.';

 assert( d[0] >= d[0] );
 assert( d[0] <= d[0] );

 std::cout << '.';

 int i;
 int j;
 int k;

 for ( i = 1; i < 8; i++ )
 {
    for ( j = 0; j+i < 8; j++)
    {
        assert( d[j] <= d[j+i] );
        assert( d[j] <  d[j+i] );
    }
    for (      ; j < 8; j++)
    {
        assert( d[j] >= d[(j+i)%8] );
        assert( d[j] >  d[(j+i)%8] );
    }
    assert( d[i] >= d[i] );
    assert( d[i] <= d[i] );
 }

 std::cout << '.';

 for (i = 0; i < 8; ++i )
     for (j = 0 ; j < 8; ++j )
         for ( k = 0; k < 8; ++k)
         {
             // std::cout << i << ' ' << j << ' ' << k << std::endl;
             if ( ccw(i,j,k) || ((j == k)&&( i != j)) )
                assert( d[i].counterclockwise_in_between(d[j],d[k]));
             else
                assert(!d[i].counterclockwise_in_between(d[j],d[k]));
             // true  if j --- i --- k along CCW rotation
             // false if i -=- j -=- k
         }

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FCT_DIRECTION_2_H
