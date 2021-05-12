// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
//
//
// Author(s)     : Stefan Schirra


#ifndef CGAL__TEST_CLS_CIRCLE_2_H
#define CGAL__TEST_CLS_CIRCLE_2_H

#include <CGAL/Bbox_2.h>
#include <cassert>

#include <boost/type_traits/is_same.hpp>

template <class K>
void _test_construct_radical_line(const K &k) {
  typedef typename K::FT                               FT;
  typedef typename K::Point_2                          Point_2;
  typedef typename K::Line_2                           Line_2;
  typedef typename K::Circle_2                         Circle_2;
  typedef typename K::Construct_circle_2               Construct_circle_2;
  typedef typename K::Construct_radical_line_2         Construct_radical_line_2;

  Construct_circle_2 theConstruct_circle_2 = k.construct_circle_2_object();
  Construct_radical_line_2 theConstruct_radical_line_2 = k.construct_radical_line_2_object();

  std::cout << "Testing radical_line(Circle,Circle)..." << std::endl;
  Circle_2 s = theConstruct_circle_2(Point_2(0,0),1);
  for(int vx=-3;vx<4;vx++) {
    for(int vy=-3;vy<4;vy++) {
      for(int vr=1;vr<6;vr++) {
        const FT x = FT(vx);
        const FT y = FT(vy);
        const FT r = FT(vr)/FT(2);
        if(x == 0 && y == 0) continue;
        Circle_2 sl_1 = theConstruct_circle_2(Point_2(x,y),r*r);
        int d2 = (vx*vx + vy*vy);
        Line_2 p = theConstruct_radical_line_2(s, sl_1);
        Line_2 global_p = CGAL::radical_line(s, sl_1);
        assert(p == global_p);
        // No intersection
        if((d2 > (r+1)*(r+1)) || (d2 < (r-1)*(r-1))) {
          assert(!do_intersect(s, p));
        }
        // Tangent, 1 Intersection
        else if((d2 == (r+1)*(r+1)) || (d2 == (r-1)*(r-1))) {
          assert(do_intersect(s, p));
        }
        // 1 Intersection Circle
        else {
          assert(do_intersect(s, p));
        }
      }
    }
  }
}


template <class R>
bool
_test_cls_circle_2(const R& )
{
 std::cout << "Testing class Circle_2";

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Circle_2  ic;
 CGAL::Circle_2<R> c0;

 const bool nonexact = boost::is_same<FT, double>::value;

 RT n0 =  0;
 RT n1 = 16;
 RT n2 = -4;
 RT n3 =  2;
 RT n4 =  5;
 RT n5 = 10;
 RT n6 =  6;
 RT n8 =  9;
 RT n9 = 36;
 RT n10=  3;
 RT n11=-11;
 RT n12=200;
 RT n13= 13;

 CGAL::Point_2<R> p0( n1, n2, -n2);  // ( 4, -1)
 CGAL::Point_2<R> p1( n6, n8, n10);  // ( 2,  3)
 CGAL::Point_2<R> p2( n2, n0,  n2);  // ( 1,  0)
 CGAL::Point_2<R> p3( n5, n5,  n4);  // ( 2,  2)
 CGAL::Point_2<R> p4( n0, n2,  n2);  // ( 0,  1)

 CGAL::Vector_2<R> vx = p2 - CGAL::ORIGIN;
 CGAL::Vector_2<R> vy = p4 - CGAL::ORIGIN;
 CGAL::Vector_2<R> v1 = p1 - CGAL::ORIGIN;

 CGAL::Circle_2<R> c1( p0, p1, p2);
 CGAL::Circle_2<R> c2( p0, p1, p3);
 CGAL::Circle_2<R> c3( p1, p0, p2);
 const CGAL::Circle_2<R> c4( p3, FT( n9 ));      // n9 = (n6)^2
 CGAL::Vector_2<R> vx6 = vx * n6;
 CGAL::Vector_2<R> vy6 = vy * n6;
 CGAL::Circle_2<R> c5( p3 - vx6, p3 + vx6, p3 + vy6);
 CGAL::Circle_2<R> c6( c3 );
 CGAL::Circle_2<R> c7( p3, n9, CGAL::POSITIVE);
 CGAL::Circle_2<R> c8( p3, n9, CGAL::NEGATIVE);
 CGAL::Circle_2<R> cc( p3 - vx6, p3 + vx6);
 CGAL::Circle_2<R> cp( p3 - vx6, p3 + vx6, CGAL::POSITIVE);
 CGAL::Circle_2<R> cn( p3 - vx6, p3 + vx6, CGAL::NEGATIVE);
 c0 = c3;

 assert( c1 == c1 );
 assert( c1 != c2 );
 assert( c3 == c0 );
 assert( c0 == c3 );
 assert( c3 == c6 );
 assert( c7 != c8 );
 assert( c4 == c7 );
 assert( cc == cp );
 assert( cn != cp );
 assert( cc != c8 );
 assert( cc == c7 );

 assert( c5.center() == p3 );
 assert( cc.center() == p3 );
 assert( c5.squared_radius() == FT( n9 ) );
 assert( c4.squared_radius() == cc.squared_radius() );
 assert( c4 == c5 );
 assert( c4 == c7 );
 assert( c4 != c8 );
 assert( cn == cp.opposite() );
 assert( c7.opposite() == c8 );
 assert( c8.opposite() == c7 );
 assert( c1.opposite() == c3 );
 assert( c3.opposite() == c1 );
 assert( c7.orientation() == CGAL::POSITIVE );
 assert( c8.orientation() == CGAL::NEGATIVE );
 assert( c5.orientation() == CGAL::POSITIVE );
 assert( cc.orientation() == CGAL::POSITIVE );
 assert( cp.orientation() == CGAL::POSITIVE );
 assert( cn.orientation() == CGAL::NEGATIVE );

 std::cout << '.';

 assert( c4.center() == p3 );
 assert( c5.center() == p3 );
 assert( c4.squared_radius() == FT( n9 ) );
 assert( c5.squared_radius() == FT( n9 ) );
 assert( c8.squared_radius() == FT( n9 ) );

 assert( c7.bounded_side( p3 + vx*n2 ) == CGAL::ON_BOUNDED_SIDE );
 assert( c7.bounded_side( p3 + vy*n11 ) == CGAL::ON_UNBOUNDED_SIDE );
 assert( c7.bounded_side( p3 - vy6 ) == CGAL::ON_BOUNDARY );
 assert( c8.bounded_side( p3 + vx*n2 ) == CGAL::ON_BOUNDED_SIDE );
 assert( c8.bounded_side( p3 + vy*n11 ) == CGAL::ON_UNBOUNDED_SIDE );
 assert( c8.bounded_side( p3 - vy6 ) == CGAL::ON_BOUNDARY );
 assert( cc.has_on_boundary( p3 + vy6) );
 assert( cc.has_on_boundary( p3 - vx6) );

 std::cout << '.';

 CGAL::Aff_transformation_2<R>
          rotate1(CGAL::ROTATION,CGAL::Direction_2<R>(n11,n13),-n2,n12),
          rotate2(CGAL::ROTATION,CGAL::Direction_2<R>(-n8, n9),-n2,n12),
          rotate3(CGAL::ROTATION,CGAL::Direction_2<R>( n5,-n1),-n2,n12),
          rotate4(CGAL::ROTATION,CGAL::Direction_2<R>(-n5,-n11),-n2,n12);
 CGAL::Point_2<R> p6 = p2.transform( rotate1 );
 CGAL::Point_2<R> p7 = p2.transform( rotate2 );
 CGAL::Point_2<R> p8 = p2.transform( rotate3 );
 CGAL::Point_2<R> p9 = p2.transform( rotate4 );
 p6 = p6 + v1;  // v1 = ( 2, 3 )
 p7 = p7 + v1;
 p8 = p8 + v1;
 p9 = p9 + v1;
 CGAL::Circle_2<R> c10 (p6, p8, p7);
 assert( c10.orientation() == CGAL::POSITIVE );
 assert( c10.opposite().orientation() == CGAL::NEGATIVE );

 assert( c10.oriented_side(c10.center() ) == CGAL::ON_POSITIVE_SIDE );
 assert( c10.oriented_side(CGAL::ORIGIN + v1 + vx/n2 ) \
         == CGAL::ON_POSITIVE_SIDE );
 assert( c10.oriented_side(CGAL::ORIGIN + v1 + vx*n2 ) \
         == CGAL::ON_NEGATIVE_SIDE );
 assert( c10.oriented_side(p9 ) == CGAL::ON_ORIENTED_BOUNDARY || nonexact);
 assert( c10.has_on_boundary(p9) || nonexact);
 assert( c10.has_on_boundary(p4 + v1) || nonexact );
 CGAL::Point_2<R> p11( n4, n4, n3) ; // (2.5, 2.5)
 CGAL::Point_2<R> p12( n5, n5, n3) ; // ( 5 ,  5 )
 assert( c10.has_on_bounded_side( p11 ) );
 assert( ! c10.has_on_bounded_side( p12 ) );
 assert( c10.has_on_unbounded_side( p12 ) );
 assert( c10.has_on_positive_side( p11 ) );
 assert( c10.has_on_negative_side( p12 ) );
 assert( c10.opposite().has_on_negative_side( p11 ) );
 assert( c10.opposite().has_on_positive_side( p12 ) );
 assert( c10.has_on_boundary( p6 ) || nonexact);
 assert( c10.has_on_boundary( p8 ) || nonexact);

 std::cout << '.';

 CGAL::Circle_2<R> c11( p0 );
 CGAL::Circle_2<R> c12( p0, CGAL::POSITIVE );
 CGAL::Circle_2<R> c13( p0, CGAL::NEGATIVE );
 assert( c11.orientation() == CGAL::POSITIVE );
 assert( c12.orientation() == CGAL::POSITIVE );
 assert( c13.orientation() == CGAL::NEGATIVE );
 assert( c11.is_degenerate() );
 assert( c12.is_degenerate() );
 assert( c13.is_degenerate() );

 std::cout << '.';

 CGAL::Bbox_2 bb = c4.bbox();
 assert(bb.xmin() <= -4.0);
 assert(bb.xmax() >= 8.0);
 assert(bb.ymin() <= -4.0);
 assert(bb.ymax() >= 8.0);

 _test_construct_radical_line(R());

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_CIRCLE_2_H
