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


#ifndef CGAL__TEST_CLS_AFF_TRANSFORMATION_3_H
#define CGAL__TEST_CLS_AFF_TRANSFORMATION_3_H

#include <CGAL/use.h>

template <class R>
bool
_test_cls_aff_transformation_3(const R& )
{
 std::cout << "Testing class Aff_transformation_3" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;
 const bool nonexact = std::is_floating_point<FT>::value;

 typename R::Aff_transformation_3 ia;
 CGAL::Aff_transformation_3<R> a1(ia);

 RT n0 =  0;
 RT n1 =-15;
 RT n2 = 44;
 RT n3 =  4;
 RT n4 =  5;
 RT n5 = 25;
 RT n6 = -2;
 RT n7 =  8;
 RT n8 = 18;
 RT n9 =  9;
 RT n10=  3;
 RT n11=-12;
 RT n12= 20;
 RT n13=  1;
 RT n14= 35;

 CGAL::Vector_3<R> vec( n3, n8, n7,  n6 );  // (-2,-9,-4)
 CGAL::Vector_3<R> tvec;
 CGAL::Point_3<R>  pnt( n8, n1, n9, n10 );  // ( 6,-5, 3)
 CGAL::Point_3<R>  tpnt;
 CGAL::Point_3<R>  pvec = CGAL::ORIGIN + vec; CGAL_USE(pvec);
 CGAL::Vector_3<R> vpnt = pnt - CGAL::ORIGIN;

 CGAL::Point_3<R>  p1(-n3, n7, n11, n3 );   // (-1, 2,-3)
 CGAL::Point_3<R>  p2( n5, n4,-n12, n4 );   // ( 5, 1,-4)
 CGAL::Point_3<R>  p3( n1, n0, n14, n4 );   // (-3, 0, 7)
 CGAL::Point_3<R>  p4( n7, n2, n8, -n6 );   // ( 4,11, 9)
 CGAL::Weighted_point_3<R> wp4( p4, n7 );

 CGAL::Direction_3<R> d0(n13, n0, n0);
 CGAL::Direction_3<R> d1(n0, n13, n0);
 CGAL::Direction_3<R> d2(n0, n0, n13);
 CGAL::Direction_3<R> dir = (p2 - p4).direction();
 CGAL::Direction_3<R> tdir;

 CGAL::Plane_3<R>  pla(p1,p2,p3);
 CGAL::Plane_3<R>  tpla;

 CGAL::Point_3<R>   tp1;
 CGAL::Point_3<R>   tp2;
 CGAL::Point_3<R>   tp3;
 CGAL::Point_3<R>   tp4;
 CGAL::Weighted_point_3<R>   twp4;
 CGAL::Segment_3<R> seg(p1,p2);
 CGAL::Segment_3<R> tseg;
 CGAL::Ray_3<R>     ray(p3,p2);
 CGAL::Ray_3<R>     tray;
 CGAL::Line_3<R>    lin(p2,p4);
 CGAL::Line_3<R>    tlin;
 CGAL::Triangle_3<R>     tri( p2,p3,p4);
 CGAL::Triangle_3<R>     ttri;
 CGAL::Tetrahedron_3<R>  tet(p1, p2, p3, p4);
 CGAL::Tetrahedron_3<R>  ttet;
 CGAL::Sphere_3<R>       unit_sphere(CGAL::ORIGIN, 1);


 CGAL::Aff_transformation_3<R> ident( CGAL::IDENTITY );
 assert( p1 == p1.transform( ident) );
 assert( p1 == p1.transform( ident.inverse() ) );

 CGAL::Aff_transformation_3<R> gat1( n7,  n9,  n8,  n2,
                                     n5, n11, n10,  n4,
                                     n3,  n6, n12,  n2,
                                                    n3 );
 assert( p1 == (p1.transform(gat1)).transform(gat1.inverse() ) || nonexact );

 CGAL::Aff_transformation_3<R> gat2( n7,  n9,  n8,  n2,
                                     n5, n11, n10,  n4,
                                     n3,  n6, n12,  n2,
                                                    n13 );

 CGAL::Aff_transformation_3<R> gat3( n4,  n6,  n7,  n0,
                                    n12,  n8,  n8,  n0,
                                    n11,  n3,  n9,  n0,
                                                    n13 );

 CGAL::Aff_transformation_3<R> scale11( CGAL::SCALING, n2, n3 );

 CGAL::Aff_transformation_3<R> gscale(n2,  n0, n0,  n0,    // n2 = 44
                                      n0,  n2, n0,  n0,    // n3 =  4
                                      n0,  n0, n2,  n0,
                                                    n3 );


 CGAL::Aff_transformation_3<R> gtrans(n10, n0,  n0,  n8,
                                      n0, n10,  n0,  n1,
                                      n0,  n0, n10,  n9,
                                                     n10 );

 CGAL::Aff_transformation_3<R> translate( CGAL::TRANSLATION, vpnt );

 CGAL::Aff_transformation_3<R> xrefl(-n4,  n0, n0,  n0,
                                      n0,  n4, n0,  n0,
                                      n0,  n0, n4,  n0,
                                                    n4 );

 CGAL::Aff_transformation_3<R> gat4( gat3);

 CGAL::Aff_transformation_3<R> gat5( n7,  n9,  n8,
                                     n5, n11, n10,
                                     n3,  n6, n12,
                                                    n13 );

 CGAL::Aff_transformation_3<R> gat6( n4,  n6,  n7,
                                    n12,  n8,  n8,
                                    n11,  n3,  n9,
                                                    n13 );


 // even
 assert( translate.is_even() );
 assert( gtrans.is_even() );
 assert( scale11.is_even() );
 assert( gscale.is_even() );
 assert( ident.is_even() );
 assert( xrefl.is_odd() );

 // translation
 assert( translate.is_translation() );
 assert( ! scale11.is_translation() );
 assert( ! gtrans.is_translation() );

 // scaling
 assert( scale11.is_scaling() );
 assert( ! translate.is_scaling() );
 assert( ! gscale.is_scaling() );

 CGAL::Aff_transformation_3<R> a[11];

 std::cout << '.';

 a[0] = ident;
 a[1] = gat1;
 a[2] = gat2;
 a[3] = gat3;
 a[4] = scale11;
 a[5] = gscale;
 a[6] = gtrans;
 a[7] = translate;
 a[8] = xrefl;
 a[9] = gat5;
 a[10]= gat6;

 CGAL::Aff_transformation_3<R> inv;

 for (int i = 0; i< 10; i++)
 {
    // std::cout << " . " << i ;
    tp1 = p1.transform( a[i] );
    tp2 = p2.transform( a[i] );
    tp3 = p3.transform( a[i] );
    tp4 = p4.transform( a[i] );
    twp4 = wp4.transform( a[i] );
    tpla = pla.transform( a[i] );
    tseg = seg.transform( a[i] );
    tray = ray.transform( a[i] );
    tlin = lin.transform( a[i] );
    ttri = tri.transform( a[i] );
    ttet = tet.transform( a[i] );
    assert( twp4.point() == tp4 );
    assert( twp4.weight() == wp4.weight() );
    assert( tpla == CGAL::Plane_3<R>( tp1, tp2, tp3) || nonexact );
    assert( tseg == CGAL::Segment_3<R>(tp1, tp2) );
    assert( tray == CGAL::Ray_3<R>(tp3, tp2) );
    assert( tlin == CGAL::Line_3<R>(tp2, tp4) );
    assert( ttri == CGAL::Triangle_3<R>(tp2, tp3, tp4) );
    assert( ttet == CGAL::Tetrahedron_3<R>(tp1, tp2, tp3, tp4) );
    inv = a[i].inverse();
    tp4  = tp4.transform(  inv );
    twp4 = twp4.transform( inv );
    tpla = tpla.transform( inv );
    tseg = tseg.transform( inv );
    tray = tray.transform( inv );
    tlin = tlin.transform( inv );
    ttri = ttri.transform( inv );
    ttet = ttet.transform( inv );
    assert( tp4  == p4 || nonexact );
    assert( twp4 == wp4 || nonexact );
    assert( tpla == pla || nonexact );
    assert( tseg == seg || nonexact );
    assert( tray == ray || nonexact );
    assert( tlin == lin || nonexact );
    assert( ttri == tri || nonexact );
    assert( ttet == tet || nonexact );
 };

 std::cout << '.';

 // ident
 assert( vec.transform(ident) == vec );
 assert( dir.transform(ident) == dir );
 assert( pnt.transform(ident) == pnt );
 assert( pla.transform(ident) == pla || nonexact );

 // scale11 and gscale
 tpnt = pnt.transform(scale11);
 tvec = vec.transform(scale11);
 tdir = dir.transform(scale11);
 tseg = seg.transform(scale11);
 assert( tseg.squared_length() == FT(11)* FT(11)* seg.squared_length() );
 assert( FT(11)* FT(11)* FT( vec*vec ) == FT( tvec*tvec ) );
 assert( vec.direction() == tvec.direction() );
 assert( dir == tdir );
 tdir = d0.transform(scale11);
 assert( d0 == tdir);
 tpnt = pnt.transform(gscale);
 tvec = vec.transform(gscale);
 tdir = dir.transform(gscale);
 tseg = seg.transform(gscale);
 assert( tseg.squared_length() == FT(11)* FT(11)* seg.squared_length() );
 assert( FT(11)* FT(11)* FT( vec*vec ) == FT( tvec*tvec ) );
 assert( seg.transform(scale11) == seg.transform(gscale) );
 assert( vec.transform(scale11) == vec.transform(gscale) );
 assert( dir.transform(scale11) == dir.transform(gscale) );
 assert( pnt.transform(scale11) == pnt.transform(gscale) );
 assert( pla.transform(scale11) == pla.transform(gscale) || nonexact );

 // translate and gtrans
 tvec = vec.transform(translate);
 tdir = dir.transform(translate);
 tp2 = p2.transform(translate);
 tp3 = p3.transform(translate);
 assert( vec == tvec );
 assert( dir == tdir );
 assert( tp2  == p2 + vpnt );
 assert( tp3  == p3 + vpnt );
 tvec = vec.transform(gtrans);
 tdir = dir.transform(gtrans);
 tp2 = p2.transform(gtrans);
 tp3 = p3.transform(gtrans);
 assert( vec == tvec );
 assert( dir == tdir );
 assert( tp2  == p2 + vpnt );
 assert( tp3  == p3 + vpnt );
 assert( vec.transform(translate) == vec.transform(gtrans) );
 assert( dir.transform(translate) == dir.transform(gtrans) );
 assert( pnt.transform(translate) == pnt.transform(gtrans) );
 assert( pla.transform(translate) == pla.transform(gtrans) || nonexact );

 // xrefl
 tdir = d0.transform(xrefl);
 assert( tdir == -d0 );
 tdir = d1.transform(xrefl);
 assert( tdir == d1 );
 tdir = d2.transform(xrefl);
 assert( tdir == d2 );

 std::cout << '.';

 // composition
 assert( pnt.transform(xrefl).transform(xrefl) == pnt );
 assert( dir.transform(xrefl).transform(xrefl) == dir );
 assert( vec.transform(xrefl).transform(xrefl) == vec );
 assert( pla.transform(xrefl).transform(xrefl) == pla || nonexact );
 CGAL::Aff_transformation_3<R> co1 = xrefl * xrefl;
 assert( pnt.transform(xrefl).transform(xrefl) == pnt.transform(co1) );
 assert( dir.transform(xrefl).transform(xrefl) == dir.transform(co1) );
 assert( vec.transform(xrefl).transform(xrefl) == vec.transform(co1) );
 assert( pla.transform(xrefl).transform(xrefl) == pla.transform(co1) );
 co1 = gat2 * gat3;
 assert( pnt.transform(gat3).transform(gat2) == pnt.transform(co1) );
 assert( dir.transform(gat3).transform(gat2) == dir.transform(co1) );
 assert( vec.transform(gat3).transform(gat2) == vec.transform(co1) );
 assert( pla.transform(gat3).transform(gat2) == pla.transform(co1) || nonexact );
 co1 = ident * gat1;
 assert( vec.transform(gat1) == vec.transform(co1) );
 assert( dir.transform(gat1) == dir.transform(co1) );
 assert( pnt.transform(gat1) == pnt.transform(co1) );
 assert( pla.transform(gat1) == pla.transform(co1) );
 co1 = gat1 * ident;
 assert( vec.transform(gat1) == vec.transform(co1) );
 assert( dir.transform(gat1) == dir.transform(co1) );
 assert( pnt.transform(gat1) == pnt.transform(co1) );
 assert( pla.transform(gat1) == pla.transform(co1) );
 co1 = gat1 * gat1.inverse() ;
 assert( vec == vec.transform(co1) || nonexact );
 assert( dir == dir.transform(co1) || nonexact );
 assert( pnt == pnt.transform(co1) || nonexact );
 assert( pla == pla.transform(co1) || nonexact );

 assert( vec.transform( gat5 ) == vec.transform( gat2 ) );
 assert( dir.transform( gat5 ) == dir.transform( gat2 ) );

 assert( pnt.transform( gat6 ) == pnt.transform( gat3 ) );
 assert( vec.transform( gat6 ) == vec.transform( gat3 ) );
 assert( dir.transform( gat6 ) == dir.transform( gat3 ) );
 assert( pla.transform( gat6 ) == pla.transform( gat3 ) );

 // access
 FT   FTone(RT(1));
 FT   FTzero(RT(0));
 assert( gat2.cartesian(0,0) == FT(n7) / FT(n13) );
 assert( gat2.cartesian(0,1) == FT(n9) / FT(n13) );
 assert( gat2.cartesian(0,2) == FT(n8) / FT(n13) );
 assert( gat2.cartesian(0,3) == FT(n2) / FT(n13) );
 assert( gat2.cartesian(1,0) == FT(n5) / FT(n13) );
 assert( gat2.cartesian(1,1) == FT(n11) / FT(n13) );
 assert( gat2.cartesian(1,2) == FT(n10) / FT(n13) );
 assert( gat2.cartesian(1,3) == FT(n4) / FT(n13) );
 assert( gat2.cartesian(2,0) == FT(n3) / FT(n13) );
 assert( gat2.cartesian(2,1) == FT(n6) / FT(n13) );
 assert( gat2.cartesian(2,2) == FT(n12) / FT(n13) );
 assert( gat2.cartesian(2,3) == FT(n2) / FT(n13) );
 assert( gat2.cartesian(3,0) == FTzero );
 assert( gat2.cartesian(3,1) == FTzero );
 assert( gat2.cartesian(3,2) == FTzero );
 assert( gat2.cartesian(3,3) == FTone );

 assert( translate.cartesian(0,0) == FTone );
 assert( translate.cartesian(0,1) == FTzero );
 assert( translate.cartesian(0,2) == FTzero );
 assert( translate.cartesian(0,3) == FT(n8) / FT(n10) );
 assert( translate.cartesian(1,0) == FTzero );
 assert( translate.cartesian(1,1) == FTone );
 assert( translate.cartesian(1,2) == FTzero );
 assert( translate.cartesian(1,3) == FT(n1) / FT(n10) );
 assert( translate.cartesian(2,0) == FTzero );
 assert( translate.cartesian(2,1) == FTzero );
 assert( translate.cartesian(2,2) == FTone );
 assert( translate.cartesian(2,3) == FT(n9) / FT(n10) );
 assert( translate.cartesian(3,0) == FTzero );
 assert( translate.cartesian(3,1) == FTzero );
 assert( translate.cartesian(3,2) == FTzero );
 assert( translate.cartesian(3,3) == FTone );

 FT fscale = FT(n2) / FT(n3);
 assert( scale11.cartesian(0,0) == fscale );
 assert( scale11.cartesian(0,1) == FTzero );
 assert( scale11.cartesian(0,2) == FTzero );
 assert( scale11.cartesian(0,3) == FTzero );
 assert( scale11.cartesian(1,0) == FTzero );
 assert( scale11.cartesian(1,1) == fscale );
 assert( scale11.cartesian(1,2) == FTzero );
 assert( scale11.cartesian(1,3) == FTzero );
 assert( scale11.cartesian(2,0) == FTzero );
 assert( scale11.cartesian(2,1) == FTzero );
 assert( scale11.cartesian(2,2) == fscale );
 assert( scale11.cartesian(2,3) == FTzero );
 assert( scale11.cartesian(3,0) == FTzero );
 assert( scale11.cartesian(3,1) == FTzero );
 assert( scale11.cartesian(3,2) == FTzero );
 assert( scale11.cartesian(3,3) == FTone );

 assert( ident.cartesian(0,0) == FTone );
 assert( ident.cartesian(0,1) == FTzero );
 assert( ident.cartesian(0,2) == FTzero );
 assert( ident.cartesian(0,3) == FTzero );
 assert( ident.cartesian(1,0) == FTzero );
 assert( ident.cartesian(1,1) == FTone );
 assert( ident.cartesian(1,2) == FTzero );
 assert( ident.cartesian(1,3) == FTzero );
 assert( ident.cartesian(2,0) == FTzero );
 assert( ident.cartesian(2,1) == FTzero );
 assert( ident.cartesian(2,2) == FTone );
 assert( ident.cartesian(2,3) == FTzero );
 assert( ident.cartesian(3,0) == FTzero );
 assert( ident.cartesian(3,1) == FTzero );
 assert( ident.cartesian(3,2) == FTzero );
 assert( ident.cartesian(3,3) == FTone );

 // m
 assert( gat2.m(0,0) == FT(n7) / FT(n13) );
 assert( gat2.m(0,1) == FT(n9) / FT(n13) );
 assert( gat2.m(0,2) == FT(n8) / FT(n13) );
 assert( gat2.m(0,3) == FT(n2) / FT(n13) );
 assert( gat2.m(1,0) == FT(n5) / FT(n13) );
 assert( gat2.m(1,1) == FT(n11) / FT(n13) );
 assert( gat2.m(1,2) == FT(n10) / FT(n13) );
 assert( gat2.m(1,3) == FT(n4) / FT(n13) );
 assert( gat2.m(2,0) == FT(n3) / FT(n13) );
 assert( gat2.m(2,1) == FT(n6) / FT(n13) );
 assert( gat2.m(2,2) == FT(n12) / FT(n13) );
 assert( gat2.m(2,3) == FT(n2) / FT(n13) );
 assert( gat2.m(3,0) == FTzero );
 assert( gat2.m(3,1) == FTzero );
 assert( gat2.m(3,2) == FTzero );
 assert( gat2.m(3,3) == FTone );

 assert( translate.m(0,0) == FTone );
 assert( translate.m(0,1) == FTzero );
 assert( translate.m(0,2) == FTzero );
 assert( translate.m(0,3) == FT(n8) / FT(n10) );
 assert( translate.m(1,0) == FTzero );
 assert( translate.m(1,1) == FTone );
 assert( translate.m(1,2) == FTzero );
 assert( translate.m(1,3) == FT(n1) / FT(n10) );
 assert( translate.m(2,0) == FTzero );
 assert( translate.m(2,1) == FTzero );
 assert( translate.m(2,2) == FTone );
 assert( translate.m(2,3) == FT(n9) / FT(n10) );
 assert( translate.m(3,0) == FTzero );
 assert( translate.m(3,1) == FTzero );
 assert( translate.m(3,2) == FTzero );
 assert( translate.m(3,3) == FTone );

 assert( scale11.m(0,0) == fscale );
 assert( scale11.m(0,1) == FTzero );
 assert( scale11.m(0,2) == FTzero );
 assert( scale11.m(0,3) == FTzero );
 assert( scale11.m(1,0) == FTzero );
 assert( scale11.m(1,1) == fscale );
 assert( scale11.m(1,2) == FTzero );
 assert( scale11.m(1,3) == FTzero );
 assert( scale11.m(2,0) == FTzero );
 assert( scale11.m(2,1) == FTzero );
 assert( scale11.m(2,2) == fscale );
 assert( scale11.m(2,3) == FTzero );
 assert( scale11.m(3,0) == FTzero );
 assert( scale11.m(3,1) == FTzero );
 assert( scale11.m(3,2) == FTzero );
 assert( scale11.m(3,3) == FTone );

 assert( ident.m(0,0) == FTone );
 assert( ident.m(0,1) == FTzero );
 assert( ident.m(0,2) == FTzero );
 assert( ident.m(0,3) == FTzero );
 assert( ident.m(1,0) == FTzero );
 assert( ident.m(1,1) == FTone );
 assert( ident.m(1,2) == FTzero );
 assert( ident.m(1,3) == FTzero );
 assert( ident.m(2,0) == FTzero );
 assert( ident.m(2,1) == FTzero );
 assert( ident.m(2,2) == FTone );
 assert( ident.m(2,3) == FTzero );
 assert( ident.m(3,0) == FTzero );
 assert( ident.m(3,1) == FTzero );
 assert( ident.m(3,2) == FTzero );
 assert( ident.m(3,3) == FTone );

 // homogeneous
 assert( FT( gat2.hm(0,0) ) / FT( gat2.hm(3,3) ) \
         == FT(n7) / FT(n13) );
 assert( FT( gat2.hm(0,1) ) / FT( gat2.hm(3,3) ) \
         == FT(n9) / FT(n13) );
 assert( FT( gat2.hm(0,2) ) / FT( gat2.hm(3,3) ) \
         == FT(n8) / FT(n13) );
 assert( FT( gat2.hm(0,3) ) / FT( gat2.hm(3,3) ) \
         == FT(n2) / FT(n13) );
 assert( FT( gat2.hm(1,0) ) / FT( gat2.hm(3,3) ) \
         == FT(n5) / FT(n13) );
 assert( FT( gat2.hm(1,1) ) / FT( gat2.hm(3,3) ) \
         == FT(n11) / FT(n13) );
 assert( FT( gat2.hm(1,2) ) / FT( gat2.hm(3,3) ) \
         == FT(n10) / FT(n13) );
 assert( FT( gat2.hm(1,3) ) / FT( gat2.hm(3,3) ) \
         == FT(n4) / FT(n13) );
 assert( FT( gat2.hm(2,0) ) / FT( gat2.hm(3,3) ) \
         == FT(n3) / FT(n13) );
 assert( FT( gat2.hm(2,1) ) / FT( gat2.hm(3,3) ) \
         == FT(n6) / FT(n13) );
 assert( FT( gat2.hm(2,2) ) / FT( gat2.hm(3,3) ) \
         == FT(n12) / FT(n13) );
 assert( FT( gat2.hm(2,3) ) / FT( gat2.hm(3,3) ) \
         == FT(n2) / FT(n13) );
 assert( FT( gat2.hm(3,0) ) / FT( gat2.hm(3,3) ) == FTzero );
 assert( FT( gat2.hm(3,1) ) / FT( gat2.hm(3,3) ) == FTzero );
 assert( FT( gat2.hm(3,2) ) / FT( gat2.hm(3,3) ) == FTzero );
 assert( FT( gat2.hm(3,3) ) / FT( gat2.hm(3,3) ) == FTone );

 assert( FT( translate.hm(0,0) ) / FT( translate.hm(3,3) ) == FTone );
 assert( FT( translate.hm(0,1) ) / FT( translate.hm(3,3) ) == FTzero );
 assert( FT( translate.hm(0,1) ) / FT( translate.hm(3,3) ) == FTzero );
 assert( FT( translate.hm(0,3) ) / FT( translate.hm(3,3) ) \
         == FT(n8) / FT(n10) );
 assert( FT( translate.hm(1,0) ) / FT( translate.hm(3,3) ) == FTzero );
 assert( FT( translate.hm(1,1) ) / FT( translate.hm(3,3) ) == FTone );
 assert( FT( translate.hm(1,2) ) / FT( translate.hm(3,3) ) == FTzero );
 assert( FT( translate.hm(1,3) ) / FT( translate.hm(3,3) ) \
         == FT(n1) / FT(n10) );
 assert( FT( translate.hm(2,0) ) / FT( translate.hm(3,3) ) == FTzero );
 assert( FT( translate.hm(2,1) ) / FT( translate.hm(3,3) ) == FTzero );
 assert( FT( translate.hm(2,2) ) / FT( translate.hm(3,3) ) == FTone );
 assert( FT( translate.hm(2,3) ) / FT( translate.hm(3,3) ) \
         == FT(n9) / FT(n10) );
 assert( FT( translate.hm(3,0) ) / FT( translate.hm(3,3) ) == FTzero );
 assert( FT( translate.hm(3,0) ) / FT( translate.hm(3,3) ) == FTzero );
 assert( FT( translate.hm(3,1) ) / FT( translate.hm(3,3) ) == FTzero );
 assert( FT( translate.hm(3,3) ) / FT( translate.hm(3,3) ) == FTone );

 assert( FT( scale11.hm(0,0) ) / FT( scale11.hm(3,3) ) == fscale );
 assert( FT( scale11.hm(0,1) ) / FT( scale11.hm(3,3) ) == FTzero );
 assert( FT( scale11.hm(0,2) ) / FT( scale11.hm(3,3) ) == FTzero );
 assert( FT( scale11.hm(0,3) ) / FT( scale11.hm(3,3) ) == FTzero );
 assert( FT( scale11.hm(1,0) ) / FT( scale11.hm(3,3) ) == FTzero );
 assert( FT( scale11.hm(1,1) ) / FT( scale11.hm(3,3) ) == fscale );
 assert( FT( scale11.hm(1,2) ) / FT( scale11.hm(3,3) ) == FTzero );
 assert( FT( scale11.hm(1,3) ) / FT( scale11.hm(3,3) ) == FTzero );
 assert( FT( scale11.hm(2,0) ) / FT( scale11.hm(3,3) ) == FTzero );
 assert( FT( scale11.hm(2,1) ) / FT( scale11.hm(3,3) ) == FTzero );
 assert( FT( scale11.hm(2,2) ) / FT( scale11.hm(3,3) ) == fscale );
 assert( FT( scale11.hm(2,3) ) / FT( scale11.hm(3,3) ) == FTzero );
 assert( FT( scale11.hm(3,0) ) / FT( scale11.hm(3,3) ) == FTzero );
 assert( FT( scale11.hm(3,1) ) / FT( scale11.hm(3,3) ) == FTzero );
 assert( FT( scale11.hm(3,2) ) / FT( scale11.hm(3,3) ) == FTzero );
 assert( FT( scale11.hm(3,3) ) / FT( scale11.hm(3,3) ) == FTone );

 // samples
 // cartesian == m
 assert( gat1.cartesian(1,2) == gat1.m(1,2) );
 assert( scale11.cartesian(0,0) == scale11.m(0,0) );
 assert( translate.cartesian(1,3) == translate.m(1,3) );
 assert( ident.cartesian(1,2) == ident.m(1,2) );
 assert( gscale.cartesian(1,1) == gscale.m(1,1) );
 assert( gat3.cartesian(2,1) == gat3.m(2,1) );

 // m == hm/hm(3,3)
 assert( gat3.m(1,2) == FT( gat3.hm(1,2)) / FT( gat3.hm(3,3) ));
 assert( gat4.m(1,0) == FT( gat4.hm(1,0)) / FT( gat4.hm(3,3) ));
 assert( gat5.m(2,0) == FT( gat5.hm(2,0)) / FT( gat5.hm(3,3) ));
 assert( translate.m(1,3) \
         == FT( translate.hm(1,3)) / FT( translate.hm(3,3) ));
 assert( scale11.m(1,1) == FT( scale11.hm(1,1)) / FT( scale11.hm(3,3) ));

 // homogeneous == hm
 assert( gat1.homogeneous(1,2) == gat1.hm(1,2) );
 assert( gat3.homogeneous(2,3) == gat3.hm(2,3) );
 assert( gat4.homogeneous(2,3) == gat4.hm(2,3) );
 assert( scale11.homogeneous(0,0) == scale11.hm(0,0) );
 assert( translate.homogeneous(1,2) == translate.hm(1,2) );

 // orthogonal_transform
 assert(unit_sphere.orthogonal_transform(ident) == unit_sphere);
 assert(unit_sphere.orthogonal_transform(translate).center() == pnt);

 //equality
 CGAL::Aff_transformation_3<R> a2(0,1,0,1,0,1,1,0,1,0,0,1),
     a3(0,1,0,1,0,1,1,0,1,0,0,1), a4(0,0,1,1,0,0,1,1,0,0,1,1);
 assert(a2 == a3);
 assert(a3 != a4);

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_CLS_AFF_TRANSFORMATION_3_H
