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
 

#ifndef CGAL__TEST_CLS_AFF_TRANSFORMATION_2_H
#define CGAL__TEST_CLS_AFF_TRANSFORMATION_2_H

template <class R>
bool
_test_cls_aff_transformation_2(const R& )
{
 std::cout << "Testing class Aff_transformation_2" ;

 typedef typename  R::RT    RT;
 typedef typename  R::FT    FT;

 typename R::Aff_transformation_2 ia;
 CGAL::Aff_transformation_2<R> a1(ia);

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
 RT n100 = 100;

 CGAL::Vector_2<R> vec( n3, n8,  n6 );      // (-2,-9)
 CGAL::Vector_2<R> tvec;
 CGAL::Point_2<R>  pnt( n8, n1, n10 );      // ( 6,-5)
 CGAL::Point_2<R>  tpnt;
 CGAL::Point_2<R>  pvec = CGAL::ORIGIN + vec;
 CGAL::Vector_2<R> vpnt = pnt - CGAL::ORIGIN;

 CGAL::Point_2<R>  p1(-n3, n7, n3 );        // (-1, 2)
 CGAL::Point_2<R>  p2( n5, n4, n4 );        // ( 5, 1)
 CGAL::Point_2<R>  p3( n1, n0, n4 );        // (-3, 0)
 CGAL::Point_2<R>  p4( n7, n2,-n6 );        // ( 4,11)

 CGAL::Direction_2<R> d0(n13, n0);
 CGAL::Direction_2<R> d1(n0, n13);
 CGAL::Direction_2<R> dir = (p2 - p4).direction();
 CGAL::Direction_2<R> tdir;

 CGAL::Point_2<R>   tp1;
 CGAL::Point_2<R>   tp2;
 CGAL::Point_2<R>   tp3;
 CGAL::Point_2<R>   tp4;
 CGAL::Segment_2<R> seg(p1,p2);
 CGAL::Segment_2<R> tseg;
 CGAL::Ray_2<R>     ray(p3,p2);
 CGAL::Ray_2<R>     tray;
 CGAL::Line_2<R>    lin(p2,p4);
 CGAL::Line_2<R>    tlin;
 CGAL::Triangle_2<R>     tri( p2,p3,p4);
 CGAL::Triangle_2<R>     ttri;

 CGAL::Circle_2<R>         circ(p2, p3, p4);
 CGAL::Circle_2<R>         tcirc;
 CGAL::Iso_rectangle_2<R>  isor(p3, p4);
 CGAL::Iso_rectangle_2<R>  tisor;

 CGAL::Aff_transformation_2<R> ident( CGAL::IDENTITY );

 CGAL::Aff_transformation_2<R> gat1( n7,  n9,  n2,
                                     n5, n11,  n4,
                                               n3 );

 CGAL::Aff_transformation_2<R> gat2( n7,  n9,  n2,
                                     n5, n11,  n4,
                                               n13 );

 CGAL::Aff_transformation_2<R> gat3( n4,  n6,  n0,
                                    n12,  n8,  n0,
                                               n13 );

 CGAL::Aff_transformation_2<R> scale11( CGAL::SCALING, n2, n3 );

 CGAL::Aff_transformation_2<R> gscale(n2,  n0, n0,
                                      n0,  n2, n0,
                                               n3 );

 CGAL::Aff_transformation_2<R> gtrans(n10, n0, n8,
                                      n0, n10, n1,
                                               n10 );

 CGAL::Aff_transformation_2<R> translate( CGAL::TRANSLATION, vpnt );

 CGAL::Aff_transformation_2<R> xrefl(-n4,  n0, n0,
                                      n0,  n4, n0,
                                               n4 );

 CGAL::Aff_transformation_2<R> gat4( gat3);

 CGAL::Aff_transformation_2<R> gat5( n7,  n9,
                                     n5, n11,
                                               n13 );

 CGAL::Aff_transformation_2<R> gat6( n4,  n6,
                                    n12,  n8,
                                               n13 );

 CGAL::Aff_transformation_2<R> rot90( CGAL::ROTATION, d1, n13, n100 );

 CGAL::Aff_transformation_2<R> rot2( CGAL::ROTATION, dir, n13, n100 );

 CGAL::Aff_transformation_2<R> rot3( CGAL::ROTATION, RT(3),RT(4),RT(5));



 CGAL::Aff_transformation_2<R> a[14];


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
 a[11]= rot90;
 a[12]= rot2;
 a[13]= rot3;

 CGAL::Aff_transformation_2<R> inv;

 for (int i = 0; i< 14; i++)
 {
    tp1 = p1.transform( a[i] );
    tp2 = p2.transform( a[i] );
    tp3 = p3.transform( a[i] );
    tp4 = p4.transform( a[i] );
    tseg = seg.transform( a[i] );
    tray = ray.transform( a[i] );
    tlin = lin.transform( a[i] );
    ttri = tri.transform( a[i] );
    tisor= isor.transform( a[i]);
    assert( tseg == CGAL::Segment_2<R>(tp1, tp2) );
    assert( tray == CGAL::Ray_2<R>(tp3, tp2) );
    assert( tlin == CGAL::Line_2<R>(tp2, tp4) );
    assert( ttri == CGAL::Triangle_2<R>(tp2, tp3, tp4) );
    assert( tisor== CGAL::Iso_rectangle_2<R>( tp3, tp4 ) );

    inv = a[i].inverse();
    tp4  = tp4.transform(  inv );
    tseg = tseg.transform( inv );
    tray = tray.transform( inv );
    tlin = tlin.transform( inv );
    ttri = ttri.transform( inv );
    assert( tp4  == p4 );
    assert( tseg == seg );
    assert( tray == ray );
    assert( tlin == lin );
    assert( ttri == tri );
 };

 std::cout << '.';

 // ident
 assert( vec.transform(ident) == vec );
 assert( dir.transform(ident) == dir );
 assert( pnt.transform(ident) == pnt );
 assert( lin.transform(ident) == lin );

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
 assert( lin.transform(scale11) == lin.transform(gscale) );

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
 assert( lin.transform(translate) == lin.transform(gtrans) );

 // xrefl
 tdir = d0.transform(xrefl);
 assert( tdir == -d0 );
 tdir = d1.transform(xrefl);
 assert( tdir == d1 );

 std::cout << '.';

 // composition
 assert( pnt.transform(xrefl).transform(xrefl) == pnt );
 assert( dir.transform(xrefl).transform(xrefl) == dir );
 assert( vec.transform(xrefl).transform(xrefl) == vec );
 assert( lin.transform(xrefl).transform(xrefl) == lin );
 CGAL::Aff_transformation_2<R> co1 = xrefl * xrefl;
 assert( pnt.transform(xrefl).transform(xrefl) == pnt.transform(co1) );
 assert( dir.transform(xrefl).transform(xrefl) == dir.transform(co1) );
 assert( vec.transform(xrefl).transform(xrefl) == vec.transform(co1) );
 assert( lin.transform(xrefl).transform(xrefl) == lin.transform(co1) );
 co1 = gat2 * gat3;
 assert( pnt.transform(gat3).transform(gat2) == pnt.transform(co1) );
 assert( dir.transform(gat3).transform(gat2) == dir.transform(co1) );
 assert( vec.transform(gat3).transform(gat2) == vec.transform(co1) );
 assert( lin.transform(gat3).transform(gat2) == lin.transform(co1) );
 co1 = ident * gat1;
 assert( vec.transform(gat1) == vec.transform(co1) );
 assert( dir.transform(gat1) == dir.transform(co1) );
 assert( pnt.transform(gat1) == pnt.transform(co1) );
 assert( lin.transform(gat1) == lin.transform(co1) );
 co1 = gat1 * ident;
 assert( vec.transform(gat1) == vec.transform(co1) );
 assert( dir.transform(gat1) == dir.transform(co1) );
 assert( pnt.transform(gat1) == pnt.transform(co1) );
 assert( lin.transform(gat1) == lin.transform(co1) );
 co1 = gat1 * gat1.inverse() ;
 assert( vec == vec.transform(co1) );
 assert( pnt == pnt.transform(co1) );
 assert( dir == dir.transform(co1) );
 assert( lin == lin.transform(co1) );

 // even
 assert( translate.is_even() );
 assert( gtrans.is_even() );
 assert( gscale.is_even() );
 assert( scale11.is_even() );
 assert( ident.is_even() );
 assert( rot90.is_even() );
 assert( rot2.is_even() );
 assert( rot3.is_even() );
 assert( xrefl.is_odd() );

 // rotation
 assert( d0.transform( rot90 ) == d1 );
 assert( d1.transform( rot90.inverse() ) == d0 );
 assert( d0.transform( rot3 ) == CGAL::Direction_2<R>( RT(4), RT(3)) );
 co1 = rot3 * rot90;
 assert( d1.transform( rot3) == d0.transform( co1 ) );
 co1 = rot2 * rot90;
 assert( d1.transform( rot2) == d0.transform( co1 ) );
 co1 = rot90 * rot2;
 assert( d1.transform( rot2) == d0.transform( co1 ) );
 co1 = rot90 * rot90 * rot90 * rot90;
 assert( vec == vec.transform(co1) );
 assert( dir == dir.transform(co1) );
 assert( pnt == pnt.transform(co1) );
 assert( lin == lin.transform(co1) );
 co1 = rot3 * rot3 * rot3.inverse();
 assert( vec.transform(rot3) == vec.transform(co1) );
 assert( dir.transform(rot3) == dir.transform(co1) );
 assert( pnt.transform(rot3) == pnt.transform(co1) );
 assert( lin.transform(rot3) == lin.transform(co1) );

 //circle
 tp2 = p2.transform( translate );
 tp3 = p3.transform( translate );
 tp4 = p4.transform( translate );
 tcirc = circ.orthogonal_transform( translate );
 assert( tcirc == CGAL::Circle_2<R>( tp2, tp3, tp4 ) );
 tp2 = p2.transform( xrefl );
 tp3 = p3.transform( xrefl );
 tp4 = p4.transform( xrefl );
 tcirc = circ.orthogonal_transform( xrefl );
 assert( tcirc == CGAL::Circle_2<R>( tp2, tp3, tp4 ) );
 tp2 = p2.transform( rot3 );
 tp3 = p3.transform( rot3 );
 tp4 = p4.transform( rot3 );
 tcirc = circ.orthogonal_transform( rot3 );
 assert( tcirc == CGAL::Circle_2<R>( tp2, tp3, tp4 ) );


 // copy
 assert( vec.transform( gat5 ) == vec.transform( gat2 ) );
 assert( dir.transform( gat5 ) == dir.transform( gat2 ) );

 assert( pnt.transform( gat6 ) == pnt.transform( gat3 ) );
 assert( vec.transform( gat6 ) == vec.transform( gat3 ) );
 assert( dir.transform( gat6 ) == dir.transform( gat3 ) );
 assert( lin.transform( gat6 ) == lin.transform( gat3 ) );

 // access
 // general form
 FT   FTone(RT(1));
 FT   FTzero(RT(0));
 assert( gat2.cartesian(0,0) == FT(n7) / FT(n13) );
 assert( gat2.cartesian(0,1) == FT(n9) / FT(n13) );
 assert( gat2.cartesian(0,2) == FT(n2) / FT(n13) );
 assert( gat2.cartesian(1,0) == FT(n5) / FT(n13) );
 assert( gat2.cartesian(1,1) == FT(n11) / FT(n13) );
 assert( gat2.cartesian(1,2) == FT(n4) / FT(n13) );
 assert( gat2.cartesian(2,0) == FTzero );
 assert( gat2.cartesian(2,1) == FTzero );
 assert( gat2.cartesian(2,2) == FTone );

 assert( gat6.cartesian(0,0) == FT(n4) / FT(n13) );
 assert( gat6.cartesian(0,1) == FT(n6) / FT(n13) );
 assert( gat6.cartesian(0,2) == FTzero );
 assert( gat6.cartesian(1,0) == FT(n12) / FT(n13) );
 assert( gat6.cartesian(1,1) == FT(n8) / FT(n13) );
 assert( gat6.cartesian(1,2) == FTzero );
 assert( gat6.cartesian(2,0) == FTzero );
 assert( gat6.cartesian(2,1) == FTzero );
 assert( gat6.cartesian(2,2) == FTone );

 // translation
 assert( translate.cartesian(0,0) == FTone );
 assert( translate.cartesian(0,1) == FTzero );
 assert( translate.cartesian(0,2) == FT(n8) / FT(n10) );
 assert( translate.cartesian(1,0) == FTzero );
 assert( translate.cartesian(1,1) == FTone );
 assert( translate.cartesian(1,2) == FT(n1) / FT(n10) );
 assert( translate.cartesian(2,0) == FTzero );
 assert( translate.cartesian(2,1) == FTzero );
 assert( translate.cartesian(2,2) == FTone );

 // rotation
 FT f3o5 = FT(RT(3)) / FT(RT(5));
 FT f4o5 = FT(RT(4)) / FT(RT(5));
 assert( rot3.cartesian(0,0) == f4o5 );
 assert( rot3.cartesian(0,1) == -f3o5 );
 assert( rot3.cartesian(0,2) == FTzero );
 assert( rot3.cartesian(1,0) == f3o5 );
 assert( rot3.cartesian(1,1) == f4o5 );
 assert( rot3.cartesian(1,2) == FTzero );
 assert( rot3.cartesian(2,0) == FTzero );
 assert( rot3.cartesian(2,1) == FTzero );
 assert( rot3.cartesian(2,2) == FTone );

 // scaling
 FT fscale = FT(n2) / FT(n3);
 assert( scale11.cartesian(0,0) == fscale );
 assert( scale11.cartesian(0,1) == FTzero );
 assert( scale11.cartesian(0,2) == FTzero );
 assert( scale11.cartesian(1,0) == FTzero );
 assert( scale11.cartesian(1,1) == fscale );
 assert( scale11.cartesian(1,2) == FTzero );
 assert( scale11.cartesian(2,0) == FTzero );
 assert( scale11.cartesian(2,1) == FTzero );
 assert( scale11.cartesian(2,2) == FTone );

 // ident
 assert( ident.cartesian(0,0) == FTone );
 assert( ident.cartesian(0,1) == FTzero );
 assert( ident.cartesian(0,2) == FTzero );
 assert( ident.cartesian(1,0) == FTzero );
 assert( ident.cartesian(1,1) == FTone );
 assert( ident.cartesian(1,2) == FTzero );
 assert( ident.cartesian(2,0) == FTzero );
 assert( ident.cartesian(2,1) == FTzero );
 assert( ident.cartesian(2,2) == FTone );

 // same with m
 // general form
 assert( gat2.m(0,0) == FT(n7) / FT(n13) );
 assert( gat2.m(0,1) == FT(n9) / FT(n13) );
 assert( gat2.m(0,2) == FT(n2) / FT(n13) );
 assert( gat2.m(1,0) == FT(n5) / FT(n13) );
 assert( gat2.m(1,1) == FT(n11) / FT(n13) );
 assert( gat2.m(1,2) == FT(n4) / FT(n13) );
 assert( gat2.m(2,0) == FTzero );
 assert( gat2.m(2,1) == FTzero );
 assert( gat2.m(2,2) == FTone );

 assert( gat6.m(0,0) == FT(n4) / FT(n13) );
 assert( gat6.m(0,1) == FT(n6) / FT(n13) );
 assert( gat6.m(0,2) == FTzero );
 assert( gat6.m(1,0) == FT(n12) / FT(n13) );
 assert( gat6.m(1,1) == FT(n8) / FT(n13) );
 assert( gat6.m(1,2) == FTzero );
 assert( gat6.m(2,0) == FTzero );
 assert( gat6.m(2,1) == FTzero );
 assert( gat6.m(2,2) == FTone );

 // translation
 assert( translate.m(0,0) == FTone );
 assert( translate.m(0,1) == FTzero );
 assert( translate.m(0,2) == FT(n8) / FT(n10) );
 assert( translate.m(1,0) == FTzero );
 assert( translate.m(1,1) == FTone );
 assert( translate.m(1,2) == FT(n1) / FT(n10) );
 assert( translate.m(2,0) == FTzero );
 assert( translate.m(2,1) == FTzero );
 assert( translate.m(2,2) == FTone );

 // rotation
 assert( rot3.m(0,0) == f4o5 );
 assert( rot3.m(0,1) == -f3o5 );
 assert( rot3.m(0,2) == FTzero );
 assert( rot3.m(1,0) == f3o5 );
 assert( rot3.m(1,1) == f4o5 );
 assert( rot3.m(1,2) == FTzero );
 assert( rot3.m(2,0) == FTzero );
 assert( rot3.m(2,1) == FTzero );
 assert( rot3.m(2,2) == FTone );

 // scaling
 assert( scale11.m(0,0) == fscale );
 assert( scale11.m(0,1) == FTzero );
 assert( scale11.m(0,2) == FTzero );
 assert( scale11.m(1,0) == FTzero );
 assert( scale11.m(1,1) == fscale );
 assert( scale11.m(1,2) == FTzero );
 assert( scale11.m(2,0) == FTzero );
 assert( scale11.m(2,1) == FTzero );
 assert( scale11.m(2,2) == FTone );

 // ident
 assert( ident.m(0,0) == FTone );
 assert( ident.m(0,1) == FTzero );
 assert( ident.m(0,2) == FTzero );
 assert( ident.m(1,0) == FTzero );
 assert( ident.m(1,1) == FTone );
 assert( ident.m(1,2) == FTzero );
 assert( ident.m(2,0) == FTzero );
 assert( ident.m(2,1) == FTzero );
 assert( ident.m(2,2) == FTone );

 // homogeneous and hm are not like Quebec (== unique)
 // general form
 assert( FT( gat2.hm(0,0) ) / FT( gat2.hm(2,2) ) \
         == FT(n7) / FT(n13) );
 assert( FT( gat2.hm(0,1) ) / FT( gat2.hm(2,2) ) \
         == FT(n9) / FT(n13) );
 assert( FT( gat2.hm(0,2) ) / FT( gat2.hm(2,2) ) \
         == FT(n2) / FT(n13) );
 assert( FT( gat2.hm(1,0) ) / FT( gat2.hm(2,2) ) \
         == FT(n5) / FT(n13) );
 assert( FT( gat2.hm(1,1) ) / FT( gat2.hm(2,2) ) \
         == FT(n11) / FT(n13) );
 assert( FT( gat2.hm(1,2) ) / FT( gat2.hm(2,2) ) \
         == FT(n4) / FT(n13) );
 assert( FT( gat2.hm(2,0) ) / FT( gat2.hm(2,2) ) == FTzero );
 assert( FT( gat2.hm(2,1) ) / FT( gat2.hm(2,2) ) == FTzero );
 assert( FT( gat2.hm(2,2) ) / FT( gat2.hm(2,2) ) == FTone );

 assert( FT( gat6.hm(0,0) ) / FT( gat6.hm(2,2) ) \
         == FT(n4) / FT(n13) );
 assert( FT( gat6.hm(0,1) ) / FT( gat6.hm(2,2) ) \
         == FT(n6) / FT(n13) );
 assert( FT( gat6.hm(0,2) ) / FT( gat6.hm(2,2) ) == FTzero );
 assert( FT( gat6.hm(1,0) ) / FT( gat6.hm(2,2) ) \
         == FT(n12) / FT(n13) );
 assert( FT( gat6.hm(1,1) ) / FT( gat6.hm(2,2) ) \
         == FT(n8) / FT(n13) );
 assert( FT( gat6.hm(1,2) ) / FT( gat6.hm(2,2) ) == FTzero );
 assert( FT( gat6.hm(2,0) ) / FT( gat6.hm(2,2) ) == FTzero );
 assert( FT( gat6.hm(2,1) ) / FT( gat6.hm(2,2) ) == FTzero );
 assert( FT( gat6.hm(2,2) ) / FT( gat6.hm(2,2) ) == FTone );

 // translation
 assert( FT( translate.hm(0,0) ) / FT( translate.hm(2,2) ) == FTone );
 assert( FT( translate.hm(0,1) ) / FT( translate.hm(2,2) ) == FTzero );
 assert( FT( translate.hm(0,2) ) / FT( translate.hm(2,2) ) \
         == FT(n8) / FT(n10) );
 assert( FT( translate.hm(1,0) ) / FT( translate.hm(2,2) ) == FTzero );
 assert( FT( translate.hm(1,1) ) / FT( translate.hm(2,2) ) == FTone );
 assert( FT( translate.hm(1,2) ) / FT( translate.hm(2,2) ) \
         == FT(n1) / FT(n10) );
 assert( FT( translate.hm(2,0) ) / FT( translate.hm(2,2) ) == FTzero );
 assert( FT( translate.hm(2,1) ) / FT( translate.hm(2,2) ) == FTzero );
 assert( FT( translate.hm(2,2) ) / FT( translate.hm(2,2) ) == FTone );

 // rotation
 assert( FT( rot3.hm(0,0) ) / FT( rot3.hm(2,2) ) == f4o5 );
 assert( FT( rot3.hm(0,1) ) / FT( rot3.hm(2,2) ) == -f3o5 );
 assert( FT( rot3.hm(0,2) ) / FT( rot3.hm(2,2) ) == FTzero );
 assert( FT( rot3.hm(1,0) ) / FT( rot3.hm(2,2) ) == f3o5 );
 assert( FT( rot3.hm(1,1) ) / FT( rot3.hm(2,2) ) == f4o5 );
 assert( FT( rot3.hm(1,2) ) / FT( rot3.hm(2,2) ) == FTzero );
 assert( FT( rot3.hm(2,0) ) / FT( rot3.hm(2,2) ) == FTzero );
 assert( FT( rot3.hm(2,1) ) / FT( rot3.hm(2,2) ) == FTzero );
 assert( FT( rot3.hm(2,2) ) / FT( rot3.hm(2,2) ) == FTone );

 // scaling
 assert( FT( scale11.hm(0,0) ) / FT( scale11.hm(2,2) ) == fscale );
 assert( FT( scale11.hm(0,1) ) / FT( scale11.hm(2,2) ) == FTzero );
 assert( FT( scale11.hm(0,2) ) / FT( scale11.hm(2,2) ) == FTzero );
 assert( FT( scale11.hm(1,0) ) / FT( scale11.hm(2,2) ) == FTzero );
 assert( FT( scale11.hm(1,1) ) / FT( scale11.hm(2,2) ) == fscale );
 assert( FT( scale11.hm(1,2) ) / FT( scale11.hm(2,2) ) == FTzero );
 assert( FT( scale11.hm(2,0) ) / FT( scale11.hm(2,2) ) == FTzero );
 assert( FT( scale11.hm(2,1) ) / FT( scale11.hm(2,2) ) == FTzero );
 assert( FT( scale11.hm(2,2) ) / FT( scale11.hm(2,2) ) == FTone );

 // ident
 assert( FT( ident.hm(0,0) ) / FT( ident.hm(2,2) ) == FTone );
 assert( FT( ident.hm(0,1) ) / FT( ident.hm(2,2) ) == FTzero );
 assert( FT( ident.hm(0,2) ) / FT( ident.hm(2,2) ) == FTzero );
 assert( FT( ident.hm(1,0) ) / FT( ident.hm(2,2) ) == FTzero );
 assert( FT( ident.hm(1,1) ) / FT( ident.hm(2,2) ) == FTone );
 assert( FT( ident.hm(1,2) ) / FT( ident.hm(2,2) ) == FTzero );
 assert( FT( ident.hm(2,0) ) / FT( ident.hm(2,2) ) == FTzero );
 assert( FT( ident.hm(2,1) ) / FT( ident.hm(2,2) ) == FTzero );
 assert( FT( ident.hm(2,2) ) / FT( ident.hm(2,2) ) == FTone );

 // samples
 // cartesian == m
 assert( gat1.cartesian(1,2) == gat1.m(1,2) );
 assert( scale11.cartesian(0,0) == scale11.m(0,0) );
 assert( translate.cartesian(1,2) == translate.m(1,2) );
 assert( rot3.cartesian(1,0) == rot3.m(1,0) );
 assert( ident.cartesian(1,2) == ident.m(1,2) );
 assert( gscale.cartesian(1,1) == gscale.m(1,1) );

 // m == hm/hm(2,2)
 assert( gat3.m(1,2) == FT( gat3.hm(1,2)) / FT( gat3.hm(2,2) ));
 assert( gat4.m(1,0) == FT( gat4.hm(1,0)) / FT( gat4.hm(2,2) ));
 assert( gat5.m(2,0) == FT( gat5.hm(2,0)) / FT( gat5.hm(2,2) ));
 assert( translate.m(1,2) \
         == FT( translate.hm(1,2)) / FT( translate.hm(2,2) ));
 assert( scale11.m(1,1) == FT( scale11.hm(1,1)) / FT( scale11.hm(2,2) ));
 assert( rot2.m(2,0) == FT( rot2.hm(2,0)) / FT( rot2.hm(2,2) ));

 // homogeneous == hm
 assert( gat1.homogeneous(1,2) == gat1.hm(1,2) );
 assert( scale11.homogeneous(0,0) == scale11.hm(0,0) );
 assert( translate.homogeneous(1,2) == translate.hm(1,2) );
 assert( rot3.homogeneous(1,0) == rot3.hm(1,0) );
 assert( ident.homogeneous(1,2) == ident.hm(1,2) );
 assert( gscale.homogeneous(1,1) == gscale.hm(1,1) );

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_CLS_AFF_TRANSFORMATION_2_H
