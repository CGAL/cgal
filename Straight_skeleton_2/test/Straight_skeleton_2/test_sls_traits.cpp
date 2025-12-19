// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#include <CGAL/test_sls_traits_types.h>


#ifdef CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE
#include <iostream>
#include <string>
void Straight_skeleton_traits_external_trace( std::string s )
{
  std::cout << s << std::endl ;
}
#endif

#include <CGAL/Straight_skeleton_builder_traits_2.h>

std::string sPrefix ;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef CGAL::Straight_skeleton_builder_traits_2<K> Traits ;

typedef Traits::FT        FT;
typedef Traits::Point_2   Point ;
typedef Traits::Segment_2 Segment ;

typedef Traits::Trisegment_2 Trisegment ;
typedef Trisegment::Self_ptr Trisegment_ptr ;

Traits sTraits ;

const int S = 5 ;

static int sid = 0;

struct Grid
{
  void Setup ( double ox, double oy, double sx, double sy )
  {
    mOX = ox;
    mOY = oy;
    mSX = sx;
    mSY = sy;

    for( int y = 0 ; y < S ; ++ y )
      for ( int x = 0 ; x < S ; ++ x )
        Set(x,y);
  }

  Point const& at( char c ) const { return mP[idx(c)] ; }

  static int idx( char c ) { return c - 'a' ; }

private :
  void Set ( int x, int y ) { mP[(y*S)+x] = Point(mOX+mSX*FT(x),mOY+mSY*FT(y)) ; }

  FT mOX ;
  FT mOY ;
  FT mSX ;
  FT mSY ;
  Point mP[S*S] ;
}
sGrid ;

using namespace CGAL::CGAL_SS_i ;

struct triple
{
  triple( char const* desc )
  {
    mP[0]         = sGrid.at(desc[0]);
    mP[1] = mP[2] = sGrid.at(desc[1]);
    mP[3] = mP[4] = sGrid.at(desc[2]);
    mP[5]         = sGrid.at(desc[3]);
  }
  triple( char const* desc0, char const* desc1 )
  {
    mP[0]         = sGrid.at(desc0[0]);
    mP[1] = mP[2] = sGrid.at(desc0[1]);
    mP[3]         = sGrid.at(desc0[2]);
    mP[4]         = sGrid.at(desc1[0]);
    mP[5]         = sGrid.at(desc1[1]);
  }

  int idx( char const* d, int i ) { return d[i] - 'a' ; }

  Trisegment_ptr trisegment() const
  {
    int sid0 = sid++;
    int sid1 = sid++;
    int sid2 = sid++;

    return sTraits.construct_ss_trisegment_2_object()(
             Segment(Point(mP[0].x(),mP[0].y()), Point(mP[1].x(),mP[1].y()), sid0),
             FT(1),
             Segment(Point(mP[2].x(),mP[2].y()), Point(mP[3].x(),mP[3].y()), sid1),
             FT(1),
             Segment(Point(mP[4].x(),mP[4].y()), Point(mP[5].x(),mP[5].y()), sid2),
             FT(1));
  }

  friend std::ostream& operator<<( std::ostream& os, Point const& aP )
  {
    return (os << '(' << aP.x() << ',' << aP.y() << ')' );
  }

  friend std::ostream& operator<<( std::ostream& os, triple const& aTriple )
  {
    return (os << "{[" << aTriple.mP[0] << "->" << aTriple.mP[1] << "]["
                       << aTriple.mP[2] << "->" << aTriple.mP[3] << "]["
                       << aTriple.mP[4] << "->" << aTriple.mP[5] << "}") ;
  }

  Point mP[6];
} ;

#include <CGAL/test_sls_traits_aux.cpp>

void test_exist_event()
{
//    u v w x y
//    p q r s t
//    k l m n o
//    f g h i j
//    a b c d e

  const int c = 10 ;

  const triple triples[c] = { triple("fabg")
                            , triple("fbci")
                            , triple("abfa")
                            , triple("fbcg")
                            , triple("abch")
                            , triple("abci")
                            , triple("abcg")
                            , triple("afgl")
                            , triple("agc","qr")
                            , triple("acga")
                            } ;

  const bool expected[c] = {true
                           ,true
                           ,true
                           ,true

                           ,true
                           ,true
                           ,true
                           ,false

                           ,true
                           ,true
                           };

  for(int i = 0 ; i < c ; ++ i )
    test_exist_event(i,sTraits,triples[i],expected[i]) ;
}

void test_compare_events()
{
//    u v w x y
//    p q r s t
//    k l m n o
//    f g h i j
//    a b c d e

  const int c = 4 ;

  const triple tripleA[c] = { triple("fabg")
                            , triple("upqv")
                            , triple("bgim")
                            , triple("abch")
                            } ;

  const triple tripleB[c] = { triple("fbci")
                            , triple("vqrw")
                            , triple("dmlh")
                            , triple("abci")
                            } ;

  const CGAL::Comparison_result expected[c] = {CGAL::SMALLER
                                              ,CGAL::EQUAL
                                              ,CGAL::LARGER
                                              ,CGAL::SMALLER
                                              };

  for(int i = 0 ; i < c ; ++ i )
    test_compare_events(i,sTraits,tripleA[i],tripleB[i],expected[i]) ;
}


typedef double* edge ;
typedef edge*   edgetriple ;



int main()
{
  using std::ldexp;
  std::cout << "Testing Straight_skeleton_traits_2\n";

  const int oc = 4 ;
  const double ox[oc]={0,10,1000,10000};
  const double oy[oc]={0, 5,5000, 5000};

  const int sc = 10 ;
  const double sx[sc]={1,ldexp(1.,-2),ldexp(1.,-8),ldexp(1.,-32),ldexp(1.,2),ldexp(1.,8),ldexp(1.,32),ldexp(1., 2),ldexp(1., 8),ldexp(1.,-8)};
  const double sy[sc]={1,ldexp(1.,-2),ldexp(1.,-8),ldexp(1.,-32),ldexp(1.,2),ldexp(1.,8),ldexp(1.,32),ldexp(1.,-2),ldexp(1.,-8),ldexp(1., 8)};

  for ( int si = 0; si < sc ; ++ si )
  {
    for ( int oi = 0 ; oi < oc ; ++ oi )
    {
      sPrefix = boost::str(boost::format("(o=%f,%f s=%f,%f) ") % ox[oi] % oy[oi] % sx[si] % sy[si]);
      sGrid.Setup(ox[oi],oy[oi],sx[si],sy[si]);
      test_exist_event   ();
      test_compare_events();
    }
  }

  if ( sReportSummary )
  {
    std::cout << sSucceeded << " cases succeeded." << std::endl << sFailed << " cases failed." << std::endl ;
  }

  return (sFailed == 0) ? EXIT_SUCCESS : EXIT_FAILURE ;
}
