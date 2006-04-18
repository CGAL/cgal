// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#include <CGAL/test_sls_traits_types.h>


#ifdef CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE
#include<iostream>
#include<string>
void Straight_skeleton_traits_external_trace( std::string s )
{
  std::cout << s << std::endl ;
}
#endif

#include <CGAL/Straight_skeleton_builder_traits_2.h>

std::string sPrefix ;

#include <CGAL/test_sls_traits_aux.C>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::FT      FT;
typedef K::Point_2 Point ;

typedef CGAL::Straight_skeleton_builder_traits_2<K> Traits ;

Traits sTraits ;


const int S = 5 ;

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

  private :

  void Set ( int x, int y ) { mP[(y*S)+x] = Point(mOX+mSX*FT(x),mOY+mSY*FT(y)) ; }

  static int idx( char c ) { return c - 'a' ; }

  FT mOX ;
  FT mOY ;
  FT mSX ;
  FT mSY ;
  Point mP[S*S] ;
}
sGrid ;

using namespace CGAL::CGAL_SS_i ;

struct triedge
{
  triedge( char const* desc )
  {
    mP[0]         = sGrid.at(desc[0]);
    mP[1] = mP[2] = sGrid.at(desc[1]);
    mP[3] = mP[4] = sGrid.at(desc[2]);
    mP[5]         = sGrid.at(desc[3]);
  }
  triedge( char const* desc0, char const* desc1 )
  {
    mP[0]         = sGrid.at(desc0[0]);
    mP[1] = mP[2] = sGrid.at(desc0[1]);
    mP[3]         = sGrid.at(desc0[2]);
    mP[4]         = sGrid.at(desc1[0]);
    mP[5]         = sGrid.at(desc1[1]);
  }

  int idx( char const* d, int i ) { return d[i] - 'a' ; }

  typedef CGAL::CGAL_SS_i::Vertex <FT> Vertex ;
  typedef CGAL::CGAL_SS_i::Edge   <FT> Edge ;
  typedef CGAL::CGAL_SS_i::Triedge<FT> Triedge ;

  Triedge triple() const
  {
    return Triedge( Edge( Vertex(mP[0].x(),mP[0].y()),  Vertex(mP[1].x(),mP[1].y()))
                  , Edge( Vertex(mP[2].x(),mP[2].y()),  Vertex(mP[3].x(),mP[3].y()))
                  , Edge( Vertex(mP[4].x(),mP[4].y()),  Vertex(mP[5].x(),mP[5].y()))
                  );
  }

  friend std::ostream& operator<<( std::ostream& os, Point const& aP )
  {
    return (os << '(' << aP.x() << ',' << aP.y() << ')' );
  }

  friend std::ostream& operator<<( std::ostream& os, triedge const& aTriedge )
  {
    return (os << "{[" << aTriedge.mP[0] << "->" << aTriedge.mP[1] << "]["
                       << aTriedge.mP[2] << "->" << aTriedge.mP[3] << "]["
                       << aTriedge.mP[4] << "->" << aTriedge.mP[5] << "}") ;
  }

  Point mP[6];
} ;


void test_exist_event()
{
//    u v w x y
//    p q r s t
//    k l m n o
//    f g h i j
//    a b c d e

  const int c = 11 ;

  const triedge triedges[c] = { triedge("fabg")
                              , triedge("fbci")
                              , triedge("abfa")
                              , triedge("fbcg")
                              , triedge("abch")
                              , triedge("abci")
                              , triedge("abcg")
                              , triedge("afgl")
                              , triedge("aght")
                              , triedge("agc","qr")
                              , triedge("acga")
                              } ;

  const bool expected[c] = {true
                           ,true
                           ,true
                           ,true

                           ,true
                           ,true
                           ,true
                           ,false
                           ,false

                           ,true
                           ,true
                           };

  for(int i = 0 ; i < c ; ++ i )
    test_exist_event(i,sTraits,triedges[i],expected[i]) ;
}

void test_compare_events()
{
//    u v w x y
//    p q r s t
//    k l m n o
//    f g h i j
//    a b c d e

  const int c = 4 ;

  const triedge triedgeA[c] = { triedge("fabg")
                              , triedge("upqv")
                              , triedge("bgim")
                              , triedge("abch")
                              } ;

  const triedge triedgeB[c] = { triedge("fbci")
                              , triedge("vqrw")
                              , triedge("dmlh")
                              , triedge("abci")
                              } ;

  const CGAL::Comparison_result expected[c] = {CGAL::SMALLER
                                              ,CGAL::EQUAL
                                              ,CGAL::LARGER
                                              ,CGAL::SMALLER
                                              };

  for(int i = 0 ; i < c ; ++ i )
    test_compare_events(i,sTraits,triedgeA[i],triedgeB[i],expected[i]) ;
}

void test_compare_sdist_to_seed1()
{
//    u v w x y
//    p q r s t
//    k l m n o
//    f g h i j
//    a b c d e

  const int c = 4 ;

  const Point seed[c] = { sGrid.at('c')
                        , sGrid.at('g')
                        , sGrid.at('g')
                        , sGrid.at('b')
                        } ;
  const triedge triedgeA[c] = { triedge("gcin")
                              , triedge("agc","tp")
                              , triedge("agc","lk")
                              , triedge("abcj")
                              } ;

  const triedge triedgeB[c] = { triedge("gcij")
                              , triedge("agc","on")
                              , triedge("agc","nm")
                              , triedge("abcg")
                              } ;

  const CGAL::Comparison_result expected[c] = {CGAL::SMALLER
                                              ,CGAL::LARGER
                                              ,CGAL::EQUAL
                                              ,CGAL::LARGER
                                              };

  for(int i = 0 ; i < c ; ++ i )
    test_compare_sdist_to_seed(i,sTraits,seed[i],triedgeA[i],triedgeB[i],expected[i]) ;
}

void test_compare_sdist_to_seed2()
{
//    u v w x y
//    p q r s t
//    k l m n o
//    f g h i j
//    a b c d e

  const int c = 1 ;

  const triedge triedgeA[c] = { triedge("fach")
                              //, triedge("upqv")
                              } ;
  const triedge triedgeB[c] = { triedge("gcin")
                              //, triedge("agc","pt")
                              } ;

  const triedge triedgeC[c] = { triedge("gcij")
                              //, triedge("agc","no")
                              } ;

  const CGAL::Comparison_result expected[c] = {CGAL::SMALLER
                                              //,CGAL::EQUAL
                                              };

  for(int i = 0 ; i < c ; ++ i )
    test_compare_sdist_to_seed(i,sTraits,triedgeA[i],triedgeB[i],triedgeC[i],expected[i]) ;
}


void test_is_inside_offset_zone()
{
//    u v w x y
//    p q r s t
//    k l m n o
//    f g h i j
//    a b c d e

  const int c = 4 ;

  const triedge triedgeA[c] = { triedge("gmi","yu")
                              , triedge("afb","yu")
                              , triedge("gmi","tp")
                              , triedge("afb","tp")
                              } ;

  const triedge triedgeB[c] = { triedge("tyup")
                              , triedge("tyup")
                              , triedge("xsrw")
                              , triedge("xsrw")
                              } ;

  const bool expected[c] = {true
                           ,false
                           ,true
                           ,false
                           };

  for(int i = 0 ; i < c ; ++ i )
    test_is_inside_offset_zone(i,sTraits,triedgeA[i],triedgeB[i],expected[i]) ;
}

typedef double* edge ;
typedef edge*   edgetriple ;



int main()
{
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
      test_compare_sdist_to_seed1();
      test_compare_sdist_to_seed2();
      test_is_inside_offset_zone();
    }
  }

  if ( sReportSummary )
  {
    std::cout << sSucceeded << " cases succeeded." << std::endl << sFailed << " cases failed." << std::endl ;
  }

  return sFailed == 0 ? 0 : 1 ;
}
