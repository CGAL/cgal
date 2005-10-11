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
// file          : test_traits_2.C
// revision      : 
// revision_date : 
// author(s)     : Fernando Cacciola (fernando.cacciola@gmail.com)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/_test_types.h>
#include <CGAL/Kernel_traits.h>

//#define EXACT_KERNEL
//#define JUST_SINGLE_CASE
//#define CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE

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

#include <CGAL/_test_traits.C>

#ifdef EXACT_KERNEL
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt K ; 
#else
typedef CGAL::Exact_predicates_inexact_constructions_kernel K ; 
#endif

typedef K::RT      RT;
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
  
  void Set ( int x, int y ) { mP[(y*S)+x] = Point(mOX+mSX*RT(x),mOY+mSY*RT(y)) ; }

  static int idx( char c ) { return c - 'a' ; }
  
  RT mOX ;
  RT mOY ;
  RT mSX ;
  RT mSY ;
  Point mP[S*S] ;
} 
sGrid ;

struct Triedge
{
  Triedge( char const* desc ) 
  {
    mP[0]         = sGrid.at(desc[0]);
    mP[1] = mP[2] = sGrid.at(desc[1]);
    mP[3] = mP[4] = sGrid.at(desc[2]);
    mP[5]         = sGrid.at(desc[3]);
  }
  Triedge( char const* desc0, char const* desc1 ) 
  {
    mP[0]         = sGrid.at(desc0[0]);
    mP[1] = mP[2] = sGrid.at(desc0[1]);
    mP[3]         = sGrid.at(desc0[2]);
    mP[4]         = sGrid.at(desc1[0]);
    mP[5]         = sGrid.at(desc1[1]);
  }

  int idx( char const* d, int i ) { return d[i] - 'a' ; }
  
  typedef tuple<Point,Point> Edge ;
  typedef tuple<Edge,Edge,Edge> Edge_triple ;
  
  Edge_triple triple() const
  {
    return make_tuple( make_tuple(mP[0],mP[1]), make_tuple(mP[2],mP[3]), make_tuple(mP[4],mP[5]));
  }                                                         
   
  friend std::ostream& operator<<( std::ostream& os, Point const& aP )
  {
    return (os << '(' << aP.x() << ',' << aP.y() << ')' );
  }
  
  friend std::ostream& operator<<( std::ostream& os, Triedge const& aTriedge )
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
  
  const Triedge triedge[c] = { Triedge("fabg")
                             , Triedge("fbci")
                             , Triedge("abfa")
                             , Triedge("fbcg")
                             , Triedge("abch")
                             , Triedge("abci")
                             , Triedge("abcg")
                             , Triedge("afgl")
                             , Triedge("aght")
                             , Triedge("agc","qr")
                             , Triedge("acga")
                             } ;    
 
  const bool expected[c] = {true
                           ,true
                           ,true
                           ,true
                           
                           ,false
                           ,false
                           ,false
                           ,false
                           ,false
                           
                           ,true
                           ,true
                           };

  for(int i = 0 ; i < c ; ++ i )
    test_exist_event(i,sTraits,triedge[i],expected[i]) ;
}

void test_compare_events()
{
//    u v w x y 
//    p q r s t
//    k l m n o
//    f g h i j   
//    a b c d e

  const int c = 3 ;
    
  const Triedge triedgeA[c] = { Triedge("fabg")
                              , Triedge("upqv")
                              , Triedge("bgim")
                              } ;    
  
  const Triedge triedgeB[c] = { Triedge("fbci")
                              , Triedge("vqrw")
                              , Triedge("bgil")
                              } ;    
                             
  const CGAL::Comparison_result expected[c] = {CGAL::SMALLER
                                              ,CGAL::EQUAL
                                              ,CGAL::LARGER
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

  const int c = 3 ;
    
  const Point seed[c] = { sGrid.at('c')
                        , sGrid.at('g')
                        } ;    
  const Triedge triedgeA[c] = { Triedge("gcin")
                              , Triedge("agc","pt")
                              , Triedge("agc","kl")
                              } ;    
  
  const Triedge triedgeB[c] = { Triedge("gcij")
                              , Triedge("agc","no")
                              , Triedge("agc","mn")
                              } ;    
                             
  const CGAL::Comparison_result expected[c] = {CGAL::SMALLER
                                              ,CGAL::LARGER
                                              ,CGAL::EQUAL
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
    
  const Triedge triedgeA[c] = { Triedge("fach")
                              //, Triedge("upqv")
                              } ;    
  const Triedge triedgeB[c] = { Triedge("gcin")
                              //, Triedge("agc","pt")
                              } ;    
  
  const Triedge triedgeC[c] = { Triedge("gcij")
                              //, Triedge("agc","no")
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
    
  const Triedge triedgeA[c] = { Triedge("gmi","yu")
                              , Triedge("afb","yu")
                              , Triedge("gmi","tp")
                              , Triedge("afb","tp")
                              } ;    
  
  const Triedge triedgeB[c] = { Triedge("tyup")
                              , Triedge("tyup")
                              , Triedge("xsrw")
                              , Triedge("xsrw")
                              } ;    
                             
  const bool expected[c] = {true
                           ,false
                           ,true
                           ,false
                           };

  for(int i = 0 ; i < c ; ++ i )
    test_is_inside_offset_zone(i,sTraits,triedgeA[i],triedgeB[i],expected[i]) ;
}
   
int main()  
{ 
  std::cout << "Testing Straight_skeleton_traits_2\n";
  
#ifdef JUST_SINGLE_CASE
  sGrid.Setup(0,0,0.1,10); 
  test_is_inside_offset_zone(-1,sTraits,Triedge("gmi","yu"),Triedge("tyup"),true) ;
#else

  const int oc = 4 ;  
  const double ox[oc]={0,10,1000,10000};
  const double oy[oc]={0, 5,5000, 5000};
  
#ifdef EXACT_KERNEL 
  const int sc = 10 ;  
  const double sx[sc]={1,1e-10,1e-5,10,1e5,1e10,10 ,0.1,1e4 ,1e-4};
  const double sy[sc]={1,1e-10,1e-5,10,1e5,1e10,0.1,10 ,1e-4,1e4 };
#else
  const int sc = 10 ;  
  const double sx[sc]={1,ldexp(1.,-2),ldexp(1.,-8),ldexp(1.,-32),ldexp(1.,2),ldexp(1.,8),ldexp(1.,32),ldexp(1., 2),ldexp(1., 8),ldexp(1.,-8)};
  const double sy[sc]={1,ldexp(1.,-2),ldexp(1.,-8),ldexp(1.,-32),ldexp(1.,2),ldexp(1.,8),ldexp(1.,32),ldexp(1.,-2),ldexp(1.,-8),ldexp(1., 8)};
#endif

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
#endif

  if ( sReportSummary )   
  {
    std::cout << sSucceeded << " cases succeeded." << std::endl << sFailed << " cases failed." << std::endl ;
  }
  
  return sFailed == 0 ? 0 : 1 ;
}
 
