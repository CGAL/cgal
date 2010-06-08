// Copyright (c) 2010 Fernando Luis Cacciola Carballal. All rights reserved.
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Straight_skeleton_2/test/Straight_skeleton_2/test_offset_builder.cpp $
// $Id: test_offset_builder.cpp 33476 2006-08-22 08:24:07Z afabri $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#include<iostream>
#include<algorithm>
#include<utility>

#include<boost/optional.hpp>

#include <CGAL/assertions_behaviour.h>

#include <CGAL/test_sls_types.h>

#include <CGAL/create_straight_skeleton_2.h>

void test_polygon_vertices_range()
{
  IPoint Pts[3]= { IPoint(1,1), IPoint(2,2), IPoint(3,3) } ;
  
  IPolygon VectorPoly(Pts, Pts+3) ;
  
  IPolygonPtr VectorPolyPtr( new IPolygon(Pts, Pts+3) ) ;
  
  CGAL::Polygon_2<IK> CGALPoly(Pts,Pts+3) ;
  
  assert( equal( CGAL_ss_i::vertices_begin(VectorPoly)   , CGAL_ss_i::vertices_end(VectorPoly)   , Pts ) ) ;
  assert( equal( CGAL_ss_i::vertices_begin(VectorPolyPtr), CGAL_ss_i::vertices_end(VectorPolyPtr), Pts ) ) ;
  assert( equal( CGAL_ss_i::vertices_begin(CGALPoly)     , CGAL_ss_i::vertices_end(CGALPoly)     , Pts ) ) ;
}

struct Test
{
  Test()
  {
    IPoint lOPts[4]= { IPoint(0,0), IPoint(10,0), IPoint(10,10), IPoint(0,10) } ;
    IPoint lHPts[4]= { IPoint(3,3), IPoint(3,6), IPoint(6,6), IPoint(6,3) } ;
    
    double lWts[4]= { 1, 1, 1, 1 } ;
    
    mOuterP = IPolygonPtr( new IPolygon(lOPts, lOPts+4) ) ;
    mOuterW = IWeightsPtr( new IWeights(lWts , lWts +4) ) ;
    
    mHolesP.push_back( IPolygonPtr( new IPolygon(lHPts, lHPts+4) ) ) ;
    mHolesW.push_back( IWeightsPtr( new IWeights(lWts , lWts +4) ) ) ;
    
    mWeight = 1.0 ;
  }
  
  void test_create_stright_skeleton_2_A()
  {
    mSls = CGAL::create_straight_skeleton_2(mOuterP->begin()
                                           ,mOuterP->end()
                                           ,mOuterW->begin()
                                           ,mOuterW->end()
                                           ,mHolesP. begin()
                                           ,mHolesP. end()
                                           ,mHolesW. begin()
                                           ,mHolesW. end()
                                           ,mMaxT
                                           ,mK
                                           );
                                                  
    mSls = CGAL::create_straight_skeleton_2(mOuterP->begin()
                                           ,mOuterP->end()
                                           ,mOuterW->begin()
                                           ,mOuterW->end()
                                           ,mHolesP. begin()
                                           ,mHolesP. end()
                                           ,mHolesW. begin()
                                           ,mHolesW. end()
                                           ,mMaxT
                                           );
                                           
    mSls = CGAL::create_straight_skeleton_2(mOuterP->begin()
                                           ,mOuterP->end()
                                           ,mOuterW->begin()
                                           ,mOuterW->end()
                                           ,mHolesP. begin()
                                           ,mHolesP. end()
                                           ,mHolesW. begin()
                                           ,mHolesW. end()
                                           );
  }
  
  void test_create_stright_skeleton_2_B()
  {
    mSls = CGAL::create_straight_skeleton_2(mOuterP->begin()
                                           ,mOuterP->end()
                                           ,mHolesP. begin()
                                           ,mHolesP. end()
                                           ,mWeight
                                           ,mMaxT
                                           ,mK
                                           );
                                                   
    mSls = CGAL::create_straight_skeleton_2(mOuterP->begin()
                                           ,mOuterP->end()
                                           ,mHolesP. begin()
                                           ,mHolesP. end()
                                           ,mWeight
                                           ,mMaxT
                                           );
                                           
    mSls = CGAL::create_straight_skeleton_2(mOuterP->begin()
                                           ,mOuterP->end()
                                           ,mHolesP. begin()
                                           ,mHolesP. end()
                                           ,mWeight
                                           );
  }
  
  void test_create_stright_skeleton_2_C()
  {
    mSls = CGAL::create_straight_skeleton_2(mOuterP->begin()
                                           ,mOuterP->end()
                                           ,mOuterW->begin()
                                           ,mOuterW->end()
                                           ,mMaxT
                                           ,mK
                                           
                                           );
    mSls = CGAL::create_straight_skeleton_2(mOuterP->begin()
                                           ,mOuterP->end()
                                           ,mOuterW->begin()
                                           ,mOuterW->end()
                                           ,mMaxT
                                           );
                                                   
    mSls = CGAL::create_straight_skeleton_2(mOuterP->begin()
                                           ,mOuterP->end()
                                           ,mOuterW->begin()
                                           ,mOuterW->end()
                                           );
  }
  
  void test_create_stright_skeleton_2_D()
  {
    mSls = CGAL::create_straight_skeleton_2(mOuterP->begin()
                                           ,mOuterP->end()
                                           ,mWeight
                                           ,mMaxT
                                           ,mK
                                           
                                           );
    mSls = CGAL::create_straight_skeleton_2(mOuterP->begin()
                                           ,mOuterP->end()
                                           ,mWeight
                                           ,mMaxT
                                           );
                                                   
    mSls = CGAL::create_straight_skeleton_2(mOuterP->begin()
                                           ,mOuterP->end()
                                           ,mWeight
                                           );
  }
  
  void test_create_stright_skeleton_2_E()
  {
    mSls = CGAL::create_straight_skeleton_2(mOuterP
                                           ,mWeight
                                           ,mMaxT
                                           ,mK
                                           
                                           );
    mSls = CGAL::create_straight_skeleton_2(mOuterP
                                           ,mWeight
                                           ,mMaxT
                                           );
                                                   
    mSls = CGAL::create_straight_skeleton_2(mOuterP
                                           ,mWeight
                                           );
                                           
    mSls = CGAL::create_straight_skeleton_2(mOuterP);
  }
  
  IPolygonPtr             mOuterP ;
  IWeightsPtr             mOuterW ;
  IBoundaries             mHolesP ;
  IWeightsPtrContainer    mHolesW ;
  IK                      mK ;
  boost::optional<double> mMaxT ;
  double                  mWeight ; 
  ISlsPtr                 mSls ;
  
} ;

int main( int argc, char const* argv[] )
{
  test_polygon_vertices_range(); 
  
  Test t ;
  
  t.test_create_stright_skeleton_2_A();
  t.test_create_stright_skeleton_2_B();
  t.test_create_stright_skeleton_2_C();
  t.test_create_stright_skeleton_2_D();
  t.test_create_stright_skeleton_2_E();
  
  return 0 ;
}


