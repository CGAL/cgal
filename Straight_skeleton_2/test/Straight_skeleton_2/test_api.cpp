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

#include <CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>

void test_polygon_vertices_range()
{
  IPoint Pts[3]= { IPoint(1,1), IPoint(2,2), IPoint(3,3) } ;
  
  IPolygon VectorPoly(Pts, Pts+3) ;
  
  IPolygonPtr VectorPolyPtr( new IPolygon(Pts, Pts+3) ) ;
  
  CGAL::Polygon_2<IK> CGALPoly(Pts,Pts+3) ;
  
  assert( equal( CGAL::CGAL_SS_i::vertices_begin(VectorPoly)   , CGAL::CGAL_SS_i::vertices_end(VectorPoly)   , Pts ) ) ;
  assert( equal( CGAL::CGAL_SS_i::vertices_begin(VectorPolyPtr), CGAL::CGAL_SS_i::vertices_end(VectorPolyPtr), Pts ) ) ;
  assert( equal( CGAL::CGAL_SS_i::vertices_begin(CGALPoly)     , CGAL::CGAL_SS_i::vertices_end(CGALPoly)     , Pts ) ) ;
}

struct Test
{
  typedef CGAL::Polygon_2<IK>            Poly ;
  typedef boost::shared_ptr<Poly>        PolyPtr ;
  typedef std::vector<PolyPtr>           PolyPtrContainer ;
  typedef CGAL::Polygon_with_holes_2<IK> PolyWH ;
  typedef boost::shared_ptr<PolyWH>      PolyWHPtr ;
  typedef std::vector<PolyWHPtr>         PolyWHPtrContainer ;
  
  Test()
  {
    IPoint lOPts[4]= { IPoint(0,0), IPoint(10,0), IPoint(10,10), IPoint(0,10) } ;
    IPoint lHPts[4]= { IPoint(3,3), IPoint(3,6), IPoint(6,6), IPoint(6,3) } ;
    
    double lWts[4]= { 1, 1, 1, 1 } ;
  
    mOuterP = Poly(lOPts, lOPts+4) ;
    mOuterW = IWeightsPtr( new IWeights(lWts , lWts +4) ) ;
    
    mHolesP.push_back( PolyPtr    ( new Poly    (lHPts, lHPts+4) ) ) ;
    mHolesW.push_back( IWeightsPtr( new IWeights(lWts , lWts +4) ) ) ;
    
    mWeight         = 1.0 ;
    mOffsetDistance = 2.0 ;
  }
  
  void test_create_straight_skeleton_2_A()
  {
    mSls = CGAL::create_straight_skeleton_2(mOuterP. vertices_begin()
                                           ,mOuterP. vertices_end()
                                           ,mOuterW->begin()
                                           ,mOuterW->end()
                                           ,mHolesP. begin()
                                           ,mHolesP. end()
                                           ,mHolesW. begin()
                                           ,mHolesW. end()
                                           ,mMaxT
                                           ,mK
                                           );
                                                  
    mSls = CGAL::create_straight_skeleton_2(mOuterP. vertices_begin()
                                           ,mOuterP. vertices_end()
                                           ,mOuterW->begin()
                                           ,mOuterW->end()
                                           ,mHolesP. begin()
                                           ,mHolesP. end()
                                           ,mHolesW. begin()
                                           ,mHolesW. end()
                                           ,mMaxT
                                           );
                                           
    mSls = CGAL::create_straight_skeleton_2(mOuterP. vertices_begin()
                                           ,mOuterP. vertices_end()
                                           ,mOuterW->begin()
                                           ,mOuterW->end()
                                           ,mHolesP. begin()
                                           ,mHolesP. end()
                                           ,mHolesW. begin()
                                           ,mHolesW. end()
                                           );
  }
  
  void test_create_straight_skeleton_2_B()
  {
    mSls = CGAL::create_straight_skeleton_2(mOuterP. vertices_begin()
                                           ,mOuterP. vertices_end()
                                           ,mHolesP. begin()
                                           ,mHolesP. end()
                                           ,mWeight
                                           ,mMaxT
                                           ,mK
                                           );
                                                   
    mSls = CGAL::create_straight_skeleton_2(mOuterP. vertices_begin()
                                           ,mOuterP. vertices_end()
                                           ,mHolesP. begin()
                                           ,mHolesP. end()
                                           ,mWeight
                                           ,mMaxT
                                           );
                                           
    mSls = CGAL::create_straight_skeleton_2(mOuterP. vertices_begin()
                                           ,mOuterP. vertices_end()
                                           ,mHolesP. begin()
                                           ,mHolesP. end()
                                           ,mWeight
                                           );
  }
  
  void test_create_straight_skeleton_2_C()
  {
    mSls = CGAL::create_straight_skeleton_2(mOuterP. vertices_begin()
                                           ,mOuterP. vertices_end()
                                           ,mOuterW->begin()
                                           ,mOuterW->end()
                                           ,mMaxT
                                           ,mK
                                           
                                           );
    mSls = CGAL::create_straight_skeleton_2(mOuterP. vertices_begin()
                                           ,mOuterP. vertices_end()
                                           ,mOuterW->begin()
                                           ,mOuterW->end()
                                           ,mMaxT
                                           );
                                                   
    mSls = CGAL::create_straight_skeleton_2(mOuterP. vertices_begin()
                                           ,mOuterP. vertices_end()
                                           ,mOuterW->begin()
                                           ,mOuterW->end()
                                           );
  }
  
  void test_create_straight_skeleton_2_D()
  {
    mSls = CGAL::create_straight_skeleton_2(mOuterP. vertices_begin()
                                           ,mOuterP. vertices_end()
                                           ,mWeight
                                           ,mMaxT
                                           ,mK
                                           
                                           );
    mSls = CGAL::create_straight_skeleton_2(mOuterP. vertices_begin()
                                           ,mOuterP. vertices_end()
                                           ,mWeight
                                           ,mMaxT
                                           );
                                                   
    mSls = CGAL::create_straight_skeleton_2(mOuterP. vertices_begin()
                                           ,mOuterP. vertices_end()
                                           ,mWeight
                                           );
  }
  
  void test_create_straight_skeleton_2_E()
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
  
  void test_create_straight_skeleton_2_PWH()
  {
    mSls = CGAL::create_straight_skeleton_2(mPolyWH, mWeight, mMaxT, mK);
    mSls = CGAL::create_straight_skeleton_2(mPolyWH, mWeight, mMaxT);
  }
  
  void test_create_interior_straight_skeleton_2_A()
  {
    mSls = CGAL::create_interior_straight_skeleton_2(mOuterP. vertices_begin()
                                                    ,mOuterP. vertices_end()
                                                    ,mHolesP. begin()
                                                    ,mHolesP. end()
                                                    ,mK
                                                    );
                                                    
    mSls = CGAL::create_interior_straight_skeleton_2(mOuterP. vertices_begin()
                                                    ,mOuterP. vertices_end()
                                                    ,mHolesP. begin()
                                                    ,mHolesP. end()
                                                   );
  }
  
  void test_create_interior_straight_skeleton_2_B()
  {
    mSls = CGAL::create_interior_straight_skeleton_2(mOuterP. vertices_begin()
                                                    ,mOuterP. vertices_end()
                                                    ,mK
                                                    );
                                                    
    mSls = CGAL::create_interior_straight_skeleton_2(mOuterP. vertices_begin()
                                                    ,mOuterP. vertices_end()
                                                   );
  }
  
  void test_create_interior_straight_skeleton_2_C()
  {
    mSls = CGAL::create_interior_straight_skeleton_2(mOuterP
                                                    ,mK
                                                    );
                                                    
    mSls = CGAL::create_interior_straight_skeleton_2(mOuterP
                                                    );
  }
  
  void test_create_interior_straight_skeleton_2_PWH()
  {
    mSls = CGAL::create_interior_straight_skeleton_2(mPolyWH, mMaxT);
    mSls = CGAL::create_interior_straight_skeleton_2(mPolyWH);
  }
  
  void test_create_exterior_straight_skeleton_2_A()
  {
    mSls = CGAL::create_exterior_straight_skeleton_2(mOuterP. vertices_begin()
                                                    ,mOuterP. vertices_end()
                                                    ,mHolesP. begin()
                                                    ,mHolesP. end()
                                                    ,mK
                                                    );
                                                    
    mSls = CGAL::create_exterior_straight_skeleton_2(mOuterP. vertices_begin()
                                                    ,mOuterP. vertices_end()
                                                    ,mHolesP. begin()
                                                    ,mHolesP. end()
                                                   );
  }
  
  void test_create_exterior_straight_skeleton_2_B()
  {
    mSls = CGAL::create_exterior_straight_skeleton_2(mOuterP. vertices_begin()
                                                    ,mOuterP. vertices_end()
                                                    ,mK
                                                    );
                                                    
    mSls = CGAL::create_exterior_straight_skeleton_2(mOuterP. vertices_begin()
                                                    ,mOuterP. vertices_end()
                                                   );
  }
  
  void test_create_exterior_straight_skeleton_2_C()
  {
    mSls = CGAL::create_exterior_straight_skeleton_2(mOuterP
                                                    ,mK
                                                    );
                                                    
    mSls = CGAL::create_exterior_straight_skeleton_2(mOuterP
                                                    );
  }
  
  void test_create_exterior_straight_skeleton_2_PWH()
  {
    mSls = CGAL::create_exterior_straight_skeleton_2(mPolyWH, mMaxT);
    mSls = CGAL::create_exterior_straight_skeleton_2(mPolyWH);
  }
  
  void test_create_offset_polygons_2_A()
  {
    mOffsetP = CGAL::create_offset_polygons_2<Poly>(mOffsetDistance, *mSls, mK);
    
    mOffsetP = CGAL::create_offset_polygons_2(mOffsetDistance, *mSls, mK);
  }
  
  void test_create_offset_polygons_2_B()
  {
    mOffsetP = CGAL::create_offset_polygons_2<Poly>(mOffsetDistance, *mSls);
    
    mOffsetP = CGAL::create_offset_polygons_2(mOffsetDistance, *mSls);
  }
  
  void test_create_interior_skeleton_and_offset_polygons_2_A()
  {
    mOffsetP = CGAL::create_interior_skeleton_and_offset_polygons_2(mOffsetDistance
                                                                   ,mOuterP
                                                                   ,mHolesP. begin()
                                                                   ,mHolesP. end()
                                                                   ,mK
                                                                   );
                                                                   
    mOffsetP = CGAL::create_interior_skeleton_and_offset_polygons_2(mOffsetDistance
                                                                   ,mOuterP
                                                                   ,mHolesP. begin()
                                                                   ,mHolesP. end()
                                                                   );
  }
  
  void test_create_interior_skeleton_and_offset_polygons_2_B()
  {
    mOffsetP = CGAL::create_interior_skeleton_and_offset_polygons_2(mOffsetDistance
                                                                   ,mOuterP
                                                                   ,mK
                                                                   );
                                                                   
    mOffsetP = CGAL::create_interior_skeleton_and_offset_polygons_2(mOffsetDistance
                                                                   ,mOuterP
                                                                   );
  }
  
  void test_create_exterior_skeleton_and_offset_polygons_2_A()
  {
    mOffsetP = CGAL::create_exterior_skeleton_and_offset_polygons_2(mOffsetDistance
                                                                   ,mOuterP
                                                                   ,mHolesP. begin()
                                                                   ,mHolesP. end()
                                                                   ,mK
                                                                   );
                                                                   
    mOffsetP = CGAL::create_exterior_skeleton_and_offset_polygons_2(mOffsetDistance
                                                                   ,mOuterP
                                                                   ,mHolesP. begin()
                                                                   ,mHolesP. end()
                                                                   );
  }
  
  void test_create_exterior_skeleton_and_offset_polygons_2_B()
  {
    mOffsetP = CGAL::create_exterior_skeleton_and_offset_polygons_2(mOffsetDistance
                                                                   ,mOuterP
                                                                   ,mK
                                                                   );
                                                                   
    mOffsetP = CGAL::create_exterior_skeleton_and_offset_polygons_2(mOffsetDistance
                                                                   ,mOuterP
                                                                   );
  }
  
  void test_create_skeleton_and_offset_polygons_with_holes_2()
  {
    mOffsetPWH = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(mOffsetDistance, mPolyWH );
                                                                   
    mOffsetPWH = CGAL::create_exterior_skeleton_and_offset_polygons_with_holes_2(mOffsetDistance, mPolyWH);
  }

  Poly                    mOuterP ;
  IWeightsPtr             mOuterW ;
  PolyPtrContainer        mHolesP ;
  IWeightsPtrContainer    mHolesW ;
  PolyWH                  mPolyWH ;
  IK                      mK ;
  boost::optional<double> mMaxT ;
  double                  mWeight ; 
  ISlsPtr                 mSls ;
  PolyPtrContainer        mOffsetP ;
  PolyWHPtrContainer      mOffsetPWH ; 
  double                  mOffsetDistance ;
  
} ;

int main( int argc, char const* argv[] )
{
  test_polygon_vertices_range(); 
  
  Test t ;
  
  t.test_create_straight_skeleton_2_A();
  t.test_create_straight_skeleton_2_B();
  t.test_create_straight_skeleton_2_C();
  t.test_create_straight_skeleton_2_D();
  t.test_create_straight_skeleton_2_E();
  t.test_create_straight_skeleton_2_PWH();
  
  t.test_create_interior_straight_skeleton_2_A();
  t.test_create_interior_straight_skeleton_2_B();
  t.test_create_interior_straight_skeleton_2_C();
  t.test_create_interior_straight_skeleton_2_PWH();
  
  t.test_create_exterior_straight_skeleton_2_A();
  t.test_create_exterior_straight_skeleton_2_B();
  t.test_create_exterior_straight_skeleton_2_C();
  t.test_create_exterior_straight_skeleton_2_PWH();
  
  t.test_create_offset_polygons_2_A();
  t.test_create_offset_polygons_2_B();
  
  t.test_create_interior_skeleton_and_offset_polygons_2_A();
  t.test_create_interior_skeleton_and_offset_polygons_2_B();
  
  t.test_create_exterior_skeleton_and_offset_polygons_2_A();
  t.test_create_exterior_skeleton_and_offset_polygons_2_B();
  
  t.test_create_skeleton_and_offset_polygons_with_holes_2();
  
  return 0 ;
}


