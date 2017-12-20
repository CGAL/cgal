// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid Surface_mesh_simplification license may use this file in
// accordance with the Surface_mesh_simplification license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_IMPL_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_IMPL_H 1

#include <CGAL/license/Surface_mesh_simplification.h>


namespace CGAL {


//
// Implementation of the strategy from:
//
//  "Fast and Memory Efficient Polygonal Symplification"
//  Peter Lindstrom, Greg Turk
//

namespace Surface_mesh_simplification
{

template<class ECM, class K>
LindstromTurkCore<ECM,K>::LindstromTurkCore( Params const& aParams, Profile const& aProfile )
  :
   mParams(aParams)
  ,mProfile(aProfile)
  ,mConstraints_n(0)
  ,mConstraints_A(NULL_MATRIX)
  ,mConstraints_b(NULL_VECTOR) 
{
  double alpha = 1.0 * CGAL_PI / 180.0;
 
  FT cos_alpha = std::cos(alpha);
  FT sin_alpha = std::sin(alpha);

  mSquared_cos_alpha = cos_alpha * cos_alpha ;
  mSquared_sin_alpha = sin_alpha * sin_alpha ;

  Extract_triangle_data();
  Extract_boundary_data();
}

template<class ECM, class K>
void LindstromTurkCore<ECM,K>::Extract_boundary_data()
{
  mBdry_data.reserve(mProfile.border_edges().size());
  for ( const_border_edge_iterator it = mProfile.border_edges().begin(), eit = mProfile.border_edges().end() ; it != eit ; ++ it )
  {
    halfedge_descriptor border_edge = *it ;
    
    halfedge_descriptor face_edge = opposite(border_edge,surface()) ;
        
    vertex_descriptor sv = source(face_edge,surface());
    vertex_descriptor tv = target(face_edge,surface());
    
    Point const& sp = get_point(sv);
    Point const& tp = get_point(tv);
    
    Vector v = tp - sp ;
    Vector n = Point_cross_product(tp,sp) ;
    
    CGAL_ECMS_LT_TRACE(1,"  Boundary edge. S:" << xyz_to_string(sp) << " T:" << xyz_to_string(tp)
                      << " V:" << xyz_to_string(v) << " N:" << xyz_to_string(n) 
                      ) ;
    
    mBdry_data.push_back( Boundary_data(sp,tp,v,n) ) ;
  }  
}

template<class ECM, class K>
void LindstromTurkCore<ECM,K>::Extract_triangle_data()
{
  mTriangle_data.reserve(mProfile.triangles().size());
  for ( const_triangle_iterator it = mProfile.triangles().begin(), eit = mProfile.triangles().end() ; it != eit ; ++ it )
  {
    Triangle const& tri = *it ;

    Point const& p0 = get_point(tri.v0);
    Point const& p1 = get_point(tri.v1);
    Point const& p2 = get_point(tri.v2);
    
    Vector v01 = p1 - p0 ;
    Vector v02 = p2 - p0 ;
    
    Vector lNormalV = cross_product(v01,v02);
    
    FT lNormalL = Point_cross_product(p0,p1) * (p2-ORIGIN);
    
    CGAL_ECMS_LT_TRACE(1,"  Extracting triangle v" << tri.v0 << "->v" << tri.v1 << "->v" << tri.v2
                      << " N:" << xyz_to_string(lNormalV) << " L:" << n_to_string(lNormalL)
                      );
    
    mTriangle_data.push_back(Triangle_data(lNormalV,lNormalL));
  }
}

template<class ECM, class K>
typename LindstromTurkCore<ECM,K>::Optional_point LindstromTurkCore<ECM,K>::compute_placement()
{
  Optional_point  rPlacementP ;
  Optional_vector lPlacementV ;
  
  CGAL_ECMS_LT_TRACE(0,"Computing LT placement for E" << mProfile.v0_v1() << " (V" << mProfile.v0() << "->V" << mProfile.v1() << ")" );
  
  //
  // Each vertex constraint is an equation of the form: Ai * v = bi
  // Where 'v' is a vector representing the vertex,
  // 'Ai' is a (row) vector and 'bi' a scalar.
  //
  // The vertex is completely determined with 3 such constraints, 
  // so is the solution to the folloing system:
  //
  //  A.r0(). * v = b0
  //  A1 * v = b1
  //  A2 * v = b2
  //
  // Which in matrix form is :  A * v = b
  //
  // (with 'A' a 3x3 matrix and 'b' a vector)
  //
  // The member variable mConstrinas contains A and b. Indidivual constraints (Ai,bi) can be added to it.
  // Once 3 such constraints have been added 'v' is directly solved a:
  //
  //  v = b*inverse(A)
  //
  // A constraint (Ai,bi) must be alpha-compatible with the previously added constraints (see Paper); if it's not, is discarded.
  //
  if ( mBdry_data.size() > 0 )
    Add_boundary_preservation_constraints(mBdry_data);
  
  if ( mConstraints_n < 3 )
    Add_volume_preservation_constraints(mTriangle_data);
    
  if ( mConstraints_n < 3 )
    Add_boundary_and_volume_optimization_constraints(mBdry_data,mTriangle_data); 
    
  if ( mConstraints_n < 3 )
    Add_shape_optimization_constraints(mProfile.link());
  
  // It might happen that there were not enough alpha-compatible constraints.
  // In that case there is simply no good vertex placement
  if ( mConstraints_n == 3 ) 
  {
    // If the matrix is singular it's inverse cannot be computed so an 'absent' value is returned.
    optional<Matrix> lOptional_Ai = inverse_matrix(mConstraints_A);
    if ( lOptional_Ai )
    {
      Matrix const& lAi = *lOptional_Ai ;
      
      CGAL_ECMS_LT_TRACE(2,"       b: " << xyz_to_string(mConstraints_b) );
      CGAL_ECMS_LT_TRACE(2,"  inv(A): " << matrix_to_string(lAi) );
      
      lPlacementV = filter_infinity(mConstraints_b * lAi) ;
      
      CGAL_ECMS_LT_TRACE(0,"  New vertex point: " << xyz_to_string(*lPlacementV) );
    }
    else
    {
      CGAL_ECMS_LT_TRACE(1,"  Can't solve optimization, singular system.");
    }  
  }
  else
  {
    CGAL_ECMS_LT_TRACE(1,"  Can't solve optimization, not enough alpha-compatible constraints.");
  }  
  
  if ( lPlacementV )
    rPlacementP = ORIGIN + (*lPlacementV) ;
    
  return rPlacementP;
}

template<class ECM, class K>
typename LindstromTurkCore<ECM,K>::Optional_FT LindstromTurkCore<ECM,K>::compute_cost( Optional_point const& aP )
{
  Optional_FT rCost ;

  if ( aP )
  {
    CGAL_ECMS_LT_TRACE(0,"Computing LT cost for E" << mProfile.v0_v1() );
    Vector lV = (*aP) - ORIGIN ;

    FT lSquaredLength = squared_distance(mProfile.p0(),mProfile.p1());
   
    FT lBdryCost   = Compute_boundary_cost(lV ,mBdry_data);
    FT lVolumeCost = Compute_volume_cost  (lV ,mTriangle_data);
    FT lShapeCost  = Compute_shape_cost   (*aP,mProfile.link());
    
    FT lTotalCost =   FT(mParams.VolumeWeight)   * lVolumeCost
                    + FT(mParams.BoundaryWeight) * lBdryCost   * lSquaredLength
                    + FT(mParams.ShapeWeight)    * lShapeCost  * lSquaredLength * lSquaredLength ;
    
    rCost = filter_infinity(lTotalCost);
    
    CGAL_ECMS_LT_TRACE(0,  "    Squared edge length: " << n_to_string(lSquaredLength)
                      << "\n    Boundary cost: " << n_to_string(lBdryCost) << " weight: " << mParams.BoundaryWeight
                      << "\n    Volume cost: " << n_to_string(lVolumeCost) << " weight: " << mParams.VolumeWeight
                      << "\n    Shape cost: " << n_to_string(lShapeCost)   << " weight: " << mParams.ShapeWeight
                      << "\n  TOTAL COST: " << n_to_string(lTotalCost) 
                      );
                      
  }
    
  return rCost;
}


template<class ECM, class K>
void LindstromTurkCore<ECM,K>::Add_boundary_preservation_constraints( Boundary_data_vector const& aBdry )
{
 
  if ( aBdry.size() > 0 )
  {
    Vector e1 = NULL_VECTOR ; 
    Vector e2 = NULL_VECTOR ;

    FT e1x = FT(0.0), e1y = FT(0.0), e1z = FT(0.0) ;
    
    for ( typename Boundary_data_vector::const_iterator it = aBdry.begin() ; it != aBdry.end() ; ++ it )
    {
      e1 = e1 + it->v ;
      e2 = e2 + it->n ;
      
      FT vx = it->v.x() ;
      FT vy = it->v.y() ;
      FT vz = it->v.z() ;
      
      e1x = e1x + vx;
      e1y = e1y + vy;
      e1z = e1z + vz;
      
      CGAL_ECMS_LT_TRACE(1,"    vx:" << n_to_string(vx) << " vy:" << n_to_string(vy) << " vz:" << n_to_string(vz) << " e1x:" 
                            << n_to_string(e1x) << " e1y:" << n_to_string(e1y) << " e1z:" << n_to_string(e1z) );
    }     
    
    Matrix H = LT_product(e1);
    
    Vector c = cross_product(e1,e2);
    
    CGAL_ECMS_LT_TRACE(1
                      ,"  Adding boundary preservation constraint."
                       << "\n      SumV:" << xyz_to_string(e1) 
                       << "\n      SumN:" << xyz_to_string(e2) 
                       << "\n      H:" << matrix_to_string(H)
                       << "\n      c:" << xyz_to_string(c)                       
                      );
    
    Add_constraint_from_gradient(H,c);
  }
}

template<class ECM, class K>
void LindstromTurkCore<ECM,K>::Add_volume_preservation_constraints( Triangle_data_vector const& aTriangles )
{
  CGAL_ECMS_LT_TRACE(1,"  Adding volume preservation constraint. " << aTriangles.size() << " triangles.");
  
  Vector lSumV = NULL_VECTOR ;
  FT     lSumL(0) ;
  
  for( typename Triangle_data_vector::const_iterator it = aTriangles.begin(), eit = aTriangles.end() ; it != eit ; ++it )
  {
    lSumV = lSumV + it->NormalV ;
    lSumL = lSumL + it->NormalL ;  
  }   
  
  CGAL_ECMS_LT_TRACE(1, "      SumV:" << xyz_to_string(lSumV) << "\n      SumL:" << n_to_string(lSumL) );
  
  Add_constraint_if_alpha_compatible(lSumV,lSumL);   

}

template<class ECM, class K>
void LindstromTurkCore<ECM,K>::Add_boundary_and_volume_optimization_constraints( Boundary_data_vector const& aBdry
                                                                            , Triangle_data_vector const& aTriangles
                                                                            )
{
  CGAL_ECMS_LT_TRACE(1,"  Adding boundary and volume optimization constraints. ");
  
  Matrix H = NULL_MATRIX ;
  Vector c = NULL_VECTOR ;

  //
  // Volume optimization  
  //
  for( typename Triangle_data_vector::const_iterator it = aTriangles.begin(), eit = aTriangles.end() ; it != eit ; ++it )
  {
    Triangle_data const& lTri = *it ;
    
    H += direct_product(lTri.NormalV,lTri.NormalV) ;
    
    c = c - ( lTri.NormalL * lTri.NormalV ) ;
  }   
  
  CGAL_ECMS_LT_TRACE(2,"      Hv:" << matrix_to_string(H) << "\n      cv:" << xyz_to_string(c) ) ;
  
  
  if ( aBdry.size() > 0 )
  {
    //
    // Boundary optimization
    //
    Matrix Hb = NULL_MATRIX ;
    Vector cb = NULL_VECTOR ;
    
    for ( typename Boundary_data_vector::const_iterator it = aBdry.begin() ; it != aBdry.end() ; ++ it )
    {
      Matrix H = LT_product(it->v);
      Vector c = cross_product(it->v,it->n);
      
      Hb += H ;
      cb = cb + c ;
    }     
    
    CGAL_ECMS_LT_TRACE(2,"      Hb:" << matrix_to_string(Hb) << "\n      cb:" << xyz_to_string(cb) ) ;
    
    //
    // Weighted average
    //
    FT lScaledBoundaryWeight = FT(9) * FT(mParams.BoundaryWeight) * squared_distance(mProfile.p0(),mProfile.p1())  ;
    
    H *= mParams.VolumeWeight ;
    c = c * mParams.VolumeWeight ;
    
    H += lScaledBoundaryWeight * Hb ;
    c = c + ( lScaledBoundaryWeight * cb ) ;
    
    CGAL_ECMS_LT_TRACE(2,"      H:" << matrix_to_string(H) << "\n      c:" << xyz_to_string(c) ) ;
    CGAL_ECMS_LT_TRACE(2,"      VolW:" << mParams.VolumeWeight << " BdryW:" << mParams.BoundaryWeight << " ScaledBdryW:" << lScaledBoundaryWeight ) ;
    
  }
  
  Add_constraint_from_gradient(H,c);
}

template<class ECM, class K>
void LindstromTurkCore<ECM,K>::Add_shape_optimization_constraints( vertex_descriptor_vector const& aLink )
{
  FT s((double)aLink.size());
  
  Matrix H (  s,0.0,0.0
           ,0.0,  s,0.0
           ,0.0,0.0,  s
           );
           
  Vector c = NULL_VECTOR ;
  
  for( typename vertex_descriptor_vector::const_iterator it = aLink.begin(), eit = aLink.end() ; it != eit ; ++it )
    c = c + (ORIGIN - get_point(*it)) ;  
           
  CGAL_ECMS_LT_TRACE(1,"  Adding shape optimization constraint. Shape vector: " << xyz_to_string(c) );
  
  Add_constraint_from_gradient(H,c);
}

template<class ECM, class K>
typename LindstromTurkCore<ECM,K>::FT
LindstromTurkCore<ECM,K>::Compute_boundary_cost( Vector const& v, Boundary_data_vector const& aBdry )
{
  FT rCost(0);
  for ( typename Boundary_data_vector::const_iterator it = aBdry.begin() ; it != aBdry.end() ; ++ it )
  {
    Vector u = (it->t - ORIGIN ) - v ;
    Vector c = cross_product(it->v,u);
    rCost += c*c;  
  }     
  return rCost / FT(4) ;
}

template<class ECM, class K>
typename LindstromTurkCore<ECM,K>::FT
LindstromTurkCore<ECM,K>::Compute_volume_cost( Vector const& v, Triangle_data_vector const& aTriangles )
{
  FT rCost(0);
  
  for( typename Triangle_data_vector::const_iterator it = aTriangles.begin(), eit = aTriangles.end() ; it != eit ; ++it )
  {
    Triangle_data const& lTri = *it ;
    
    FT lF = lTri.NormalV * v - lTri.NormalL ;
    
    rCost += ( lF * lF ) ;
    
  }
    
  return rCost / FT(36) ;
}

template<class ECM, class K>
typename LindstromTurkCore<ECM,K>::FT
LindstromTurkCore<ECM,K>::Compute_shape_cost( Point const& p, vertex_descriptor_vector const& aLink )
{
  FT rCost(0);
  
  for( typename vertex_descriptor_vector::const_iterator it = aLink.begin(), eit = aLink.end() ; it != eit ; ++it )
    rCost += squared_distance(p,get_point(*it)) ;  
    
  return rCost ;
}

template<class ECM, class K>
void LindstromTurkCore<ECM,K>::Add_constraint_if_alpha_compatible( Vector const& Ai, FT const& bi )
{
  CGAL_ECMS_LT_TRACE(3,"    Adding new constraints if alpha-compatible.\n      Ai: " << xyz_to_string(Ai) << "\n      bi:" << n_to_string(bi) << ")" );

  FT slai = Ai*Ai ;
  
  CGAL_ECMS_LT_TRACE(3,"\n      slai: " << n_to_string(slai) << ")" );
  
  if ( is_finite(slai) )
  {
    FT l = CGAL_NTS sqrt( slai ) ;
  
    CGAL_ECMS_LT_TRACE(3,"      l: " << n_to_string(l) );
    
    if ( !CGAL_NTS is_zero(l) )
    {
      Vector Ain = Ai / l ;
      FT     bin = bi / l ;
      
      CGAL_ECMS_LT_TRACE(3,"      Ain: " << xyz_to_string(Ain) << " bin:" << n_to_string(bin) );
      
      bool lAddIt = true ;
      
      if ( mConstraints_n == 1 )
      {
        FT d01 = mConstraints_A.r0() * Ai  ;
        
        FT sla0 = mConstraints_A.r0() * mConstraints_A.r0() ;
        FT sd01 = d01 * d01 ;
    
        FT max = sla0 * slai * mSquared_cos_alpha ;     
        
        CGAL_ECMS_LT_TRACE(3,"      Second constraint. d01: " << n_to_string(d01) << " sla0:" << n_to_string(sla0) << " sd01:" << n_to_string(sd01) << " max:" << n_to_string(max) );
        
        if ( sd01 > max )
          lAddIt = false ;
      }
      else if ( mConstraints_n == 2 )
      {
        Vector N = cross_product(mConstraints_A.r0(),mConstraints_A.r1());
        
        FT dc012 = N * Ai ;
        
        FT slc01  = N * N ;
        FT sdc012 = dc012 * dc012;       
    
        FT min = slc01 * slai * mSquared_sin_alpha ;
        
        CGAL_ECMS_LT_TRACE(3,"      Third constraint. N: " << xyz_to_string(N) << " dc012:" << n_to_string(dc012) << " slc01:" << n_to_string(slc01)
                          << " sdc012:" << n_to_string(sdc012) << " min:" << n_to_string(min) );
        
        if ( sdc012 <= min )
          lAddIt = false ;
      }
      
      if ( lAddIt )
      {
        switch ( mConstraints_n )
        {
          case 0 :
            mConstraints_A.r0() = Ain ;
            mConstraints_b = Vector(bin,mConstraints_b.y(),mConstraints_b.z());
            break ;
          case 1 :
            mConstraints_A.r1() = Ain ;
            mConstraints_b = Vector(mConstraints_b.x(),bin,mConstraints_b.z());
            break ;
          case 2 :
            mConstraints_A.r2() = Ain ;
            mConstraints_b = Vector(mConstraints_b.x(),mConstraints_b.y(),bin);
            break ;
        }
        
        CGAL_ECMS_LT_TRACE(3,"      Accepting # " << mConstraints_n << " A:" << matrix_to_string(mConstraints_A) << " b:" << xyz_to_string(mConstraints_b) ) ;
        
        ++ mConstraints_n ;
        
      }
      else
      {
        CGAL_ECMS_LT_TRACE(3,"      INCOMPATIBLE. Discarded" ) ;
      }
    }
    else
    {
      CGAL_ECMS_LT_TRACE(3,"        l is ZERO. Discarded" );
    }
  }
  else
  {
    CGAL_ECMS_LT_TRACE(3,"      OVERFLOW. Discarded" ) ;
  }
}

template<class V>
int index_of_max_component ( V const& v )
{
  typedef typename Kernel_traits<V>::Kernel::FT FT ;
  
  int i = 0 ;
  
  FT max = v.x();

  if ( max < v.y() )
  {
    max = v.y();
    i = 1 ;
  }  
  if ( max < v.z() )
  {
    max = v.z();
    i = 2 ;
  }  
  return i ;
}

template<class ECM, class K>
void LindstromTurkCore<ECM,K>::Add_constraint_from_gradient ( Matrix const& H, Vector const& c )
{
  CGAL_ECMS_LT_TRACE(3,"    Adding constraint from gradient. Current n=" << mConstraints_n ) ;
  
  CGAL_precondition(mConstraints_n >= 0 && mConstraints_n<=2 );
  
  switch(mConstraints_n)
  {
    case 0 :
    
      Add_constraint_if_alpha_compatible(H.r0(),-c.x());
      Add_constraint_if_alpha_compatible(H.r1(),-c.y());
      Add_constraint_if_alpha_compatible(H.r2(),-c.z());
      
      break;  
      
    case 1 :
      {
        Vector const& A0 = mConstraints_A.r0();
        
        CGAL_assertion( A0 != NULL_VECTOR ) ;
        
        Vector AbsA0( CGAL::abs(A0.x())
                      , CGAL::abs(A0.y())
                      , CGAL::abs(A0.z())
                    );
           
        Vector Q0;      
        switch ( index_of_max_component(AbsA0) ) 
        {
          // Since A0 is guaranteed to be non-zero, the denominators here are known to be non-zero too.
          
          case 0: Q0 = Vector(- A0.z()/A0.x(),0              ,1              ); break;
          case 1: Q0 = Vector(0              ,- A0.z()/A0.y(),1              ); break;
          case 2: Q0 = Vector(1              ,0              ,- A0.x()/A0.z()); break;

           default : Q0 = NULL_VECTOR ; // This should never happen!
        }
        
        CGAL_ECMS_LT_TRACE(3,"      Q0:" << xyz_to_string(Q0) ) ;
        
        CGAL_assertion( Q0 != NULL_VECTOR ) ;
        
        Vector Q1 = cross_product(A0,Q0);
    
        Vector A1 = H * Q0 ;
        Vector A2 = H * Q1 ;
        
        FT b1 = - ( Q0 * c ) ;
        FT b2 = - ( Q1 * c ) ;
        
        CGAL_ECMS_LT_TRACE(3,"      Q1:" << xyz_to_string(Q1) << "\n      A1: " << xyz_to_string(A1) << "\n      A2:" << xyz_to_string(A2) 
                          << "\n      b1:" << n_to_string(b1) << "\n      b2:" << n_to_string(b2) ) ;
        
        Add_constraint_if_alpha_compatible(A1,b1);
        Add_constraint_if_alpha_compatible(A2,b2);
        
      }
      break ;
    
    case 2:
      {
    
        Vector Q = cross_product(mConstraints_A.r0(),mConstraints_A.r1());
        
        Vector A2 = H * Q ;
        
        FT b2 = - ( Q * c ) ;
        
        CGAL_ECMS_LT_TRACE(3,"      Q:" << xyz_to_string(Q) << "\n      A2: " << xyz_to_string(A2) << "\n      b2:" << n_to_string(b2) ) ;

        Add_constraint_if_alpha_compatible(A2,b2);
        
      }
      break ;
      
  }
}

} // namespace Surface_mesh_simplification

} //namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROMTURK_CORE_IMPL_H //
// EOF //
 
