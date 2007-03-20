// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid Surface_mesh_simplification license may use this file in
// accordance with the Surface_mesh_simplification license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_IMPL_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_IMPL_H 1

CGAL_BEGIN_NAMESPACE


//
// Implementation of the strategy from:
//
//  "Fast and Memory Efficient Polygonal Symplification"
//  Peter Lindstrom, Greg Turk
//

namespace Surface_mesh_simplification
{

template<class ECM>
LindstromTurkCore<ECM>::LindstromTurkCore( Params const& aParams, Profile const& aProfile )
  :
   mParams(aParams)
  ,mProfile(aProfile)
{
  Extract_triangle_data();
  Extract_boundary_data();
}

template<class ECM>
void LindstromTurkCore<ECM>::Extract_boundary_data()
{
  for ( const_border_edge_iterator it = mProfile.border_edges().begin(), eit = mProfile.border_edges().end() ; it != eit ; ++ it )
  {
    edge_descriptor border_edge = *it ;
    
    edge_descriptor face_edge = opposite_edge(border_edge,surface()) ;
        
    vertex_descriptor sv = source(face_edge,surface());
    vertex_descriptor tv = target(face_edge,surface());
    
    Point const& sp = get_point(sv);
    Point const& tp = get_point(tv);
    
    Vector v = tp - sp ;
    Vector n = Point_cross_product(tp,sp) ;
    
    CGAL_ECMS_LT_TRACE(1,"Boundary edge. S:" << xyz_to_string(sp) << " T:" << xyz_to_string(tp)
                      << " V:" << xyz_to_string(v) << " N:" << xyz_to_string(n) 
                      ) ;
    
    mBdry_data.push_back( Boundary_data(sp,tp,v,n) ) ;
  }  
}

template<class ECM>
void LindstromTurkCore<ECM>::Extract_triangle_data()
{
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
    
    CGAL_ECMS_LT_TRACE(1,"Extracting triangle v" << tri.v0->id() << "->v" << tri.v1->id() << "->v" << tri.v2->id()
                      << " N:" << xyz_to_string(lNormalV) << " L:" << lNormalL
                      );
    
    mTriangle_data.push_back(Triangle_data(lNormalV,lNormalL));
  }
}

template<class ECM>
typename LindstromTurkCore<ECM>::Optional_point LindstromTurkCore<ECM>::compute_placement()
{
  Optional_point  rPlacementP ;
  Optional_vector lPlacementV ;
  
  CGAL_ECMS_LT_TRACE(0,"Computing LT data for E" << mProfile.v0v1()->id() << " (V" << mProfile.v0()->id() << "->V" << mProfile.v1()->id() << ")" );
  
  //
  // Each vertex constrian is an equation of the form: Ai * v = bi
  // Where 'v' is a vector representing the vertex,
  // 'Ai' is a (row) vector and 'bi' a scalar.
  //
  // The vertex is completely determined with 3 such constrians, 
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
  // The member variable mConstrinas contains A and b. Indidivual constrians (Ai,bi) can be added to it.
  // Once 3 such constrians have been added 'v' is directly solved a:
  //
  //  v = b*inverse(A)
  //
  // A constrian (Ai,bi) must be alpha-compatible with the previously added constrians (see Paper); if it's not, is discarded.
  //
  if ( mBdry_data.size() > 0 )
    Add_boundary_preservation_constrians(mBdry_data);
  
  if ( mConstrians.n < 3 )
    Add_volume_preservation_constrians(mTriangle_data);
    
  if ( mConstrians.n < 3 )
    Add_boundary_and_volume_optimization_constrians(mBdry_data,mTriangle_data); 
    
  if ( mConstrians.n < 3 )
    Add_shape_optimization_constrians(mProfile.link());
  
  // It might happen that there were not enough alpha-compatible constrians.
  // In that case there is simply no good vertex placement
  if ( mConstrians.n == 3 ) 
  {
    // If the matrix is singular it's inverse cannot be computed so an 'absent' value is returned.
    optional<Matrix> lOptional_Ai = inverse_matrix(mConstrians.A);
    if ( lOptional_Ai )
    {
      Matrix const& lAi = *lOptional_Ai ;
      
      CGAL_ECMS_LT_TRACE(2,"     b: " << xyz_to_string(mConstrians.b) );
      CGAL_ECMS_LT_TRACE(2,"inv(A): " << matrix_to_string(lAi) );
      
      lPlacementV = mConstrians.b * lAi ;
      
      CGAL_ECMS_LT_TRACE(0,"New vertex point: " << xyz_to_string(*lPlacementV) );
    }
    else
      CGAL_ECMS_LT_TRACE(1,"Can't solve optimization, singular system.");
  }
  else
    CGAL_ECMS_LT_TRACE(1,"Can't solve optimization, not enough alpha-compatible constrians.");
  
  if ( lPlacementV )
    rPlacementP = filter_infinity( ORIGIN + (*lPlacementV) );
    
  return rPlacementP;
}

template<class ECM>
typename LindstromTurkCore<ECM>::Optional_FT LindstromTurkCore<ECM>::compute_cost( Optional_point const& aP )
{
  Optional_FT rCost ;

  if ( aP )
  {
    Vector lV = (*aP) - ORIGIN ;

    FT lSquaredLength = squared_distance(mProfile.p0(),mProfile.p1());
   
    FT lBdryCost   = Compute_boundary_cost(lV ,mBdry_data);
    FT lVolumeCost = Compute_volume_cost  (lV ,mTriangle_data);
    FT lShapeCost  = Compute_shape_cost   (*aP,mProfile.link());
    
    FT lTotalCost =   FT(mParams.VolumeWeight)   * lVolumeCost
                    + FT(mParams.BoundaryWeight) * lBdryCost   * lSquaredLength
                    + FT(mParams.ShapeWeight)    * lShapeCost  * lSquaredLength * lSquaredLength ;
    
    rCost = filter_infinity(lTotalCost);
    
    CGAL_ECMS_LT_TRACE(0,"\nSquared edge length: " << lSquaredLength 
                      << "\nBoundary cost: " << lBdryCost << " weight: " << mParams.BoundaryWeight
                      << "\nVolume cost: " << lVolumeCost << " weight: " << mParams.VolumeWeight
                      << "\nShape cost: " << lShapeCost  << " weight: " << mParams.ShapeWeight
                      << "\nTOTAL COST: " << lTotalCost 
                      );
                      
  }
    
  return rCost;
}


template<class ECM>
void LindstromTurkCore<ECM>::Add_boundary_preservation_constrians( Boundary_data_vector const& aBdry )
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
      
      CGAL_ECMS_LT_TRACE(1,"vx:" << vx << " vy:" << vy << " vz:" << vz << " e1x:" << e1x << " e1y:" << e1y << " e1z:" << e1z );
    }     
    
    Matrix H = LT_product(e1);
    
    Vector c = cross_product(e1,e2);
    
    CGAL_ECMS_LT_TRACE(1
                      ,"Adding boundary preservation constrians. SumV:" << xyz_to_string(e1) 
                       << " SumN:" << xyz_to_string(e2) 
                       << "\nH:" << matrix_to_string(H)
                       << "\nc:" << xyz_to_string(c)                       
                      );
    
    mConstrians.Add_from_gradient(H,c);
  }
}

template<class ECM>
void LindstromTurkCore<ECM>::Add_volume_preservation_constrians( Triangle_data_vector const& aTriangles )
{
  CGAL_ECMS_LT_TRACE(1,"Adding volume preservation constrians. " << aTriangles.size() << " triangles.");
  
  Vector lSumV = NULL_VECTOR ;
  FT     lSumL(0) ;
  
  for( typename Triangle_data_vector::const_iterator it = aTriangles.begin(), eit = aTriangles.end() ; it != eit ; ++it )
  {
    lSumV = lSumV + it->NormalV ;
    lSumL = lSumL + it->NormalL ;  
  }   
  
  CGAL_ECMS_LT_TRACE(1, " SumV:" << xyz_to_string(lSumV) << " SumL:" << lSumL );
  
  mConstrians.Add_if_alpha_compatible(lSumV,lSumL);   

}

template<class ECM>
void LindstromTurkCore<ECM>::Add_boundary_and_volume_optimization_constrians( Boundary_data_vector const& aBdry
                                                                            , Triangle_data_vector const& aTriangles
                                                                            )
{
  CGAL_ECMS_LT_TRACE(1,"Adding boundary and volume optimization constrians. ");
  
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
  
  CGAL_ECMS_LT_TRACE(2,"Hv:" << matrix_to_string(H) << "\ncv:" << xyz_to_string(c) ) ;
  
  
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
    
    CGAL_ECMS_LT_TRACE(2,"Hb:" << matrix_to_string(Hb) << "\ncb:" << xyz_to_string(cb) ) ;
    
    //
    // Weighted average
    //
    FT lScaledBoundaryWeight = FT(9) * FT(mParams.BoundaryWeight) * squared_distance(mProfile.p0(),mProfile.p1())  ;
    
    H *= mParams.VolumeWeight ;
    c = c * mParams.VolumeWeight ;
    
    H += lScaledBoundaryWeight * Hb ;
    c = c + ( lScaledBoundaryWeight * cb ) ;
    
    CGAL_ECMS_LT_TRACE(2,"H:" << matrix_to_string(H) << "\nc:" << xyz_to_string(c) ) ;
    CGAL_ECMS_LT_TRACE(2,"VolW:" << mParams.VolumeWeight << " BdryW:" << mParams.BoundaryWeight << " ScaledBdryW:" << lScaledBoundaryWeight ) ;
    
  }
  
  mConstrians.Add_from_gradient(H,c);
}

template<class ECM>
void LindstromTurkCore<ECM>::Add_shape_optimization_constrians( vertex_descriptor_vector const& aLink )
{
  FT s((double)aLink.size());
  
  Matrix H (  s,0.0,0.0
           ,0.0,  s,0.0
           ,0.0,0.0,  s
           );
           
  Vector c = NULL_VECTOR ;
  
  for( typename vertex_descriptor_vector::const_iterator it = aLink.begin(), eit = aLink.end() ; it != eit ; ++it )
    c = c + (ORIGIN - get_point(*it)) ;  
           
  CGAL_ECMS_LT_TRACE(1,"Adding shape optimization constrians. Shape vector: " << xyz_to_string(c) );
  
  mConstrians.Add_from_gradient(H,c);
}

template<class ECM>
typename LindstromTurkCore<ECM>::FT
LindstromTurkCore<ECM>::Compute_boundary_cost( Vector const& v, Boundary_data_vector const& aBdry )
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

template<class ECM>
typename LindstromTurkCore<ECM>::FT
LindstromTurkCore<ECM>::Compute_volume_cost( Vector const& v, Triangle_data_vector const& aTriangles )
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

template<class ECM>
typename LindstromTurkCore<ECM>::FT
LindstromTurkCore<ECM>::Compute_shape_cost( Point const& p, vertex_descriptor_vector const& aLink )
{
  FT rCost(0);
  
  for( typename vertex_descriptor_vector::const_iterator it = aLink.begin(), eit = aLink.end() ; it != eit ; ++it )
    rCost += squared_distance(p,get_point(*it)) ;  
    
  return rCost ;
}

template<class ECM>
void LindstromTurkCore<ECM>::Constrians::Add_if_alpha_compatible( Vector const& Ai, FT const& bi )
{
  FT slai = Ai*Ai ;
  
  CGAL_ECMS_LT_TRACE(3,"[constrians] Adding new if alpha-compatble.\nslai: " << slai );
  
  FT l = CGAL_NTS sqrt(slai) ;
  
  CGAL_ECMS_LT_TRACE(3,"[constrians] Adding new if alpha-compatble.\nslai: " << slai );
  
  Vector Ain ;
  FT     bin ;
  
  if ( !CGAL_NTS is_zero(l) )
  {
    Ain = Ai / l ;
    bin = bi / l ;
  }
  else
  {
    CGAL_ECMS_LT_TRACE(3,"[constrians] l is ZERO." );
    
    Ain = Vector( big_value(), big_value(), big_value() ) ;
    bin = big_value() ;
  }
  
  CGAL_ECMS_LT_TRACE(3,"[constrians] Ain: " << xyz_to_string(Ain) << " bin:" << bin );
  
  bool lAddIt = true ;
  
  if ( n == 1 )
  {
    FT d01 = A.r0() * Ai  ;
    
    FT sla0 = A.r0() * A.r0() ;
    FT sd01 = d01 * d01 ;

    FT max = sla0 * slai * squared_cos_alpha() ;     
    
    CGAL_ECMS_LT_TRACE(3,"[constrians] Second constrain. d01: " << d01 << " sla0:" << sla0 << " sd01:" << sd01 << " max:" << max );
    
    if ( sd01 > max )
      lAddIt = false ;
  }
  else if ( n == 2 )
  {
    Vector N = cross_product(A.r0(),A.r1());
    
    FT dc012 = N * Ai ;
    
    FT slc01  = N * N ;
    FT sdc012 = dc012 * dc012;       

    FT min = slc01 * slai * squared_sin_alpha() ;
    
    CGAL_ECMS_LT_TRACE(3,"[constrians] Thirds constrain. N: " << xyz_to_string(N) << " dc012:" << dc012 << " slc01:" << slc01 
                       << " sdc012:" << sdc012 << " min:" << min );
    
    if ( sdc012 <= min )
      lAddIt = false ;
  }
  
  if ( lAddIt )
  {
    switch ( n )
    {
      case 0 :
        A.r0() = Ain ;
        b = Vector(bin,b.y(),b.z());
        break ;
      case 1 :
        A.r1() = Ain ;
        b = Vector(b.x(),bin,b.z());
        break ;
      case 2 :
        A.r2() = Ain ;
        b = Vector(b.x(),b.y(),bin);
        break ;
    }
    
    CGAL_ECMS_LT_TRACE(3,"[constrains] Accepting # " << n << " A:" << matrix_to_string(A) << " b:" << xyz_to_string(b) ) ;
    
    ++ n ;
    
  }
  else
  {
    CGAL_ECMS_LT_TRACE(3,"[constrains] INCOMPATIBLE. Discarded" ) ;
  }
}

template<class V>
int index_of_max_component ( V const& v )
{
  typedef typename Kernel_traits<V>::Kernel::FT FT ;
  
  int i = 0 ;
  
  FT max = v.x();
  FT y = v.y();

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

template<class ECM>
void LindstromTurkCore<ECM>::Constrians::Add_from_gradient ( Matrix const& H, Vector const& c )
{
  CGAL_ECMS_LT_TRACE(3,"[constrains] Adding from gradient. Current n=" << n ) ;
  
  CGAL_precondition(n >= 0 && n<=2 );
  
  switch(n)
  {
    case 0 :
    
      Add_if_alpha_compatible(H.r0(),-c.x());
      Add_if_alpha_compatible(H.r1(),-c.y());
      Add_if_alpha_compatible(H.r2(),-c.z());
      
      break;  
      
    case 1 :
      {
        Vector const& A0 = A.r0();
        
        CGAL_assertion( A0 != NULL_VECTOR ) ;
        
        Vector AbsA0( CGAL_NTS abs(A0.x())
                    , CGAL_NTS abs(A0.y())
                    , CGAL_NTS abs(A0.z())
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
        
        CGAL_ECMS_LT_TRACE(3,"[constrains] Q0:" << xyz_to_string(Q0) ) ;
        
        CGAL_assertion( Q0 != NULL_VECTOR ) ;
        
        Vector Q1 = cross_product(A0,Q0);
    
        Vector A1 = H * Q0 ;
        Vector A2 = H * Q1 ;
        FT b1 = - ( Q0 * c ) ;
        FT b2 = - ( Q1 * c ) ;
        
        CGAL_ECMS_LT_TRACE(3,"[constrains] Q1:" << xyz_to_string(Q1) << " A1: " << xyz_to_string(A1) << " A2:" << xyz_to_string(A2) << "\nb1:" << b1 << " b2:" << b2 ) ;
        
        Add_if_alpha_compatible(A1,b1);
        Add_if_alpha_compatible(A2,b2);
        
      }
      break ;
    
    case 2:
      {
    
        Vector Q = cross_product(A.r0(),A.r1());
        
        Vector A2 = H * Q ;
        
        FT b2 = - ( Q * c ) ;
        
        CGAL_ECMS_LT_TRACE(3,"[constrains] Q:" << xyz_to_string(Q) << " A2: " << xyz_to_string(A2) << " b2:" << b2 ) ;
        Add_if_alpha_compatible(A2,b2);
        
      }
      break ;
      
  }
}

} // namespace Surface_mesh_simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROMTURK_CORE_IMPL_H //
// EOF //
 

