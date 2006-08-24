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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_IMPL_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_IMPL_H 1

CGAL_BEGIN_NAMESPACE

//
// Implementation of the strategy from:
//
//  "Fast and Memory Efficient Polygonal Symplification"
//  Peter Lindstrom, Greg Turk
//

namespace Triangulated_surface_mesh { namespace Simplification { namespace Edge_collapse  
{

template<class TSM>
LindstromTurkCore<TSM>::LindstromTurkCore( Params const&          aParams
                                         , edge_descriptor const& aP_Q
                                         , TSM&                   aSurface 
                                         , bool                   aComputeCost
                                         )
  :
   mParams(aParams)
  ,mP_Q(aP_Q)
  ,mSurface(aSurface)    
  ,mComputeCost(aComputeCost)
  
  ,mP  ( source       (aP_Q,aSurface) )
  ,mQ  ( target       (aP_Q,aSurface) )
  ,mQ_P( opposite_edge(aP_Q,aSurface) )
{
}

template<class TSM>
typename LindstromTurkCore<TSM>::result_type LindstromTurkCore<TSM>::compute()
{
  Optional_FT    lOptionalCost ;
  Optional_point lOptionalP ;
  
  CGAL_TSMS_LT_TRACE(2,"Computing LT data for E" << mP_Q->ID << " (V" << mP->ID << "->V" << mQ->ID << ")" );
  
  // Volume preservation and optimization constrians are based on the normals to the triangles in the star of the collapsing egde 
  // Triangle shape optimization constrians are based on the link of the collapsing edge (the cycle of vertices around the edge)
  Triangles lTriangles;
  Link      lLink;
  
  lTriangles.reserve(16);
  lLink     .reserve(16);
  
  Extract_triangles_and_link(lTriangles,lLink);

#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_LT_TRACE
  std::ostringstream ss ; 
  for( typename Link::const_iterator it = lLink.begin(), eit = lLink.end() ; it != eit ; ++it )
    ss << "v" << (*it)->ID << " " ;
  std::string s = ss.str(); 
  CGAL_TSMS_LT_TRACE(3,"Link: " << s );
#endif
  
  BoundaryEdges lBdry ;
  
  Extract_boundary_edges(lBdry);
    
  Point const& lP = get_point(mP) ;
  Point const& lQ = get_point(mQ) ;
  
  Optional_vector lOptionalV ;
  
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
  if ( lBdry.size() > 0 )
    Add_boundary_preservation_constrians(lBdry);
  
  if ( mConstrians.n < 3 )
    Add_volume_preservation_constrians(lTriangles);
    
  if ( mConstrians.n < 3 )
    Add_boundary_and_volume_optimization_constrians(lBdry,lTriangles); 
    
  if ( mConstrians.n < 3 )
    Add_shape_optimization_constrians(lLink);
  
  // It might happen that there were not enough alpha-compatible constrians.
  // In that case there is simply no good vertex placement
  if ( mConstrians.n == 3 ) 
  {
    // If the matrix is singular it's inverse cannot be computed so an 'absent' value is returned.
    optional<Matrix> lOptional_Ai = inverse_matrix(mConstrians.A);
    if ( lOptional_Ai )
    {
      Matrix const& lAi = *lOptional_Ai ;
      
      lOptionalV = mConstrians.b * lAi ;
      
      CGAL_TSMS_LT_TRACE(1,"New vertex point: " << xyz_to_string(*lOptionalV) );
    }
    else
      CGAL_TSMS_LT_TRACE(1,"Can't solve optimization, singular system.");
  }
  else
    CGAL_TSMS_LT_TRACE(1,"Can't solve optimization, not enough alpha-compatible constrians.");
  
  if ( lOptionalV )
  {
    lOptionalP = Optional_point( ORIGIN + (*lOptionalV) );
    
    if ( mComputeCost )
    {
      FT lSquaredLength = squared_distance(lP,lQ);
      
      CGAL_TSMS_LT_TRACE(1,"Squared edge length: " << lSquaredLength ) ;
      
      FT lBdryCost   = Compute_boundary_cost(*lOptionalV,lBdry);
      FT lVolumeCost = Compute_volume_cost  (*lOptionalV,lTriangles);
      FT lShapeCost  = Compute_shape_cost   (*lOptionalP,lLink);
      
      FT lTotalCost =   FT(mParams.VolumeWeight)   * lVolumeCost
                      + FT(mParams.BoundaryWeight) * lBdryCost   * lSquaredLength
                      + FT(mParams.ShapeWeight)    * lShapeCost  * lSquaredLength * lSquaredLength ;
      
      lOptionalCost = Optional_FT(lTotalCost);
      
      CGAL_TSMS_LT_TRACE(1,"\nSquared edge length: " << lSquaredLength 
                        << "\nBoundary cost: " << lBdryCost 
                        << "\nVolume cost: " << lVolumeCost 
                        << "\nShape cost: " << lShapeCost 
                        << "\nTOTAL COST: " << lTotalCost 
                        );
    }
  }
    
  return make_tuple(lOptionalCost,lOptionalP);
}


//
// Caches the "local boundary", that is, the sequence of 3 border edges: o->p, p->q, q->e 
//
template<class TSM>
void LindstromTurkCore<TSM>::Extract_boundary_edge( edge_descriptor edge, BoundaryEdges& rBdry )
{
  edge_descriptor face_edge = is_border(edge) ? opposite_edge(edge,mSurface) : edge ;
      
  vertex_descriptor sv = source(face_edge,mSurface);
  vertex_descriptor tv = target(face_edge,mSurface);
  
  Point const& sp = get_point(sv);
  Point const& tp = get_point(tv);
  
  Vector v = tp - sp ;
  Vector n = Point_cross_product(tp,sp) ;
  
  CGAL_TSMS_LT_TRACE(3,"Boundary edge. S:" << xyz_to_string(sp) << " T:" << xyz_to_string(tp)
                    << " V:" << xyz_to_string(v) << " N:" << xyz_to_string(n) 
                    ) ;
  
  rBdry.push_back( BoundaryEdge(sp,tp,v,n) ) ;
        
}

template<class TSM>
void LindstromTurkCore<TSM>::Extract_boundary_edges( vertex_descriptor const& v
                                                  , edge_descriptor_vector&  rCollected
                                                  , BoundaryEdges&           rBdry
                                                  )
{
  in_edge_iterator eb, ee ; 
  for ( tie(eb,ee) = in_edges(v,mSurface) ; eb != ee ; ++ eb )
  {
    edge_descriptor edge = *eb ;
    
    if ( is_undirected_edge_a_border(edge) && std::find(rCollected.begin(),rCollected.end(),edge) == rCollected.end() )
    {
      Extract_boundary_edge(edge,rBdry);
      rCollected.push_back(edge);
      rCollected.push_back(opposite_edge(edge,mSurface));
    }  
  }
}

template<class TSM>
void LindstromTurkCore<TSM>::Extract_boundary_edges( BoundaryEdges& rBdry )
{
  edge_descriptor_vector lCollected ;
  Extract_boundary_edges(mP,lCollected,rBdry);
  Extract_boundary_edges(mQ,lCollected,rBdry);
}

//
// Calculates the normal of the triangle (v0,v1,v2) (both vector and its length as (v0xv1).v2)
//
template<class TSM>
typename LindstromTurkCore<TSM>::Triangle LindstromTurkCore<TSM>::Get_triangle( vertex_descriptor const& v0
                                                                            , vertex_descriptor const& v1
                                                                            , vertex_descriptor const& v2 
                                                                            )
{
  Point const& p0 = get_point(v0);
  Point const& p1 = get_point(v1);
  Point const& p2 = get_point(v2);
  
  Vector v01 = p1 - p0 ;
  Vector v02 = p2 - p0 ;
  
  Vector lNormalV = cross_product(v01,v02);
  
  FT lNormalL = Point_cross_product(p0,p1) * (p2-ORIGIN);
  
  CGAL_TSMS_LT_TRACE(3,"Extracting triangle v" << v0->ID << "->v" << v1->ID << "->v" << v2->ID 
                    << " N:" << xyz_to_string(lNormalV) << " L:" << lNormalL
                    );
  
  return Triangle(lNormalV,lNormalL);
}                              

//
// If (v0,v1,v2) is a finite triangular facet of the mesh, that is, NONE of these vertices are boundary vertices,
// the triangle (properly oriented) is added to rTriangles.
// The triangle is encoded as its normal, calculated using the actual facet orientation [(v0,v1,v2) or (v0,v2,v1)]
//
template<class TSM>
void LindstromTurkCore<TSM>::Extract_triangle( vertex_descriptor const& v0
                                            , vertex_descriptor const& v1
                                            , vertex_descriptor const& v2 
                                            , edge_descriptor   const& e02
                                            , Triangles&               rTriangles
                                            )
{
  // The 3 vertices are obtained by circulating ccw around v0, that is, e02 = next_ccw(e01).
  // Since these vertices are NOT obtained by circulating the face, the actual triangle orientation is unspecified.
  
  // The triangle is oriented v0->v2->v1 if the next edge that follows e02 (which is the edge v0->v2) is v2->v1.
  if ( target(next_edge(e02,mSurface),mSurface) == v1 ) 
  {
    // The triangle is oriented v0->v2->v1.
    // In this case e02 is an edge of the facet.
    // If this facet edge is a border edge then this triangle is not in the mesh .
    if ( !is_border(e02) )
      rTriangles.push_back(Get_triangle(v0,v2,v1) ) ;
  }
  else
  {
    // The triangle is oriented v0->v1->v2.
    // In this case, e20 and not e02, is an edge of the facet.
    // If this facet edge is a border edge then this triangle is not in the mesh .
    if ( !is_border(opposite_edge(e02,mSurface)) )
      rTriangles.push_back(Get_triangle(v0,v1,v2) ) ;
  }
}

//
// Extract all triangles (its normals) and vertices (the link) around the collpasing edge p_q
//
template<class TSM>
void LindstromTurkCore<TSM>::Extract_triangles_and_link( Triangles& rTriangles, Link& rLink )
{
  // 
  // Extract around mP, CCW
  //  
  vertex_descriptor v0 = mP;
  vertex_descriptor v1 = mQ;
  
  edge_descriptor e02 = mP_Q;
  
  do
  {
    e02 = next_edge_ccw(e02,mSurface);
    
    vertex_descriptor v2 = target(e02,mSurface);
  
    if ( v2 != mQ )
    {
      CGAL_expensive_assertion ( std::find(rLink.begin(),rLink.end(),v2) == rLink.end() ) ;
      rLink.push_back(v2) ;
    }
      
    Extract_triangle(v0,v1,v2,e02,rTriangles);
    
    v1 = v2 ;
  }
  while ( e02 != mP_Q ) ;
  
  // 
  // Extract around mQ, CCW
  //  
  
  v0 = mQ;
  
  e02 = next_edge_ccw(mQ_P,mSurface);
  
  v1 = target(e02,mSurface); 
  
  // This could have been added to the link while circulating around mP
  if ( v1 != mP && std::find(rLink.begin(),rLink.end(),v1) == rLink.end() )
    rLink.push_back(v1) ;
  
  e02 = next_edge_ccw(e02,mSurface);
  
  do
  {
    vertex_descriptor v2 = target(e02,mSurface);

    // Any of the vertices found around mP can be reached again around mQ, but we can't duplicate them here.
    if ( v2 != mP && std::find(rLink.begin(),rLink.end(),v2) == rLink.end() )
      rLink.push_back(v2) ;
    
    Extract_triangle(v0,v1,v2,e02,rTriangles);
    
    v1 = v2 ;
     
    e02 = next_edge_ccw(e02,mSurface);
    
  }
  while ( e02 != mQ_P ) ;
}

template<class TSM>
void LindstromTurkCore<TSM>::Add_boundary_preservation_constrians( BoundaryEdges const& aBdry )
{
  
  if ( aBdry.size() > 0 )
  {
    Vector e1 = NULL_VECTOR ; 
    Vector e2 = NULL_VECTOR ;
    
    for ( typename BoundaryEdges::const_iterator it = aBdry.begin() ; it != aBdry.end() ; ++ it )
    {
      e1 = e1 + it->v ;
      e2 = e2 + it->n ;
    }     
  
    CGAL_TSMS_LT_TRACE(2,"Adding boundary preservation constrians. SumV=" << xyz_to_string(e1) << " SumN=" << xyz_to_string(e2) );
    
    Matrix H = LT_product(e1);
    
    Vector c = cross_product(e1,e2);
    
    mConstrians.Add_from_gradient(H,c);
  }
}

template<class TSM>
void LindstromTurkCore<TSM>::Add_volume_preservation_constrians( Triangles const& aTriangles )
{
  CGAL_TSMS_LT_TRACE(2,"Adding volume preservation constrians. " << aTriangles.size() << " triangles.");
  
  Vector lSumV = NULL_VECTOR ;
  FT     lSumL(0) ;
  
  for( typename Triangles::const_iterator it = aTriangles.begin(), eit = aTriangles.end() ; it != eit ; ++it )
  {
    lSumV = lSumV + it->NormalV ;
    lSumL = lSumL + it->NormalL ;  
  }   
  
  mConstrians.Add_if_alpha_compatible(lSumV,lSumL);   

}

template<class TSM>
void LindstromTurkCore<TSM>::Add_boundary_and_volume_optimization_constrians( BoundaryEdges const& aBdry, Triangles const& aTriangles )
{
  CGAL_TSMS_LT_TRACE(2,"Adding boundary and volume optimization constrians. ");
  
  Matrix H = NULL_MATRIX ;
  Vector c = NULL_VECTOR ;

  //
  // Volume optimization  
  //
  for( typename Triangles::const_iterator it = aTriangles.begin(), eit = aTriangles.end() ; it != eit ; ++it )
  {
    Triangle const& lTri = *it ;
    
    H += direct_product(lTri.NormalV,lTri.NormalV) ;
    
    c = c - ( lTri.NormalL * lTri.NormalV ) ;
  }   
  
  CGAL_TSMS_LT_TRACE(3,"Hv:" << matrix_to_string(H) << "\n cv:" << xyz_to_string(c) ) ;
  
  
  if ( aBdry.size() > 0 )
  {
    //
    // Boundary optimization
    //
    Matrix Hb = NULL_MATRIX ;
    Vector cb = NULL_VECTOR ;
    
    for ( typename BoundaryEdges::const_iterator it = aBdry.begin() ; it != aBdry.end() ; ++ it )
    {
      Matrix H = LT_product(it->v);
      Vector c = cross_product(it->v,it->n);
      
      Hb += H ;
      cb = cb + c ;
    }     
    
    CGAL_TSMS_LT_TRACE(3,"Hb:" << matrix_to_string(Hb) << "\n cb:" << xyz_to_string(cb) ) ;
    
    //
    // Weighted average
    //
    FT lScaledBoundaryWeight = FT(9) * mParams.BoundaryWeight * squared_distance ( get_point(mP), get_point(mQ) )  ;
    
    H *= mParams.VolumeWeight ;
    c = c * mParams.VolumeWeight ;
    
    H += lScaledBoundaryWeight * Hb ;
    c = c + ( lScaledBoundaryWeight * cb ) ;
    
    CGAL_TSMS_LT_TRACE(3,"VolW=" << mParams.VolumeWeight << " BdryW=" << mParams.BoundaryWeight << " ScaledBdryW=" << lScaledBoundaryWeight ) ;
    
  }
  
  mConstrians.Add_from_gradient(H,c);
}

template<class TSM>
void LindstromTurkCore<TSM>::Add_shape_optimization_constrians( Link const& aLink )
{
  FT s(aLink.size());
  
  Matrix H (s,0,0
           ,0,s,0
           ,0,0,s
           );
           
  Vector c = NULL_VECTOR ;
  
  for( typename Link::const_iterator it = aLink.begin(), eit = aLink.end() ; it != eit ; ++it )
    c = c + (ORIGIN - get_point(*it)) ;  
           
  CGAL_TSMS_LT_TRACE(2,"Adding shape optimization constrians: Shape vector: " << xyz_to_string(c) );
  
  mConstrians.Add_from_gradient(H,c);
}

template<class TSM>
typename LindstromTurkCore<TSM>::FT
LindstromTurkCore<TSM>::Compute_boundary_cost( Vector const& v, BoundaryEdges const& aBdry )
{
  FT rCost(0);
  for ( typename BoundaryEdges::const_iterator it = aBdry.begin() ; it != aBdry.end() ; ++ it )
  {
    Vector u = (it->t - ORIGIN ) - v ;
    Vector c = cross_product(it->v,u);
    rCost += c*c;  
  }     
  return rCost / FT(4) ;
}

template<class TSM>
typename LindstromTurkCore<TSM>::FT
LindstromTurkCore<TSM>::Compute_volume_cost( Vector const& v, Triangles const& aTriangles )
{
  FT rCost(0);
  
  for( typename Triangles::const_iterator it = aTriangles.begin(), eit = aTriangles.end() ; it != eit ; ++it )
  {
    Triangle const& lTri = *it ;
    
    FT lF = lTri.NormalV * v - lTri.NormalL ;
    
    rCost += ( lF * lF ) ;
    
  }
    
  return rCost / FT(36) ;
}

template<class TSM>
typename LindstromTurkCore<TSM>::FT
LindstromTurkCore<TSM>::Compute_shape_cost( Point const& p, Link const& aLink )
{
  FT rCost(0);
  
  for( typename Link::const_iterator it = aLink.begin(), eit = aLink.end() ; it != eit ; ++it )
    rCost += squared_distance(p,get_point(*it)) ;  
    
  return rCost ;
}

template<class TSM>
void LindstromTurkCore<TSM>::Constrians::Add_if_alpha_compatible( Vector const& Ai, FT const& bi )
{
  double slai = to_double(Ai*Ai) ;
  if ( slai > 0.0 )
  {
    double l = CGAL_NTS sqrt(slai) ;
    
    Vector Ain = Ai / l ;
    FT     bin = bi / l ;
    
    bool lAddIt = true ;
    
    if ( n == 1 )
    {
      FT d01 = A.r0() * Ai  ;
      
      double sla0 = to_double(A.r0() * A.r0()) ;
      double sd01 = to_double(d01 * d01) ;
      
      if ( sd01 > ( sla0 * slai * squared_cos_alpha() ) )
        lAddIt = false ;
    }
    else if ( n == 2 )
    {
      Vector N = cross_product(A.r0(),A.r1());
      
      FT dc012 = N * Ai ;
      
      double slc01  = to_double(N * N) ;
      double sdc012 = to_double(dc012 * dc012);       
      
      if ( sdc012 <= slc01 * slai * squared_sin_alpha() )
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
      ++ n ;
      
      CGAL_TSMS_LT_TRACE(1,"Constrains.A:" << matrix_to_string(A) << "\nConstrains.b:" << xyz_to_string(b) ) ;
    }
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

template<class TSM>
void LindstromTurkCore<TSM>::Constrians::Add_from_gradient ( Matrix const& H, Vector const& c )
{
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
        
        Vector A02( A0.x()*A0.x()
                  , A0.y()*A0.y()
                  , A0.z()*A0.z()
                  );
           
        Vector Q0;      
        switch ( index_of_max_component(A02) ) 
        {
          case 0: Q0 = Vector(- A0.z()/A0.x(),0              ,1              ); break;
          case 1: Q0 = Vector(0              ,- A0.z()/A0.y(),1              ); break;
          case 2: Q0 = Vector(1              ,0              ,- A0.x()/A0.z()); break;

           default : Q0 = NULL_VECTOR ; // This should never happen!
        }
        
        CGAL_assertion( Q0 != NULL_VECTOR ) ;
        
        Vector Q1 = cross_product(A0,Q0);
    
        Vector A1 = H * Q0 ;
        Vector A2 = H * Q1 ;
        FT b1 = - ( Q0 * c ) ;
        FT b2 = - ( Q1 * c ) ;
        
        Add_if_alpha_compatible(A1,b1);
        Add_if_alpha_compatible(A2,b2);
        
      }
      break ;
    
    case 2:
      {
    
        Vector Q = cross_product(A.r0(),A.r1());
        
        Vector A2 = H * Q ;
        
        FT b2 = - ( Q * c ) ;
        
        Add_if_alpha_compatible(A2,b2);
        
      }
      break ;
      
  }
}

} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROMTURK_CORE_IMPL_H //
// EOF //
 

