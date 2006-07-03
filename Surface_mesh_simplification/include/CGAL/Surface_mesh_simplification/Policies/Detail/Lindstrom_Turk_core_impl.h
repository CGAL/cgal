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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_LINDSTROM_TURK_CORE_IMPL_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_LINDSTROM_TURK_CORE_IMPL_H 1

CGAL_BEGIN_NAMESPACE

struct G
{
  ~G()
  {
    std::cout << "PFixed=" << PFixed << " QFixed=" << QFixed << " Boundary=" << Boundary << std::endl ;
  }
  
  int PFixed ;
  int QFixed ;
  int Boundary ;
} ;
G g ;

//
// Implementation of the strategy from:
//
//  "Fast and Memory Efficient Polygonal Symplification"
//  Peter Lindstrom, Greg Turk
//

namespace Triangulated_surface_mesh { namespace Simplification 
{

template<class CD>
LindstromTurkCore<CD>::LindstromTurkCore( Params const&            aParams
                                        , vertex_descriptor const& aP
                                        , vertex_descriptor const& aQ
                                        , bool                     aIsPFixed
                                        , bool                     aIsQFixed
                                        , edge_descriptor const&   aP_Q
                                        , edge_descriptor const&   aQ_P
                                        , TSM&                     aSurface 
                                        )
  :
   mParams(aParams)
  ,mP(aP)
  ,mQ(aQ)
  ,mIsPFixed(aIsPFixed)
  ,mIsQFixed(aIsQFixed)
  ,mP_Q(aP_Q)
  ,mQ_P(aQ_P)
  ,mSurface(aSurface)    
{
}

template<class CD>
void LindstromTurkCore<CD>::compute( Collapse_data& rData  )
{

  Optional_FT    lOptionalCost ;
  Optional_point lOptionalP ;
  
  if ( !mIsPFixed || !mIsQFixed )
  {
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
    
    // If the collapsing edge is a boundary edge, the "local boundary" is cached in a Boundary object.    
    OptionalBoundary lBdry ;
    
    if ( is_undirected_edge_a_border(mP_Q) )
      lBdry = Extract_boundary();
      
    Point const& lP = get_point(mP) ;
    Point const& lQ = get_point(mQ) ;
    
    Optional_vector lOptionalV ;
    
    if ( mIsPFixed )
    {
      lOptionalV = Optional_vector(lP - ORIGIN);
    }  
    else if ( mIsQFixed )
    {
      lOptionalV = Optional_vector(lQ - ORIGIN);
    }  
    else
    {
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
      if ( lBdry)
        Add_boundary_preservation_constrians(*lBdry);
      
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
          
        }
        else
          CGAL_TSMS_LT_TRACE(1,"Can't solve optimization, singular system.");
      }
      else
        CGAL_TSMS_LT_TRACE(1,"Can't solve optimization, not enough alpha-compatible constrians.");
    }  
    
    if ( lOptionalV )
    {
      //
      // lP is the optimized new vertex position. Now we compute the collapsing cost:
      //
      lOptionalP = Optional_point( ORIGIN + (*lOptionalV) );
      
      FT lSquaredLength = squared_distance(lP,lQ);
      
      CGAL_TSMS_LT_TRACE(1,"Squared edge length: " << lSquaredLength ) ;
      
      FT lBdryCost = 0;
      if ( lBdry )
        lBdryCost = Compute_boundary_cost(*lOptionalV,*lBdry);
      
      FT lVolumeCost = Compute_volume_cost(*lOptionalV,lTriangles);
      
      FT lShapeCost = Compute_shape_cost(*lOptionalP,lLink);
      
      FT lTotalCost =   mParams.VolumeWeight   * lVolumeCost
                      + mParams.BoundaryWeight * lBdryCost   * lSquaredLength
                      + mParams.ShapeWeight    * lShapeCost  * lSquaredLength * lSquaredLength ;
      
      lOptionalCost = Optional_FT(lTotalCost);
      
      CGAL_TSMS_LT_TRACE(1,"New vertex point: " << xyz_to_string(*lOptionalP) << "\nTotal cost: " << lTotalCost );
                        
      CGAL_TSMS_LT_TRACE(2,"\nSquared edge length: " << lSquaredLength 
                        << "\nBoundary cost: " << lBdryCost 
                        << "\nVolume cost: " << lVolumeCost 
                        << "\nShape cost: " << lShapeCost 
                        );
    }
    
  }
  else
    CGAL_TSMS_LT_TRACE(1,"The edge is a fixed edge.");
  
    
  rData = Collapse_data(mP,mQ,mIsPFixed,mIsQFixed,mP_Q,mSurface,lOptionalCost,lOptionalP) ;
}

//
// Caches the "local boundary", that is, the sequence of 3 border edges: o->p, p->q, q->e 
//
template<class CD>
typename LindstromTurkCore<CD>::OptionalBoundary LindstromTurkCore<CD>::Extract_boundary()
{
  // Since p_q is a boundary edge, one of the previous edges (ccw or cw) is the previous boundary edge
  // Likewise, one of the next edges (ccw or cw) is the next boundary edge.
  edge_descriptor p_pt = next_edge_ccw(mP_Q,mSurface);
  edge_descriptor p_pb = next_edge_cw (mP_Q,mSurface);
  edge_descriptor q_qt = next_edge_cw (mQ_P,mSurface);
  edge_descriptor q_qb = next_edge_ccw(mQ_P,mSurface);
  
  edge_descriptor border_1 = mP_Q;
  edge_descriptor border_0 = is_undirected_edge_a_border(p_pt) ? p_pt : p_pb ;
  edge_descriptor border_2 = is_undirected_edge_a_border(q_qt) ? q_qt : q_qb ;
  
  CGAL_assertion(is_undirected_edge_a_border(border_0));
  CGAL_assertion(is_undirected_edge_a_border(border_2));

  // opposite(border0)->border1->border2 is the local boundary
  
  vertex_descriptor ov = target(border_0,mSurface);
  vertex_descriptor rv = target(border_2,mSurface);
  
  // o->p->q->r is the local boundary
  
  Point const& o = get_point(ov);
  Point const& p = get_point(mP);
  Point const& q = get_point(mQ);
  Point const& r = get_point(rv);
  
  //
  // The "Boundary" object caches the boundary as displacement vectors since the code uses that.
  //
  
  Vector op  = p - o ;
  Vector opN = Point_cross_product(p,o);
  
  Vector pq  = q - p ;
  Vector pqN = Point_cross_product(q,p);
  
  Vector qr  = r - q ;
  Vector qrN = Point_cross_product(r,q);
  
  return OptionalBoundary(Boundary(op,opN,pq,pqN,qr,qrN)) ;
}

//
// Calculates the normal of the triangle (v0,v1,v2) (both vector and its length as (v0xv1).v2)
//
template<class CD>
typename LindstromTurkCore<CD>::Triangle LindstromTurkCore<CD>::Get_triangle( vertex_descriptor const& v0
                                                                            , vertex_descriptor const& v1
                                                                            , vertex_descriptor const& v2 
                                                                            )
{
  CGAL_TSMS_LT_TRACE(3,"Extracting triangle v" << v0->ID << "->v" << v1->ID << "->v" << v2->ID );
  
  Point const& p0 = get_point(v0);
  Point const& p1 = get_point(v1);
  Point const& p2 = get_point(v2);
  
  Vector v01 = p1 - p0 ;
  Vector v02 = p2 - p0 ;
  
  Vector lNormalV = cross_product(v01,v02);
  
  FT lNormalL = Point_cross_product(p0,p1) * (p2-ORIGIN);
  
  return Triangle(lNormalV,lNormalL);
}                              

//
// If (v0,v1,v2) is a finite triangular facet of the mesh, that is, NONE of these vertices are boundary vertices,
// the triangle (properly oriented) is added to rTriangles.
// The triangle is encoded as its normal, calculated using the actual facet orientation [(v0,v1,v2) or (v0,v2,v1)]
//
template<class CD>
void LindstromTurkCore<CD>::Extract_triangle( vertex_descriptor const& v0
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
template<class CD>
void LindstromTurkCore<CD>::Extract_triangles_and_link( Triangles& rTriangles, Link& rLink )
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
  
  v1 = target(e02,mSurface); // This was added to the link while circulating around mP
  
  e02 = next_edge_ccw(e02,mSurface);
  
  do
  {
    vertex_descriptor v2 = target(e02,mSurface);

    // Any of the vertices found around mP can be reached again around mQ, but we can't duplicate them here.
    if ( std::find(rLink.begin(),rLink.end(),v2) == rLink.end() )
      rLink.push_back(v2) ;
    
    Extract_triangle(v0,v1,v2,e02,rTriangles);
    
    v1 = v2 ;
     
    e02 = next_edge_ccw(e02,mSurface);
    
  }
  while ( e02 != mQ_P ) ;
}

template<class CD>
void LindstromTurkCore<CD>::Add_boundary_preservation_constrians( Boundary const& aBdry )
{
  CGAL_TSMS_LT_TRACE(2,"Adding boundary preservation constrians. ");
  
  Vector e1 = aBdry.op  + aBdry.pq  + aBdry.qr ;
  Vector e3 = aBdry.opN + aBdry.pqN + aBdry.qrN ;

  Matrix H = LT_product(e1);
  
  Vector c = cross_product(e1,e3);
  
  mConstrians.Add_from_gradient(H,c);
}

template<class CD>
void LindstromTurkCore<CD>::Add_volume_preservation_constrians( Triangles const& aTriangles )
{
  CGAL_TSMS_LT_TRACE(2,"Adding volume preservation constrians. " << aTriangles.size() << " triangles.");
  
  Vector lSumV = NULL_VECTOR ;
  FT     lSumL(0) ;
  
  for( typename Triangles::const_iterator it = aTriangles.begin(), eit = aTriangles.end() ; it != eit ; ++it )
  {
    CGAL_TSMS_LT_TRACE(2,"V:" << xyz_to_string(it->NormalV) << ", L:" << it->NormalL);
    
    lSumV = lSumV + it->NormalV ;
    lSumL = lSumL + it->NormalL ;  
  }   
  
  
  mConstrians.Add_if_alpha_compatible(lSumV,lSumL);   

}

template<class CD>
void LindstromTurkCore<CD>::Add_boundary_and_volume_optimization_constrians( OptionalBoundary const& aBdry, Triangles const& aTriangles )
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
  
  
  if ( aBdry )
  {
    //
    // Boundary optimization
    //
    Matrix Hb = LT_product(aBdry->op) + LT_product(aBdry->pq) + LT_product(aBdry->qr) ;
    
    Vector cb =  cross_product(aBdry->op,aBdry->opN) + cross_product(aBdry->pq,aBdry->pqN) + cross_product(aBdry->qr,aBdry->qrN);
    
    //
    // Weighted average
    //
    FT lBoundaryWeight = ( FT(9) * mParams.BoundaryWeight * squared_distance ( get_point(mP), get_point(mQ) ) ) / FT(10) ;
    
    H *= mParams.VolumeWeight ;
    c = c * mParams.VolumeWeight ;
    
    H += lBoundaryWeight * Hb ;
    c = c + ( lBoundaryWeight * cb ) ;
    
  }
  
  mConstrians.Add_from_gradient(H,c);
}

template<class CD>
void LindstromTurkCore<CD>::Add_shape_optimization_constrians( Link const& aLink )
{
  CGAL_TSMS_LT_TRACE(2,"Add shape optimization constrians. ");
  
  FT s(aLink.size());
  
  Matrix H (s,0,0
           ,0,s,0
           ,0,0,s
           );
           
  Vector c = NULL_VECTOR ;
  
  for( typename Link::const_iterator it = aLink.begin(), eit = aLink.end() ; it != eit ; ++it )
    c = c + (ORIGIN - get_point(*it)) ;  
           
  mConstrians.Add_from_gradient(H,c);
}

template<class CD>
typename LindstromTurkCore<CD>::FT
LindstromTurkCore<CD>::Compute_boundary_cost( Vector const& v, Boundary const& aBdry )
{
  FT rCost(0);
  return rCost ;
}

template<class CD>
typename LindstromTurkCore<CD>::FT
LindstromTurkCore<CD>::Compute_volume_cost( Vector const& v, Triangles const& aTriangles )
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

template<class CD>
typename LindstromTurkCore<CD>::FT
LindstromTurkCore<CD>::Compute_shape_cost( Point const& p, Link const& aLink )
{
  FT rCost(0);
  
  for( typename Link::const_iterator it = aLink.begin(), eit = aLink.end() ; it != eit ; ++it )
    rCost += squared_distance(p,get_point(*it)) ;  
    
  return rCost ;
}

/*
template<class CD>
void LindstromTurkCore<CD>::Constrians::Add_if_alpha_compatible( Vector const& Ai, FT const& bi )
{
  double slai = to_double(Ai*Ai) ;
  if ( slai > 0.0 )
  {
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
          A.r0() = Ai ;
          b = Vector(bi,b.y(),b.z());
          break ;
        case 1 :
          A.r1() = Ai ;
          b = Vector(b.x(),bi,b.z());
          break ;
        case 2 :
          A.r2() = Ai ;
          b = Vector(b.x(),b.y(),bi);
          break ;
      }
      ++ n ;
    }
  }
}
*/

template<class CD>
void LindstromTurkCore<CD>::Constrians::Add_if_alpha_compatible( Vector const& Ai, FT const& bi )
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
    }
  }
  
  CGAL_TSMS_LT_TRACE(1,"Constrains.A:" << matrix_to_string(A) << "\nConstrains.b:" << xyz_to_string(b) ) ;
  
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

template<class CD>
void LindstromTurkCore<CD>::Constrians::Add_from_gradient ( Matrix const& H, Vector const& c )
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

           default : Q0 = NULL_VECTOR ; // This should never happen ;
        }
        
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

} } // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_LINDSTROMTURK_CORE_IMPL_H //
// EOF //
 

