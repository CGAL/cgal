// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
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
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_PROFILE_IMPL_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_PROFILE_IMPL_H 1

namespace CGAL {

namespace Surface_mesh_simplification
{

template<class ECM>

template<class VertexIdxMap
        ,class EdgeIdxMap
        ,class EdgeIsBorderMap
        >
Edge_profile<ECM>::Edge_profile ( edge_descriptor  const& aV0V1
                                , ECM&                    aSurface
                                , VertexIdxMap     const& vertex_idx
                                , EdgeIdxMap       const& edge_idx
                                , EdgeIsBorderMap  const& is_border
                                )
  :
   mV0V1(aV0V1)
  ,mSurface(boost::addressof(aSurface))
  
{
  mV1V0 = opposite_edge(v0_v1(),surface());
  
  mV0 = source(v0_v1(),surface());
  mV1 = target(v0_v1(),surface());
  
  CGAL_assertion( mV0 != mV1 );
  
  mP0 = get(vertex_point,surface(),mV0);
  mP1 = get(vertex_point,surface(),mV1);
  
  mIsBorderV0V1 = is_border[v0_v1()];
  mIsBorderV1V0 = is_border[v1_v0()];
  
  if ( left_face_exists() ) 
  {
    CGAL_SURF_SIMPL_TEST_assertion( !mV0V1->is_border() ) ;

    mVLV0 = prev_edge(v0_v1(),surface());
    mV1VL = next_edge(v0_v1(),surface());
    mVL   = target(v1_vL(),surface());
    
    CGAL_SURF_SIMPL_TEST_assertion( mV0 != mVL );
    CGAL_SURF_SIMPL_TEST_assertion( mVL == source(vL_v0(),surface()) );
  }
  else
  {
    CGAL_SURF_SIMPL_TEST_assertion( mV0V1->is_border() ) ;
  }
  
  if ( right_face_exists() )
  {
    CGAL_SURF_SIMPL_TEST_assertion( !mV1V0->is_border() ) ;

    mV0VR = next_edge(v1_v0(),surface());
    mVRV1 = prev_edge(v1_v0(),surface());
    mVR   = target(v0_vR(),surface());
    
    CGAL_SURF_SIMPL_TEST_assertion( mV0 != mVR );
    CGAL_SURF_SIMPL_TEST_assertion( mVR == source(vR_v1(),surface()) );
  }
  else
  {
    CGAL_SURF_SIMPL_TEST_assertion( mV1V0->is_border() ) ;
  }
    
  Extract_triangles_and_link(vertex_idx,is_border);
  Extract_borders(edge_idx,is_border);
}

template<class ECM>
template<class EdgeIdxMap, class EdgeIsBorderMap>
void Edge_profile<ECM>::Extract_borders( vertex_descriptor const& v
                                       , IdxSet&                  rCollected 
                                       , EdgeIdxMap        const& edge_idx
                                       , EdgeIsBorderMap   const& is_border
                                       )
{
  in_edge_iterator eb, ee ; 
  for ( boost::tie(eb,ee) = in_edges(v,surface()) ; eb != ee ; ++ eb )
  {
    edge_descriptor edge     = *eb ;
    edge_descriptor opp_edge = opposite_edge(edge,surface()) ;
    
    bool is_edge_border     = is_border[edge] ;
    bool is_opp_edge_border = is_border[opp_edge] ;
    
    if ( is_edge_border || is_opp_edge_border )
    {
      std::size_t eidx = edge_idx[edge];
      bool lNotCollected = rCollected.find(eidx) == rCollected.end() ;
      if ( lNotCollected )
      {  
        rCollected.insert(eidx);
        rCollected.insert(edge_idx[opp_edge]);
        
        edge_descriptor border_edge = is_edge_border ? edge : opp_edge ;
      
        mBorderEdges.push_back(border_edge) ;
      }  
    }  
  }
}

template<class ECM>
template<class EdgeIdxMap, class EdgeIsBorderMap>
void Edge_profile<ECM>::Extract_borders( EdgeIdxMap const& edge_idx, EdgeIsBorderMap  const& is_border)
{
  IdxSet lCollected ;
  Extract_borders(mV0,lCollected,edge_idx,is_border);
  Extract_borders(mV1,lCollected,edge_idx,is_border);
}

//
// If (v0,v1,v2) is a finite triangular facet of the mesh, that is, NONE of these vertices are boundary vertices,
// the triangle, properly oriented, is added to mTriangles.
//
template<class ECM>
template<class EdgeIsBorderMap>
void Edge_profile<ECM>::Extract_triangle( vertex_descriptor const& v0
                                        , vertex_descriptor const& v1
                                        , vertex_descriptor const& v2 
                                        , edge_descriptor   const& e02
                                        , EdgeIsBorderMap   const& is_border
                                        )
{
  // The 3 vertices are obtained by circulating ccw around v0, that is, e02 = next_ccw(e01).
  // Since these vertices are NOT obtained by circulating the face, the actual triangle orientation is unspecified.
  
  // The triangle is oriented v0->v2->v1 if the next edge that follows e02 (which is the edge v0->v2) is v2->v1.
  if ( target(next_edge(e02,surface()),surface()) == v1 ) 
  {
    // The triangle is oriented v0->v2->v1.
    // In this case e02 is an edge of the facet.
    // If this facet edge is a border edge then this triangle is not in the mesh .
    if ( !is_border[e02] )
      mTriangles.push_back(Triangle(v0,v2,v1) ) ;
  }
  else
  {
    // The triangle is oriented v0->v1->v2.
    // In this case, e20 and not e02, is an edge of the facet.
    // If this facet edge is a border edge then this triangle is not in the mesh .
    if ( !is_border[opposite_edge(e02,surface())] )
      mTriangles.push_back(Triangle(v0,v1,v2) ) ;
  }
}

//
// Extract all triangles (its normals) and vertices (the link) around the collpasing edge p_q
//
template<class ECM>
template<class VertexIdxMap, class EdgeIsBorderMap>
void Edge_profile<ECM>::Extract_triangles_and_link( VertexIdxMap const&    vertex_idx
                                                  , EdgeIsBorderMap const& is_border 
                                                  )
{
  
  IdxSet lCollected ;
  
  // 
  // Extract around mV0, CCW
  //  
  vertex_descriptor v0 = mV0;
  vertex_descriptor v1 = mV1;
  
  edge_descriptor e02 = mV0V1;
  
  do
  {
    e02 = opposite_edge(prev_edge(e02,surface()), surface());
    vertex_descriptor v2 = target(e02,surface());
  
    if ( v2 != mV1 )
    {
      mLink.push_back(v2) ;
      CGAL_assertion_code( bool lInserted = ) lCollected.insert(vertex_idx[v2]).second ;
      CGAL_assertion(lInserted);
    }
      
    Extract_triangle(v0,v1,v2,e02,is_border);
    
    v1 = v2 ;
  }
  while ( e02 != mV0V1 ) ;
  
  // 
  // Extract around mV1, CCW
  //  
  
  v0 = mV1;
  
  e02 = opposite_edge(prev_edge(mV1V0,surface()), surface());
  
  v1 = target(e02,surface()); 

    
  // This could have been added to the link while circulating around mP
  if ( v1 != mV0 && lCollected.find(vertex_idx[v1]) == lCollected.end() )
    mLink.push_back(v1) ;
  
  e02 = opposite_edge(prev_edge(e02,surface()), surface());
  
  do
  {
    vertex_descriptor v2 = target(e02,surface());

    // Any of the vertices found around mP can be reached again around mQ, but we can't duplicate them here.
    if ( v2 != mV0 && lCollected.find(vertex_idx[v2]) == lCollected.end() )
      mLink.push_back(v2) ;
    
    Extract_triangle(v0,v1,v2,e02,is_border);
    
    v1 = v2 ;
     
    e02 = opposite_edge(prev_edge(e02,surface()), surface());
  }
  while ( e02 != mV1V0 ) ;
}

} // namespace Surface_mesh_simplification

} //namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_PROFILE_IMPL_H
// EOF //
 
