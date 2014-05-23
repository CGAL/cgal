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
                                , VertexIdxMap     const& 
                                , EdgeIdxMap       const&
                                , EdgeIsBorderMap  const&
                                , bool has_border
                                )
  :
   mV0V1(aV0V1)
  ,mSurface(boost::addressof(aSurface))
  
{
  CGAL_PROFILER("Edge_profile constructor calls");

  mLink.reserve(12);
  mTriangles.reserve(16);
  mV1V0 = opposite_edge(v0_v1(),surface_mesh());
  
  mV0 = source(v0_v1(),surface_mesh());
  mV1 = target(v0_v1(),surface_mesh());
  
  CGAL_assertion( mV0 != mV1 );
  
  mP0 = get(vertex_point,surface_mesh(),mV0);
  mP1 = get(vertex_point,surface_mesh(),mV1);
  
  mIsBorderV0V1 = is_border(v0_v1());
  mIsBorderV1V0 = is_border(v1_v0());
  
  if ( left_face_exists() ) 
  {
    CGAL_SURF_SIMPL_TEST_assertion( !mV0V1->is_border() ) ;

    mVLV0 = prev_edge(v0_v1(),surface_mesh());
    mV1VL = next_edge(v0_v1(),surface_mesh());
    mVL   = target(v1_vL(),surface_mesh());
    
    CGAL_SURF_SIMPL_TEST_assertion( mV0 != mVL );
    CGAL_SURF_SIMPL_TEST_assertion( mVL == source(vL_v0(),surface_mesh()) );
  }
  else
  {
    CGAL_SURF_SIMPL_TEST_assertion( mV0V1->is_border() ) ;
  }
  
  if ( right_face_exists() )
  {
    CGAL_SURF_SIMPL_TEST_assertion( !mV1V0->is_border() ) ;

    mV0VR = next_edge(v1_v0(),surface_mesh());
    mVRV1 = prev_edge(v1_v0(),surface_mesh());
    mVR   = target(v0_vR(),surface_mesh());
    
    CGAL_SURF_SIMPL_TEST_assertion( mV0 != mVR );
    CGAL_SURF_SIMPL_TEST_assertion( mVR == source(vR_v1(),surface_mesh()) );
  }
  else
  {
    CGAL_SURF_SIMPL_TEST_assertion( mV1V0->is_border() ) ;
  }
  
  if(has_border){
    Extract_borders();
  }
}


template<class ECM>
void Edge_profile<ECM>::Extract_borders()
{
  edge_descriptor e = mV0V1;
  edge_descriptor oe = opposite_edge(e, surface_mesh());
  bool b;
  if((b = is_border(e)) || is_border(oe)){
    mBorderEdges.push_back(b?e:oe);
  }
  e = next_edge(oe,surface_mesh());
  oe = opposite_edge(e,surface_mesh());
  while(e != mV0V1){
    if((b = is_border(e)) || is_border(oe)){
      mBorderEdges.push_back(b?e:oe);
    }
    e = next_edge(oe,surface_mesh());
    oe = opposite_edge(e,surface_mesh());
  }
  e = opposite_edge(next_edge(e,surface_mesh()),surface_mesh());
  oe = opposite_edge(e,surface_mesh());
    while(e != mV0V1){
    if((b = is_border(e)) || is_border(oe)){
      mBorderEdges.push_back(b?e:oe);
    } 
    e = opposite_edge(next_edge(e,surface_mesh()),surface_mesh());
    oe = opposite_edge(e,surface_mesh());
    }
}



// Extract all triangles (its normals) and vertices (the link) around the collapsing edge p_q
//
template<class ECM>
void Edge_profile<ECM>::Extract_triangles_and_link()
{
  #ifdef CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK
  std::set<vertex_descriptor> vertex_already_inserted;
  #endif
  // look at the two faces or holes adjacent to edge (v0,v1)
  // and at the opposite vertex if it exists
  edge_descriptor endleft = next_edge(v1_v0(), surface_mesh());
  edge_descriptor endright = next_edge(v0_v1(), surface_mesh());

  if( left_face_exists() )
    mTriangles.push_back(Triangle(v0(),v1(),vL()) ) ;
  if( right_face_exists() )
    mTriangles.push_back(Triangle(v1(),v0(),vR()) ) ;

  // counterclockwise around v0
  edge_descriptor e02 = opposite_edge(prev_edge(v0_v1(),surface_mesh()), surface_mesh());
  vertex_descriptor v, v2 =target(e02,surface_mesh());
  while(e02 != endleft) {
    #ifdef CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK
    if (vertex_already_inserted.insert(v2).second)
    #endif
    mLink.push_back(v2);
    bool is_b = is_border(e02);
    e02 = opposite_edge(prev_edge(e02,surface_mesh()), surface_mesh());
    v = target(e02,surface_mesh());
    if(! is_b){
      mTriangles.push_back(Triangle(v,v0(),v2) ) ;
    }
    v2 = v;
  }
  if(v != vR() && (v!= vertex_descriptor())){
    #ifdef CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK
    if (vertex_already_inserted.insert(v).second)
    #endif
    mLink.push_back(v);
  }

  // counterclockwise around v1
  e02 = opposite_edge(prev_edge(v1_v0(),surface_mesh()), surface_mesh());
  v2 = target(e02,surface_mesh());
  while(e02 != endright) {
    #ifdef CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK
    if (vertex_already_inserted.insert(v2).second)
    #endif
    mLink.push_back(v2);
    bool is_b = is_border(e02);
    e02 = opposite_edge(prev_edge(e02,surface_mesh()), surface_mesh());
    v = target(e02,surface_mesh());
    if(! is_b){
      mTriangles.push_back(Triangle(v,v1(),v2) ) ;
    }
    v2 = v;
  }
  if(v != vL() && (v!= vertex_descriptor())){
    #ifdef CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK
    if (vertex_already_inserted.insert(v).second)
    #endif
    mLink.push_back(v);
  }
  
}

} // namespace Surface_mesh_simplification

} //namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_PROFILE_IMPL_H
// EOF //
 
