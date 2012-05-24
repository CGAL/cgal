// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_PROFILE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_PROFILE_H 1

#include <vector>
#include <set>

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>

namespace CGAL {

namespace Surface_mesh_simplification
{

template<class ECM_>
class Edge_profile
{
public:

  typedef ECM_ ECM ;
  
  typedef boost::graph_traits<ECM const> ConstGraphTraits ;
  typedef boost::graph_traits<ECM>       GraphTraits ;
  
  typedef typename ConstGraphTraits::vertex_descriptor const_vertex_descriptor ;
  typedef typename ConstGraphTraits::edge_descriptor   const_edge_descriptor ;
  
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor ;
  typedef typename GraphTraits::edge_descriptor   edge_descriptor ;

  typedef typename halfedge_graph_traits<ECM>::Point Point ;
    
public:

  struct Triangle
  {
    Triangle() {}
    
    Triangle( vertex_descriptor const& a_v0
            , vertex_descriptor const& a_v1
            , vertex_descriptor const& a_v2
            ) 
            : v0(a_v0), v1(a_v1), v2(a_v2)
    {
      CGAL_SURF_SIMPL_TEST_assertion( handle_assigned(v0) && handle_assigned(v1) && handle_assigned(v2) ) ;
      CGAL_SURF_SIMPL_TEST_assertion( v0 != v1 && v1 != v2 && v2 != v0 ) ;
    }
    
    vertex_descriptor v0 ;
    vertex_descriptor v1 ;
    vertex_descriptor v2 ;
  } ;
  
  typedef std::vector<vertex_descriptor> vertex_descriptor_vector ;
  typedef std::vector<edge_descriptor>   edge_descriptor_vector ;
  
  typedef std::vector<Triangle> Triangle_vector ;
  
public :
  
  template<class VertexIdxMap
          ,class EdgeIdxMap
          ,class EdgeIsBorderMap
          >
  Edge_profile ( edge_descriptor  const& aV0V1
               , ECM&                    aSurface
               , VertexIdxMap     const& aVertex_index_map
               , EdgeIdxMap       const& aEdge_index_map
               , EdgeIsBorderMap  const& aEdge_is_border_map
               ) ;
     
public :

  edge_descriptor const& v0_v1() const { return mV0V1; }
  edge_descriptor const& v1_v0() const { return mV1V0; }
  
  vertex_descriptor const& v0() const { return mV0; }
  vertex_descriptor const& v1() const { return mV1; }

  // These are null if v0v1 is a border (thius there is no face to its left)  
  vertex_descriptor const&  vL() const { return mVL; } 
  edge_descriptor const& v1_vL() const { return mV1VL; }
  edge_descriptor const& vL_v0() const { return mVLV0; }
  
  // These are null if v1v0 is a border (thius there is no face to its left)  
  vertex_descriptor const&  vR() const { return mVR; } 
  edge_descriptor const& v0_vR() const { return mV0VR; }
  edge_descriptor const& vR_v1() const { return mVRV1; }

  Triangle_vector const& triangles() const { return mTriangles ; }
  
  // The cycle of vertices around the edge  
  vertex_descriptor_vector const& link() const { return mLink ; }
  
  edge_descriptor_vector const& border_edges() const { return mBorderEdges ; }

  ECM& surface() const { return *mSurface ; } 
 
  
public :

  Point const& p0() const { return mP0; }  
  Point const& p1() const { return mP1; }  
  
  bool is_v0_v1_a_border() const { return mIsBorderV0V1 ; }
  bool is_v1_v0_a_border() const { return mIsBorderV1V0 ; }
  
  bool left_face_exists () const { return !mIsBorderV0V1 ; }
  bool right_face_exists() const { return !mIsBorderV1V0 ; }
  
private:

  typedef typename GraphTraits::in_edge_iterator  in_edge_iterator ;
  
  typedef std::set<std::size_t> IdxSet;

private:
   
  template<class EdgeIdxMap, class EdgeIsBorderMap>
  void Extract_borders( vertex_descriptor const& v
                      , IdxSet&                  rCollected 
                      , EdgeIdxMap        const& edge_idx
                      , EdgeIsBorderMap   const& is_border
                      ) ;
   
   template<class EdgeIdxMap, class EdgeIsBorderMap>
   void Extract_borders( EdgeIdxMap const& edge_idx, EdgeIsBorderMap  const& is_border) ;
   
   template<class EdgeIsBorderMap>
   void Extract_triangle( vertex_descriptor const& v0
                        , vertex_descriptor const& v1
                        , vertex_descriptor const& v2 
                        , edge_descriptor   const& e02
                        , EdgeIsBorderMap   const& is_border
                        ) ;
                        
   template<class VertexIdxMap, class EdgeIsBorderMap>
   void Extract_triangles_and_link( VertexIdxMap const& vertex_idx, EdgeIsBorderMap const& is_border ) ;
    
private:
 
  edge_descriptor mV0V1;
  edge_descriptor mV1V0;

  bool mIsBorderV0V1 ;
  bool mIsBorderV1V0 ;
  
  vertex_descriptor mV0;
  vertex_descriptor mV1;
  
  Point mP0 ;
  Point mP1 ;
  
  vertex_descriptor mVL;
  vertex_descriptor mVR;

  edge_descriptor mV1VL;
  edge_descriptor mVLV0;
  edge_descriptor mV0VR;
  edge_descriptor mVRV1;
  
  vertex_descriptor_vector mLink ;
  edge_descriptor_vector   mBorderEdges ;
  Triangle_vector          mTriangles ;
  
  ECM* mSurface ;
  
} ;
  
} // namespace Surface_mesh_simplification

} //namespace CGAL

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile_impl.h>

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_PROFILE_H
// EOF //
 
