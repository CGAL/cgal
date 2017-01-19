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

#include <CGAL/license/Surface_mesh_simplification.h>


#include <vector>
#include <set>

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>

namespace CGAL {

namespace Surface_mesh_simplification
{

  template<class ECM_, class VertexPointMap_ = typename boost::property_map<ECM_, CGAL::vertex_point_t>::type>
class Edge_profile
{
public:

  typedef ECM_ ECM ;
  typedef VertexPointMap_ VertexPointMap;
  typedef boost::graph_traits<ECM>       GraphTraits ;
  
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor ;
  typedef typename GraphTraits::face_descriptor face_descriptor ;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor ;


  //typedef typename boost::property_map<ECM, CGAL::vertex_point_t>::type Vertex_point_pmap;
  typedef typename boost::property_traits<VertexPointMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  typedef typename Kernel::FT FT;

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
  typedef std::vector<halfedge_descriptor>   halfedge_descriptor_vector ;
  
  typedef std::vector<Triangle> Triangle_vector ;
  
public :
  
  template<class VertexIdxMap
          ,class EdgeIdxMap
          >
  Edge_profile ( halfedge_descriptor  const& aV0V1
               , ECM&                    aSurface
               , VertexIdxMap     const& aVertex_index_map
               , VertexPointMap   const& aVertex_point_map
               , EdgeIdxMap       const& aEdge_index_map
               , bool has_border 
               ) ;
     
public :

  halfedge_descriptor const& v0_v1() const { return mV0V1; }
  halfedge_descriptor const& v1_v0() const { return mV1V0; }
  
  vertex_descriptor const& v0() const { return mV0; }
  vertex_descriptor const& v1() const { return mV1; }

  // These are null if v0v1 is a border (thius there is no face to its left)  
  vertex_descriptor const&  vL() const { return mVL; } 
  halfedge_descriptor const& v1_vL() const { return mV1VL; }
  halfedge_descriptor const& vL_v0() const { return mVLV0; }
  
  // These are null if v1v0 is a border (thius there is no face to its left)  
  vertex_descriptor const&  vR() const { return mVR; } 
  halfedge_descriptor const& v0_vR() const { return mV0VR; }
  halfedge_descriptor const& vR_v1() const { return mVRV1; }

  Triangle_vector const& triangles() const {

    if(mTriangles.empty()){
      const_cast<Edge_profile*>(this)->Extract_triangles_and_link();
    }
    CGAL_HISTOGRAM_PROFILER("triangles.size()", mTriangles.size());
    return mTriangles ; }
  
  // The cycle of vertices around the edge  
  vertex_descriptor_vector const& link() const {
    CGAL_PROFILER("link calls");
    if(mLink.empty()){

      const_cast<Edge_profile*>(this)->Extract_triangles_and_link();
    }
    return mLink ; }
  
  halfedge_descriptor_vector const& border_edges() const {
    return mBorderEdges ; 
  }
  ECM& surface() const { return *mSurface ; } 
  ECM& surface_mesh() const { return *mSurface ; } 
 
  VertexPointMap vertex_point_map() const { return mvpm ; }
  
public :

  Point const& p0() const { return mP0; }  
  Point const& p1() const { return mP1; }  
  
  bool is_v0_v1_a_border() const { return mIsBorderV0V1 ; }
  bool is_v1_v0_a_border() const { return mIsBorderV1V0 ; }
  
  bool left_face_exists () const { return !mIsBorderV0V1 ; }
  bool right_face_exists() const { return !mIsBorderV1V0 ; }
  
private:

  //  typedef typename GraphTraits::in_edge_iterator  in_edge_iterator ;
  
  bool is_border(halfedge_descriptor e) const
  {
    return face(e,*mSurface) == boost::graph_traits<ECM>::null_face();
  }
   

  void Extract_borders() ;
   
  void Extract_triangles_and_link() ;
    
private:

  halfedge_descriptor mV0V1;
  halfedge_descriptor mV1V0;

  bool mIsBorderV0V1 ;
  bool mIsBorderV1V0 ;
  
  vertex_descriptor mV0;
  vertex_descriptor mV1;
  
  Point mP0 ;
  Point mP1 ;
  
  vertex_descriptor mVL;
  vertex_descriptor mVR;

  halfedge_descriptor mV1VL;
  halfedge_descriptor mVLV0;
  halfedge_descriptor mV0VR;
  halfedge_descriptor mVRV1;
  
  vertex_descriptor_vector mLink ;
  halfedge_descriptor_vector   mBorderEdges ;
  Triangle_vector          mTriangles ;
  
  ECM* mSurface ;
  VertexPointMap mvpm;  
} ;
  
} // namespace Surface_mesh_simplification

} //namespace CGAL

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile_impl.h>

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_PROFILE_H
// EOF //
 
