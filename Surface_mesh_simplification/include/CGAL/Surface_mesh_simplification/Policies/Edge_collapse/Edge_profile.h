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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_PROFILE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_PROFILE_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>

#include <vector>
#include <set>

namespace CGAL {
namespace Surface_mesh_simplification {

template<class TM_,
         class VertexPointMap_ = typename boost::property_map<TM_, CGAL::vertex_point_t>::type>
class Edge_profile
{
public:
  typedef TM_                                                                 TM;
  typedef VertexPointMap_                                                     VertexPointMap;
  typedef boost::graph_traits<TM>                                             GraphTraits;

  typedef typename GraphTraits::vertex_descriptor                             vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor                           halfedge_descriptor;
  typedef typename GraphTraits::face_descriptor                               face_descriptor;

  typedef std::vector<vertex_descriptor>                                      vertex_descriptor_vector;
  typedef std::vector<halfedge_descriptor>                                    halfedge_descriptor_vector;

  typedef typename boost::property_traits<VertexPointMap>::value_type         Point;
  typedef typename Kernel_traits<Point>::Kernel                               Kernel;
  typedef typename Kernel::FT                                                 FT;

public:
  struct Triangle
  {
    Triangle() {}

    Triangle(const vertex_descriptor a_v0,
             const vertex_descriptor a_v1,
             const vertex_descriptor a_v2)
      : v0(a_v0), v1(a_v1), v2(a_v2)
    {
      CGAL_SURF_SIMPL_TEST_assertion(handle_assigned(v0) && handle_assigned(v1) && handle_assigned(v2));
      CGAL_SURF_SIMPL_TEST_assertion(v0 != v1 && v1 != v2 && v2 != v0);
    }

    vertex_descriptor v0;
    vertex_descriptor v1;
    vertex_descriptor v2;
  };

  typedef std::vector<Triangle> Triangle_vector;

public :
  template<class VertexIdxMap,
           class EdgeIdxMap>
  Edge_profile(const halfedge_descriptor& aV0V1,
               TM& aSurface,
               const VertexIdxMap& aVertex_index_map,
               const VertexPointMap& aVertex_point_map,
               const EdgeIdxMap& aEdge_index_map,
               bool has_border);

public :
  const halfedge_descriptor v0_v1() const { return mV0V1; }
  const halfedge_descriptor v1_v0() const { return mV1V0; }

  const vertex_descriptor v0() const { return mV0; }
  const vertex_descriptor v1() const { return mV1; }

  // These are null if v0v1 is a border (thius there is no face to its left)
  const vertex_descriptor vL() const { return mVL; }
  const halfedge_descriptor v1_vL() const { return mV1VL; }
  const halfedge_descriptor vL_v0() const { return mVLV0; }

  // These are null if v1v0 is a border (thius there is no face to its left)
  const vertex_descriptor vR() const { return mVR; }
  const halfedge_descriptor v0_vR() const { return mV0VR; }
  const halfedge_descriptor vR_v1() const { return mVRV1; }

  const Triangle_vector& triangles() const
  {
    if(mTriangles.empty())
      const_cast<Edge_profile*>(this)->extract_triangles_and_link();

    CGAL_HISTOGRAM_PROFILER("triangles.size()", mTriangles.size());
    return mTriangles;
  }

  // The cycle of vertices around the edge
  const vertex_descriptor_vector& link() const
  {
    CGAL_PROFILER("link calls");
    if(mLink.empty())
      const_cast<Edge_profile*>(this)->extract_triangles_and_link();

    return mLink;
  }

  const halfedge_descriptor_vector& border_edges() const { return mBorderEdges; }
  TM& surface() const { return *mSurface; }
  TM& surface_mesh() const { return *mSurface; }
  VertexPointMap vertex_point_map() const { return mvpm; }

public :
  const Point& p0() const { return mP0; }
  const Point& p1() const { return mP1; }

  bool is_v0_v1_a_border() const { return mIsBorderV0V1; }
  bool is_v1_v0_a_border() const { return mIsBorderV1V0; }

  bool left_face_exists () const { return !mIsBorderV0V1; }
  bool right_face_exists() const { return !mIsBorderV1V0; }

private:
  bool is_border(halfedge_descriptor e) const {
    return face(e,*mSurface) == boost::graph_traits<TM>::null_face();
  }

  void extract_borders();
  void extract_triangles_and_link();

private:
  halfedge_descriptor mV0V1;
  halfedge_descriptor mV1V0;

  bool mIsBorderV0V1;
  bool mIsBorderV1V0;

  vertex_descriptor mV0;
  vertex_descriptor mV1;

  Point mP0;
  Point mP1;

  vertex_descriptor mVL;
  vertex_descriptor mVR;

  halfedge_descriptor mV1VL;
  halfedge_descriptor mVLV0;
  halfedge_descriptor mV0VR;
  halfedge_descriptor mVRV1;

  vertex_descriptor_vector mLink;
  halfedge_descriptor_vector mBorderEdges;
  Triangle_vector mTriangles;

  TM* mSurface;
  VertexPointMap mvpm;
};

template<class TM, class VertexPointMap>
template<class VertexIdxMap, class EdgeIdxMap>
Edge_profile<TM, VertexPointMap>::
Edge_profile(const halfedge_descriptor& aV0V1,
             TM& aSurface,
             const VertexIdxMap&,
             const VertexPointMap& aVertex_point_map,
             const EdgeIdxMap&,
             bool has_border)
  :
    mV0V1(aV0V1),
    mSurface(boost::addressof(aSurface)),
    mvpm(aVertex_point_map)
{
  CGAL_PROFILER("Edge_profile constructor calls");

  mLink.reserve(12);
  mTriangles.reserve(16);
  mV1V0 = opposite(v0_v1(), surface_mesh());

  mV0 = source(v0_v1(), surface_mesh());
  mV1 = target(v0_v1(), surface_mesh());

  CGAL_assertion(mV0 != mV1);

  mP0 = get(vertex_point_map(),mV0);
  mP1 = get(vertex_point_map(),mV1);

  mIsBorderV0V1 = is_border(v0_v1());
  mIsBorderV1V0 = is_border(v1_v0());

  if(left_face_exists())
  {
    CGAL_SURF_SIMPL_TEST_assertion(! is_border(mV0V1));

    mVLV0 = prev(v0_v1(), surface_mesh());
    mV1VL = next(v0_v1(), surface_mesh());
    mVL = target(v1_vL(), surface_mesh());

    CGAL_SURF_SIMPL_TEST_assertion(mV0 != mVL);
    CGAL_SURF_SIMPL_TEST_assertion(mVL == source(vL_v0(), surface_mesh()));
  }
  else
  {
    CGAL_SURF_SIMPL_TEST_assertion(is_border(mV0V1));
  }

  if(right_face_exists())
  {
    CGAL_SURF_SIMPL_TEST_assertion(! is_border(mV1V0));

    mV0VR = next(v1_v0(), surface_mesh());
    mVRV1 = prev(v1_v0(), surface_mesh());
    mVR = target(v0_vR(), surface_mesh());

    CGAL_SURF_SIMPL_TEST_assertion(mV0 != mVR);
    CGAL_SURF_SIMPL_TEST_assertion(mVR == source(vR_v1(), surface_mesh()));
  }
  else
  {
    CGAL_SURF_SIMPL_TEST_assertion(is_border(mV1V0));
  }

  if(has_border)
    extract_borders();
}

template<class TM, class VertexPointMap>
void
Edge_profile<TM,VertexPointMap>::
extract_borders()
{
  halfedge_descriptor e = mV0V1;
  halfedge_descriptor oe = opposite(e, surface_mesh());

  bool b;
  if((b = is_border(e)) || is_border(oe))
    mBorderEdges.push_back(b?e:oe);

  e = next(oe, surface_mesh());
  oe = opposite(e, surface_mesh());
  while(e != mV0V1)
  {
    if((b = is_border(e)) || is_border(oe))
      mBorderEdges.push_back(b?e:oe);

    e = next(oe, surface_mesh());
    oe = opposite(e, surface_mesh());
  }

  e = opposite(next(e, surface_mesh()), surface_mesh());
  oe = opposite(e, surface_mesh());
  while(e != mV0V1)
  {
    if((b = is_border(e)) || is_border(oe))
      mBorderEdges.push_back(b?e:oe);

    e = opposite(next(e, surface_mesh()), surface_mesh());
    oe = opposite(e, surface_mesh());
  }
}

// Extract all triangles (its normals) and vertices (the link) around the collapsing edge p_q
template<class TM, class VertexPointMap>
void
Edge_profile<TM,VertexPointMap>::
extract_triangles_and_link()
{
#ifdef CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK
  std::set<vertex_descriptor> vertex_already_inserted;
#endif
  // look at the two faces or holes adjacent to edge (v0,v1)
  // and at the opposite vertex if it exists
  halfedge_descriptor endleft = next(v1_v0(), surface_mesh());
  halfedge_descriptor endright = next(v0_v1(), surface_mesh());

  if(left_face_exists())
    mTriangles.push_back(Triangle(v0(), v1(), vL()));
  if(right_face_exists())
    mTriangles.push_back(Triangle(v1(), v0(), vR()));

  // counterclockwise around v0
  halfedge_descriptor e02 = opposite(prev(v0_v1(), surface_mesh()), surface_mesh());
  vertex_descriptor v = target(e02, surface_mesh()), v2 = v;

  while(e02 != endleft)
  {
#ifdef CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK
    if(vertex_already_inserted.insert(v2).second)
#endif
      mLink.push_back(v2);

    bool is_b = is_border(e02);
    e02 = opposite(prev(e02, surface_mesh()), surface_mesh());
    v = target(e02, surface_mesh());
    if(!is_b)
      mTriangles.push_back(Triangle(v, v0(), v2));

    v2 = v;
  }

  e02 = opposite(prev(v1_v0(), surface_mesh()), surface_mesh());
  if(target(e02, surface_mesh()) != v) // add the vertex if it is not added in the following loop
  {
#ifdef CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK
    if(vertex_already_inserted.insert(v).second)
#endif
      mLink.push_back(v);
  }

  // counterclockwise around v1
  v2 = target(e02, surface_mesh());
  v = v2;
  while(e02 != endright)
  {
#ifdef CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK
    if(vertex_already_inserted.insert(v2).second)
#endif
      mLink.push_back(v2);

    bool is_b = is_border(e02);
    e02 = opposite(prev(e02, surface_mesh()), surface_mesh());
    v = target(e02, surface_mesh());

    if(!is_b)
      mTriangles.push_back(Triangle(v,v1(),v2));

    v2 = v;
  }

  if(mLink.empty() || //handle link of an isolated triangle
     target(opposite(prev(v0_v1(), surface_mesh()), surface_mesh()), surface_mesh())!=v)
  {
#ifdef CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK
    if(vertex_already_inserted.insert(v).second)
#endif
      mLink.push_back(v);
  }

  CGAL_assertion(!mLink.empty());
}

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_PROFILE_H
