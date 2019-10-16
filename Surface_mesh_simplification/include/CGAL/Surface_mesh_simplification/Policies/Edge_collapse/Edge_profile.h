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

#include <CGAL/Surface_mesh_simplification/internal/Common.h>

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

    Triangle(const vertex_descriptor v0_,
             const vertex_descriptor v1_,
             const vertex_descriptor v2_)
      : v0(v0_), v1(v1_), v2(v2_)
    {
      CGAL_assertion(handle_assigned(v0) && handle_assigned(v1) && handle_assigned(v2));
      CGAL_assertion(v0 != v1 && v1 != v2 && v2 != v0);
    }

    vertex_descriptor v0;
    vertex_descriptor v1;
    vertex_descriptor v2;
  };

  typedef std::vector<Triangle>                                               Triangle_vector;

public :
  template<class VertexIndexMap, class HalfedgeIndexMap>
  Edge_profile(const halfedge_descriptor h_v0v1,
               TM& tmesh,
               const VertexIndexMap& vim,
               const VertexPointMap& vpm,
               const HalfedgeIndexMap& him,
               bool has_border);

public :
  const halfedge_descriptor v0_v1() const { return m_h_v0v1; }
  const halfedge_descriptor v1_v0() const { return m_h_v1v0; }

  const vertex_descriptor v0() const { return m_v0; }
  const vertex_descriptor v1() const { return m_v1; }

  // These are null if v0v1 is a border (thius there is no face to its left)
  const vertex_descriptor vL() const { return m_vL; }
  const halfedge_descriptor v1_vL() const { return m_v1vL; }
  const halfedge_descriptor vL_v0() const { return m_vLv0; }

  // These are null if v1v0 is a border (thius there is no face to its left)
  const vertex_descriptor vR() const { return m_vR; }
  const halfedge_descriptor v0_vR() const { return m_v0vR; }
  const halfedge_descriptor vR_v1() const { return m_vRv1; }

  const Triangle_vector& triangles() const
  {
    if(m_triangles.empty())
      const_cast<Edge_profile*>(this)->extract_triangles_and_link();

    CGAL_HISTOGRAM_PROFILER("triangles.size()", m_triangles.size());
    return m_triangles;
  }

  // The cycle of vertices around the edge
  const vertex_descriptor_vector& link() const
  {
    CGAL_PROFILER("link calls");
    if(m_link_vertices.empty())
      const_cast<Edge_profile*>(this)->extract_triangles_and_link();

    return m_link_vertices;
  }

  const halfedge_descriptor_vector& border_edges() const { return m_border_edges; }
  const TM& surface() const { return m_tm; }
  const TM& surface_mesh() const { return m_tm; }
  VertexPointMap vertex_point_map() const { return m_vpm; }

public :
  const Point& p0() const { return m_p0; }
  const Point& p1() const { return m_p1; }

  bool is_v0_v1_a_border() const { return m_is_v0v1_border; }
  bool is_v1_v0_a_border() const { return m_is_v1v0_border; }

  bool left_face_exists() const { return !m_is_v0v1_border; }
  bool right_face_exists() const { return !m_is_v1v0_border; }

private:
  void extract_borders();
  void extract_triangles_and_link();

private:
  halfedge_descriptor m_h_v0v1;
  halfedge_descriptor m_h_v1v0;

  bool m_is_v0v1_border;
  bool m_is_v1v0_border;

  vertex_descriptor m_v0;
  vertex_descriptor m_v1;

  Point m_p0;
  Point m_p1;

  vertex_descriptor m_vL;
  vertex_descriptor m_vR;

  halfedge_descriptor m_v1vL;
  halfedge_descriptor m_vLv0;
  halfedge_descriptor m_v0vR;
  halfedge_descriptor m_vRv1;

  vertex_descriptor_vector m_link_vertices;
  halfedge_descriptor_vector m_border_edges;
  Triangle_vector m_triangles;

  const TM& m_tm;
  VertexPointMap m_vpm;
};

template<class TM, class VertexPointMap>
template<class VertexIndexMap, class HalfedgeIndexMap>
Edge_profile<TM, VertexPointMap>::
Edge_profile(const halfedge_descriptor h_v0v1,
             TM& tmesh,
             const VertexIndexMap& /*vim*/,
             const VertexPointMap& vpm,
             const HalfedgeIndexMap& /*him*/,
             bool has_border)
  :
    m_h_v0v1(h_v0v1),
    m_tm(tmesh),
    m_vpm(vpm)
{
  CGAL_PROFILER("Edge_profile constructor calls");

  m_link_vertices.reserve(12);
  m_triangles.reserve(16);
  m_h_v1v0 = opposite(v0_v1(), surface_mesh());

  m_v0 = source(v0_v1(), surface_mesh());
  m_v1 = target(v0_v1(), surface_mesh());

  CGAL_assertion(m_v0 != m_v1);

  m_p0 = get(vertex_point_map(), m_v0);
  m_p1 = get(vertex_point_map(), m_v1);

  m_is_v0v1_border = is_border(v0_v1(), surface_mesh());
  m_is_v1v0_border = is_border(v1_v0(), surface_mesh());

  if(left_face_exists())
  {
    CGAL_assertion(! is_border(m_h_v0v1, surface_mesh()));

    m_vLv0 = prev(v0_v1(), surface_mesh());
    m_v1vL = next(v0_v1(), surface_mesh());
    m_vL = target(v1_vL(), surface_mesh());

    CGAL_assertion(m_v0 != m_vL);
    CGAL_assertion(m_vL == source(vL_v0(), surface_mesh()));
  }
  else
  {
    CGAL_assertion(is_border(m_h_v0v1, surface_mesh()));
  }

  if(right_face_exists())
  {
    CGAL_assertion(! is_border(m_h_v1v0, surface_mesh()));

    m_v0vR = next(v1_v0(), surface_mesh());
    m_vRv1 = prev(v1_v0(), surface_mesh());
    m_vR = target(v0_vR(), surface_mesh());

    CGAL_assertion(m_v0 != m_vR);
    CGAL_assertion(m_vR == source(vR_v1(), surface_mesh()));
  }
  else
  {
    CGAL_assertion(is_border(m_h_v1v0, surface_mesh()));
  }

  if(has_border)
    extract_borders();
}

template<class TM, class VertexPointMap>
void
Edge_profile<TM,VertexPointMap>::
extract_borders()
{
  halfedge_descriptor e = m_h_v0v1;
  halfedge_descriptor oe = opposite(e, surface_mesh());

  bool b;
  if((b = is_border(e, surface_mesh())) || is_border(oe, surface_mesh()))
    m_border_edges.push_back(b?e:oe);

  e = next(oe, surface_mesh());
  oe = opposite(e, surface_mesh());
  while(e != m_h_v0v1)
  {
    if((b = is_border(e, surface_mesh())) || is_border(oe, surface_mesh()))
      m_border_edges.push_back(b?e:oe);

    e = next(oe, surface_mesh());
    oe = opposite(e, surface_mesh());
  }

  e = opposite(next(e, surface_mesh()), surface_mesh());
  oe = opposite(e, surface_mesh());
  while(e != m_h_v0v1)
  {
    if((b = is_border(e, surface_mesh())) || is_border(oe, surface_mesh()))
      m_border_edges.push_back(b?e:oe);

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
    m_triangles.push_back(Triangle(v0(), v1(), vL()));
  if(right_face_exists())
    m_triangles.push_back(Triangle(v1(), v0(), vR()));

  // counterclockwise around v0
  halfedge_descriptor e02 = opposite(prev(v0_v1(), surface_mesh()), surface_mesh());
  vertex_descriptor v = target(e02, surface_mesh()), v2 = v;

  while(e02 != endleft)
  {
#ifdef CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK
    if(vertex_already_inserted.insert(v2).second)
#endif
      m_link_vertices.push_back(v2);

    bool is_b = is_border(e02, surface_mesh());
    e02 = opposite(prev(e02, surface_mesh()), surface_mesh());
    v = target(e02, surface_mesh());
    if(!is_b)
      m_triangles.push_back(Triangle(v, v0(), v2));

    v2 = v;
  }

  e02 = opposite(prev(v1_v0(), surface_mesh()), surface_mesh());
  if(target(e02, surface_mesh()) != v) // add the vertex if it is not added in the following loop
  {
#ifdef CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK
    if(vertex_already_inserted.insert(v).second)
#endif
      m_link_vertices.push_back(v);
  }

  // counterclockwise around v1
  v2 = target(e02, surface_mesh());
  v = v2;
  while(e02 != endright)
  {
#ifdef CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK
    if(vertex_already_inserted.insert(v2).second)
#endif
      m_link_vertices.push_back(v2);

    bool is_b = is_border(e02, surface_mesh());
    e02 = opposite(prev(e02, surface_mesh()), surface_mesh());
    v = target(e02, surface_mesh());

    if(!is_b)
      m_triangles.push_back(Triangle(v,v1(),v2));

    v2 = v;
  }

  if(m_link_vertices.empty() || //handle link of an isolated triangle
     target(opposite(prev(v0_v1(), surface_mesh()), surface_mesh()), surface_mesh())!=v)
  {
#ifdef CGAL_SMS_EDGE_PROFILE_ALWAYS_NEED_UNIQUE_VERTEX_IN_LINK
    if(vertex_already_inserted.insert(v).second)
#endif
      m_link_vertices.push_back(v);
  }

  CGAL_assertion(!m_link_vertices.empty());
}

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_PROFILE_H
