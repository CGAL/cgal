// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_3_SUPPORT_PLANE_H
#define CGAL_KSR_3_SUPPORT_PLANE_H

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// CGAL includes.
#include <CGAL/Surface_mesh.h>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR_3/Intersection_graph.h>

namespace CGAL {
namespace KSR_3 {

template<typename GeomTraits>
class Support_plane {

public:
  using Kernel = GeomTraits;

  using FT        = typename Kernel::FT;
  using Point_2   = typename Kernel::Point_2;
  using Point_3   = typename Kernel::Point_3;
  using Vector_2  = typename Kernel::Vector_2;
  using Vector_3  = typename Kernel::Vector_3;
  using Segment_2 = typename Kernel::Segment_2;
  using Segment_3 = typename Kernel::Segment_3;
  using Line_2    = typename Kernel::Line_2;
  using Line_3    = typename Kernel::Line_3;
  using Plane_3   = typename Kernel::Plane_3;

  using Mesh = CGAL::Surface_mesh<Point_2>;
  using Intersection_graph = KSR_3::Intersection_graph<Kernel>;

  using IVertex = typename Intersection_graph::Vertex_descriptor;
  using IEdge   = typename Intersection_graph::Edge_descriptor;

  using Vertex_index   = typename Mesh::Vertex_index;
  using Face_index     = typename Mesh::Face_index;
  using Edge_index     = typename Mesh::Edge_index;
  using Halfedge_index = typename Mesh::Halfedge_index;

  using V_vector_map   = typename Mesh::template Property_map<Vertex_index, Vector_2>;
  using V_ivertex_map  = typename Mesh::template Property_map<Vertex_index, IVertex>;
  using V_iedge_map    = typename Mesh::template Property_map<Vertex_index, IEdge>;
  using V_bool_map     = typename Mesh::template Property_map<Vertex_index, bool>;
  using E_iedge_map    = typename Mesh::template Property_map<Edge_index, IEdge>;
  using F_index_map    = typename Mesh::template Property_map<Face_index, KSR::size_t>;
  using F_uint_map     = typename Mesh::template Property_map<Face_index, unsigned int>;
  using V_original_map = typename Mesh::template Property_map<Vertex_index, bool>;
  using V_time_map     = typename Mesh::template Property_map<Vertex_index, FT>;

private:
  struct Data {
    Plane_3 plane;
    Mesh mesh;
    V_vector_map direction;
    V_ivertex_map v_ivertex_map;
    V_iedge_map v_iedge_map;
    V_bool_map v_active_map;
    E_iedge_map e_iedge_map;
    F_index_map input_map;
    F_uint_map k_map;
    V_original_map v_original_map;
    V_time_map v_time_map;
    std::set<IEdge> iedges;
  };

  std::shared_ptr<Data> m_data;

public:
  Support_plane() { }

  template<typename PointRange>
  Support_plane(const PointRange& points) :
  m_data(std::make_shared<Data>()) {

    // Newell's method.
    Vector_3 normal = CGAL::NULL_VECTOR;
    const std::size_t n = points.size();
    for (std::size_t i = 0; i < n; ++i) {
      const auto& pa = points[i];
      const auto& pb = points[(i + 1) % n];
      const FT x = normal.x() + (pa.y() - pb.y()) * (pa.z() + pb.z());
      const FT y = normal.y() + (pa.z() - pb.z()) * (pa.x() + pb.x());
      const FT z = normal.z() + (pa.x() - pb.x()) * (pa.y() + pb.y());
      normal = Vector_3(x, y, z);
    }
    CGAL_assertion_msg(normal != CGAL::NULL_VECTOR, "ERROR: polygon is flat!");

    m_data->plane = Plane_3(points[0], KSR::normalize(normal));
    m_data->direction = m_data->mesh.template add_property_map<Vertex_index, Vector_2>(
      "v:direction", CGAL::NULL_VECTOR).first;
    m_data->v_ivertex_map = m_data->mesh.template add_property_map<Vertex_index, IVertex>(
      "v:ivertex", Intersection_graph::null_ivertex()).first;
    m_data->v_iedge_map = m_data->mesh.template add_property_map<Vertex_index, IEdge>(
      "v:iedge", Intersection_graph::null_iedge()).first;
    m_data->v_active_map = m_data->mesh.template add_property_map<Vertex_index, bool>(
      "v:active", true).first;
    m_data->e_iedge_map = m_data->mesh.template add_property_map<Edge_index, IEdge>(
      "e:iedge", Intersection_graph::null_iedge()).first;
    m_data->input_map = m_data->mesh.template add_property_map<Face_index, KSR::size_t>(
      "f:input", KSR::no_element()).first;
    m_data->k_map = m_data->mesh.template add_property_map<Face_index, unsigned int>(
      "f:k", 0).first;
    m_data->v_original_map = m_data->mesh.template add_property_map<Vertex_index, bool>(
      "v:original", false).first;
    m_data->v_time_map = m_data->mesh.template add_property_map<Vertex_index, FT>(
      "v:time", FT(0)).first;
  }

  const std::array<Vertex_index, 4>
  add_bbox_polygon(
    const std::array<Point_2, 4>& points,
    const std::array<IVertex, 4>& ivertices) {

    std::array<Vertex_index, 4> vertices;
    for (std::size_t i = 0; i < 4; ++i) {
      const auto vi = m_data->mesh.add_vertex(points[i]);
      m_data->v_ivertex_map[vi] = ivertices[i];
      vertices[i] = vi;
    }

    const auto fi = m_data->mesh.add_face(vertices);
    CGAL_assertion(fi != Mesh::null_face());
    m_data->input_map[fi] = KSR::no_element();
    return vertices;
  }

  const KSR::size_t add_input_polygon(
    const std::vector<Point_2>& points,
    const Point_2& centroid,
    const KSR::size_t input_idx) {

    std::vector<Vertex_index> vertices;
    const std::size_t n = points.size();
    CGAL_assertion(n >= 3);
    vertices.reserve(n);

    FT sum_length = FT(0);
    std::vector<Vector_2> directions;
    directions.reserve(n);

    for (const auto& point : points) {
      directions.push_back(Vector_2(centroid, point));
      const FT length = static_cast<FT>(
        CGAL::sqrt(CGAL::to_double(CGAL::abs(directions.back() * directions.back()))));
      sum_length += length;
    }
    CGAL_assertion(directions.size() == n);
    sum_length /= static_cast<FT>(n);

    for (std::size_t i = 0; i < n; ++i) {
      const auto& point = points[i];
      const auto vi = m_data->mesh.add_vertex(point);
      m_data->direction[vi] = directions[i] / sum_length;
      m_data->v_original_map[vi] = true;
      vertices.push_back(vi);
    }

    const auto fi = m_data->mesh.add_face(vertices);
    CGAL_assertion(fi != Mesh::null_face());
    m_data->input_map[fi] = input_idx;
    return static_cast<KSR::size_t>(fi);
  }

  // OTHER
  const Plane_3& plane() const { return m_data->plane; }

  const Mesh& mesh() const { return m_data->mesh; }
  Mesh& mesh() { return m_data->mesh; }

  const Point_2& get_point(const Vertex_index& vertex_index) const {
    return m_data->mesh.point(vertex_index);
  }

  void set_point (const Vertex_index& vertex_index, const Point_2& point)
  {
    m_data->mesh.point(vertex_index) = point;
  }

  void set_last_event_time(const Vertex_index& vertex_index, const FT time) {
    m_data->v_time_map[vertex_index] = time;
  }

  const FT last_event_time(const Vertex_index& vertex_index) const {
    return m_data->v_time_map[vertex_index];
  }

  Vertex_index prev (const Vertex_index& vertex_index) const
  {
    return m_data->mesh.source(m_data->mesh.halfedge(vertex_index));
  }
  Vertex_index next (const Vertex_index& vertex_index) const
  {
    return m_data->mesh.target(m_data->mesh.next(m_data->mesh.halfedge(vertex_index)));
  }

  Face_index face (const Vertex_index& vertex_index) const
  {
    Face_index out = m_data->mesh.face (m_data->mesh.halfedge(vertex_index));
    if (out == Face_index())
      out = m_data->mesh.face (m_data->mesh.opposite(m_data->mesh.halfedge(vertex_index)));
    CGAL_assertion (out != Face_index());
    return out;
  }

  std::pair<Face_index, Face_index> faces (const Vertex_index& vertex_index) const
  {
    for (Halfedge_index hi : halfedges_around_target (halfedge(vertex_index, m_data->mesh), m_data->mesh))
      if (has_iedge (m_data->mesh.edge(hi)))
        return std::make_pair (m_data->mesh.face (hi),
                               m_data->mesh.face (m_data->mesh.opposite(hi)));
    CGAL_assertion_msg (false, "No constrained edge found");
    return std::make_pair (Face_index(), Face_index());
  }

  Point_2 point_2 (const Vertex_index& vertex_index, FT time) const
  { return m_data->mesh.point(vertex_index) + time * m_data->direction[vertex_index]; }

  Point_3 point_3 (const Vertex_index& vertex_index, FT time) const
  { return to_3d (point_2 (vertex_index, time)); }

  Segment_2 segment_2 (const Edge_index& edge_index, FT time) const
  {
    return Segment_2 (point_2 (m_data->mesh.source (m_data->mesh.halfedge(edge_index)), time),
                      point_2 (m_data->mesh.target (m_data->mesh.halfedge(edge_index)), time));
  }

  Segment_3 segment_3 (const Edge_index& edge_index, FT time) const
  {
    return Segment_3 (point_3 (m_data->mesh.source (m_data->mesh.halfedge(edge_index)), time),
                      point_3 (m_data->mesh.target (m_data->mesh.halfedge(edge_index)), time));
  }

  void set_iedge (const Vertex_index& a, const Vertex_index& b,
                  const IEdge& iedge) const
  {
    Halfedge_index hi = m_data->mesh.halfedge (a, b);
    CGAL_assertion (hi != Halfedge_index());
    Edge_index ei = m_data->mesh.edge(hi);
    m_data->e_iedge_map[ei] = iedge;
  }

  void set_ivertex (const Vertex_index& vertex,
                    const IVertex& ivertex) const
  {
    m_data->v_ivertex_map[vertex] = ivertex;
  }

  void set_iedge (const Vertex_index& vertex,
                  const IEdge& iedge) const
  {
    m_data->v_iedge_map[vertex] = iedge;
  }

  void set_iedge (const Edge_index& edge,
                  const IEdge& iedge) const
  {
    m_data->e_iedge_map[edge] = iedge;
  }

  const IEdge& iedge (const Edge_index& edge_index) const
  {
    return m_data->e_iedge_map[edge_index];
  }

  const IEdge& iedge (const Vertex_index& vertex_index) const
  {
    return m_data->v_iedge_map[vertex_index];
  }

  const IVertex& ivertex (const Vertex_index& vertex_index) const
  {
    return m_data->v_ivertex_map[vertex_index];
  }

  bool has_iedge (const Edge_index& edge_index) const
  {
    return (m_data->e_iedge_map[edge_index] != Intersection_graph::null_iedge());
  }
  bool has_iedge (const Vertex_index& vertex_index) const
  {
    return (m_data->v_iedge_map[vertex_index] != Intersection_graph::null_iedge());
  }
  bool has_ivertex (const Vertex_index& vertex_index) const
  {
    return (m_data->v_ivertex_map[vertex_index] != Intersection_graph::null_ivertex());
  }

  const Vector_2& direction (const Vertex_index& vertex_index) const { return m_data->direction[vertex_index]; }
  Vector_2& direction (const Vertex_index& vertex_index) { return m_data->direction[vertex_index]; }
  FT speed (const Vertex_index& vertex_index) const
  { return CGAL::approximate_sqrt (m_data->direction[vertex_index].squared_length()); }

  const KSR::size_t& input (const Face_index& face_index) const { return m_data->input_map[face_index]; }
  KSR::size_t& input (const Face_index& face_index) { return m_data->input_map[face_index]; }

  const bool is_original(const Vertex_index& vertex_index) const { return m_data->v_original_map[vertex_index]; }

  const unsigned int& k (const Face_index& face_index) const { return m_data->k_map[face_index]; }
  unsigned int& k (const Face_index& face_index) { return m_data->k_map[face_index]; }

  bool is_active (const Vertex_index& vertex_index) const { return m_data->v_active_map[vertex_index]; }
  void set_active (const Vertex_index& vertex_index, bool value) { m_data->v_active_map[vertex_index] = value; }

  bool is_frozen (const Vertex_index& vertex_index) const
  {
    return (m_data->direction[vertex_index] == CGAL::NULL_VECTOR);
  }

  const std::set<IEdge>& iedges() const { return m_data->iedges; }
  std::set<IEdge>& iedges() { return m_data->iedges; }

  Point_2 to_2d (const Point_3& point) const
  {
    return m_data->plane.to_2d (point);
  }

  Line_2 to_2d (const Line_3& line) const
  {
    return Line_2 (m_data->plane.to_2d(line.point()),
                   m_data->plane.to_2d(line.point() + line.to_vector()));
  }

  Segment_2 to_2d (const Segment_3& segment) const
  {
    return Segment_2 (m_data->plane.to_2d(segment.source()),
                      m_data->plane.to_2d(segment.target()));
  }

  Vector_3 to_3d (const Vector_2& vec) const
  {
    return Vector_3 (m_data->plane.to_3d (Point_2(0,0)),
                     m_data->plane.to_3d (Point_2(0,0) + vec));
  }

  Point_3 to_3d (const Point_2& point) const { return m_data->plane.to_3d (point); }

  Edge_index edge (const Vertex_index& v0, const Vertex_index& v1)
  {
    // std::cout << int(v0) << " : " << int(v1) << std::endl;
    // std::cout << int(m_data->mesh.halfedge(v0, v1)) << std::endl;
    // std::cout << int(m_data->mesh.halfedge(v1, v0)) << std::endl;
    return m_data->mesh.edge (m_data->mesh.halfedge (v0, v1));
  }

  Edge_index add_edge (const Vertex_index& v0, const Vertex_index& v1,
                       const IEdge& iedge)
  {
    Edge_index out = m_data->mesh.edge (m_data->mesh.add_edge(v0,v1));
    m_data->e_iedge_map[out] = iedge;
    return out;
  }

  Vertex_index add_vertex (const Point_2& point)
  {
    return m_data->mesh.add_vertex(point);
  }

  Vertex_index duplicate_vertex (const Vertex_index& v)
  {
    Vertex_index vi = m_data->mesh.add_vertex (m_data->mesh.point(v));
    m_data->direction[vi] = m_data->direction[v];
    m_data->v_ivertex_map[vi] = m_data->v_ivertex_map[v];
    m_data->v_iedge_map[vi] = m_data->v_iedge_map[v];
    return vi;
  }

  void remove_vertex (const Vertex_index& v)
  {
    m_data->mesh.remove_vertex(v);
  }

  Edge_index split_vertex (const Vertex_index& ei)
  {
    return m_data->mesh.edge
      (CGAL::Euler::split_vertex (m_data->mesh.halfedge (ei),
                                  m_data->mesh.opposite(m_data->mesh.next(m_data->mesh.halfedge (ei))),
                                  m_data->mesh));
  }

  Vertex_index split_edge (const Vertex_index& v0, const Vertex_index& v1)
  {
    return m_data->mesh.target
      (CGAL::Euler::split_edge (m_data->mesh.halfedge (v0, v1), m_data->mesh));
  }


};

template<typename Kernel>
bool operator==(
  const Support_plane<Kernel>& a,
  const Support_plane<Kernel>& b) {

  using FT = typename Kernel::FT;

  const auto& planea = a.plane();
  const auto& planeb = b.plane();

  const auto va = planea.orthogonal_vector();
  const auto vb = planeb.orthogonal_vector();

  const FT sq_dist_to_plane = CGAL::squared_distance(planea.point(), planeb);

  const FT ptol = KSR::point_tolerance<FT>();
  const FT vtol = KSR::vector_tolerance<FT>();
  const FT sq_ptol = ptol * ptol;

  // Are the planes parallel?
  if (CGAL::abs(va * vb) < vtol) return false;
  // Are the planes coplanar?
  return (sq_dist_to_plane < sq_ptol);
}

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_SUPPORT_LINE_H
