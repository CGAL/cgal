// Copyright (c) 2019 GeometryFactory Sarl (France).
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

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>
#include <CGAL/KSR_3/Intersection_graph.h>
#include <CGAL/Surface_mesh.h>

namespace CGAL
{

namespace KSR_3
{

template <typename GeomTraits>
class Support_plane
{
public:
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Triangle_2 Triangle_2;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::Line_3 Line_3;
  typedef typename Kernel::Plane_3 Plane_3;

  typedef KSR_3::Intersection_graph<Kernel> Intersection_graph;
  typedef typename Intersection_graph::Vertex_descriptor Intersection_vertex;
  typedef typename Intersection_graph::Edge_descriptor Intersection_edge;

  typedef CGAL::Surface_mesh<Point_2> Mesh;
  typedef typename Mesh::Vertex_index Vertex_index;
  typedef typename Mesh::Edge_index Edge_index;
  typedef typename Mesh::Halfedge_index Halfedge_index;
  typedef typename Mesh::Face_index Face_index;

  typedef std::tuple<Vertex_index, Edge_index, Face_index> Locate_type;

  typedef typename Mesh::template Property_map<Vertex_index, Vector_2> V_vector_map;
  typedef typename Mesh::template Property_map<Vertex_index, Intersection_vertex> V_intersection_map;
  typedef typename Mesh::template Property_map<Edge_index, Intersection_edge> E_intersection_map;
  typedef typename Mesh::template Property_map<Face_index, KSR::size_t> F_index_map;


private:

  struct Data
  {
    Plane_3 plane;
    Mesh mesh;
    V_vector_map direction;
    V_intersection_map v_intersection;
    E_intersection_map e_intersection;
    F_index_map input;
    std::set<Intersection_edge> intersection_edges;
  };

  std::shared_ptr<Data> m_data;

public:

  Support_plane () { }

  template <typename PointRange>
  Support_plane (const PointRange& points)
    : m_data (new Data())
  {
    // Compute support plane
    Vector_3 normal = CGAL::NULL_VECTOR;

    //Newell's method
    for (std::size_t i = 0; i < points.size(); ++ i)
    {
      const Point_3& pa = points[i];
      const Point_3& pb = points[(i+1) % points.size()];
      FT x = normal.x() + (pa.y()-pb.y())*(pa.z()+pb.z());
      FT y = normal.y() + (pa.z()-pb.z())*(pa.x()+pb.x());
      FT z = normal.z() + (pa.x()-pb.x())*(pa.y()+pb.y());
      normal = Vector_3 (x,y,z);
    }
    CGAL_assertion_msg (normal != CGAL::NULL_VECTOR, "Polygon is flat");

    m_data->plane = Plane_3 (points[0], KSR::normalize(normal));
    m_data->direction = m_data->mesh.template add_property_map<Vertex_index, Vector_2>("v:direction", CGAL::NULL_VECTOR).first;
    m_data->v_intersection = m_data->mesh.template add_property_map<Vertex_index, Intersection_vertex>
      ("v:intersection", Intersection_graph::null_vertex()).first;
    m_data->e_intersection = m_data->mesh.template add_property_map<Edge_index, Intersection_edge>
      ("e:intersection", Intersection_graph::null_edge()).first;

    bool okay;
    std::tie (m_data->input, okay) = m_data->mesh.template add_property_map<Face_index, KSR::size_t>("f:input", KSR::no_element());
    CGAL_assertion(okay);
  }

  const Plane_3& plane() const { return m_data->plane; }

  const Mesh& mesh() const { return m_data->mesh; }
  Mesh& mesh() { return m_data->mesh; }

  Point_2 point (const Vertex_index& vertex_index, FT time) const
  { return m_data->mesh.point(vertex_index) + time * m_data->direction[vertex_index]; }
  
  Point_3 point_3 (const Vertex_index& vertex_index, FT time) const
  { return to_3d (point (vertex_index, time)); }

  Segment_2 segment_2 (const Edge_index& edge_index, FT time) const
  {
    return Segment_2 (point (m_data->mesh.source (m_data->mesh.halfedge(edge_index)), time),
                      point (m_data->mesh.target (m_data->mesh.halfedge(edge_index)), time));
  }

  Segment_3 segment_3 (const Edge_index& edge_index, FT time) const
  {
    return Segment_3 (point_3 (m_data->mesh.source (m_data->mesh.halfedge(edge_index)), time),
                      point_3 (m_data->mesh.target (m_data->mesh.halfedge(edge_index)), time));
  }

  void set_intersection_edge (const Vertex_index& a, const Vertex_index& b,
                              const Intersection_edge& intersection_edge) const
  {
    Halfedge_index hi = m_data->mesh.halfedge (a, b);
    CGAL_assertion (hi != Halfedge_index());
    Edge_index ei = m_data->mesh.edge(hi);
    m_data->e_intersection[ei] = intersection_edge;
  }

  void set_intersection_vertex (const Vertex_index& vertex, 
                                const Intersection_vertex& intersection_vertex) const
  {
    m_data->v_intersection[vertex] = intersection_vertex;
  }

  bool has_intersection_edge (const Edge_index& edge_index) const
  {
    return (m_data->e_intersection[edge_index] != Intersection_graph::null_edge());
  }

  const std::set<Intersection_edge>& intersection_edges() const { return m_data->intersection_edges; }
  std::set<Intersection_edge>& intersection_edges() { return m_data->intersection_edges; }

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

  std::array<Vertex_index, 4>
  add_bbox_polygon (const std::array<Point_2, 4>& points,
                    const std::array<Intersection_vertex, 4>& intersection_vertices)
  {
    std::array<Vertex_index, 4> vertices;
    for (std::size_t i = 0; i < 4; ++ i)
    {
      Vertex_index vi = m_data->mesh.add_vertex(points[i]);
      m_data->v_intersection[vi] = intersection_vertices[i];
      vertices[i] = vi;
    }
    
    Face_index fi = m_data->mesh.add_face (vertices);
    m_data->input[fi] = KSR::no_element();

    return vertices;
  }

  KSR::size_t add_polygon (const std::vector<Point_2>& points, const Point_2& centroid,
                           KSR::size_t input_idx)
  {
    std::vector<Vertex_index> vertices;
    vertices.reserve (points.size());
    for (const Point_2& p : points)
    {
      Vertex_index vi = m_data->mesh.add_vertex(p);
      m_data->direction[vi] = KSR::normalize (Vector_2 (centroid, p));
      vertices.push_back (vi);
    }
    
    Face_index fi = m_data->mesh.add_face (vertices);
    m_data->input[fi] = input_idx;

    return KSR::size_t(fi);
  }

  Edge_index edge (const Vertex_index& v0, const Vertex_index& v1)
  {
    for (Halfedge_index hi : halfedges_around_target (halfedge(v0, m_data->mesh), m_data->mesh))
      if (target(hi, m_data->mesh) == v1)
        return m_data->mesh.edge(hi);
    return Edge_index();
  }

  Edge_index add_edge (const Vertex_index& v0, const Vertex_index& v1,
                       const Intersection_edge& intersection_edge)
  {
    Edge_index out = m_data->mesh.edge (m_data->mesh.add_edge(v0,v1));
    m_data->e_intersection[out] = intersection_edge;
    return out;
  }

  Vertex_index add_vertex (const Point_2& point)
  {
    return m_data->mesh.add_vertex(point);
  }

  Vertex_index split_edge (const Edge_index& ei)
  {
    return m_data->mesh.target (CGAL::Euler::split_edge (m_data->mesh.halfedge (ei), m_data->mesh));
  }

  KSR::vector<Edge_index> intersected_edges (const Segment_2& segment) const
  {
    KSR::vector<Edge_index> out;
    for (Edge_index ei : m_data->mesh.edges())
    {
      Segment_2 seg (m_data->mesh.point (m_data->mesh.source (m_data->mesh.halfedge(ei))),
                     m_data->mesh.point (m_data->mesh.target (m_data->mesh.halfedge(ei))));
      if (CGAL::do_intersect (segment, seg))
        out.push_back (ei);
    }
    return out;
  }


};

template <typename Kernel>
bool operator== (const Support_plane<Kernel>& a, const Support_plane<Kernel>& b)
{
  const typename Kernel::Plane_3& va = a.plane();
  const typename Kernel::Plane_3& vb = b.plane();

  if (CGAL::abs(va.orthogonal_vector() * vb.orthogonal_vector()) < CGAL_KSR_SAME_VECTOR_TOLERANCE)
    return false;

  return (CGAL::approximate_sqrt(CGAL::squared_distance (vb.point(), va)) < CGAL_KSR_SAME_POINT_TOLERANCE);
}


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_SUPPORT_LINE_H
