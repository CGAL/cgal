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
  typedef typename Intersection_graph::Vertex_descriptor IVertex;
  typedef typename Intersection_graph::Edge_descriptor IEdge;

  typedef CGAL::Surface_mesh<Point_2> Mesh;
  typedef typename Mesh::Vertex_index Vertex_index;
  typedef typename Mesh::Edge_index Edge_index;
  typedef typename Mesh::Halfedge_index Halfedge_index;
  typedef typename Mesh::Face_index Face_index;

  typedef std::tuple<Vertex_index, Edge_index, Face_index> Locate_type;

  typedef typename Mesh::template Property_map<Vertex_index, Vector_2> V_vector_map;
  typedef typename Mesh::template Property_map<Vertex_index, IVertex> V_ivertex_map;
  typedef typename Mesh::template Property_map<Vertex_index, IEdge> V_iedge_map;
  typedef typename Mesh::template Property_map<Vertex_index, bool> V_bool_map;
  typedef typename Mesh::template Property_map<Edge_index, IEdge> E_iedge_map;
  typedef typename Mesh::template Property_map<Face_index, KSR::size_t> F_index_map;
  typedef typename Mesh::template Property_map<Face_index, unsigned int> F_uint_map;


private:

  struct Data
  {
    Plane_3 plane;
    Mesh mesh;
    V_vector_map direction;
    V_ivertex_map v_ivertex_map;
    V_iedge_map v_iedge_map;
    V_bool_map v_active_map;
    E_iedge_map e_iedge_map;
    F_index_map input_map;
    F_uint_map k_map;
    std::set<IEdge> iedges;
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
    m_data->v_ivertex_map = m_data->mesh.template add_property_map<Vertex_index, IVertex>
      ("v:ivertex", Intersection_graph::null_ivertex()).first;
    m_data->v_iedge_map = m_data->mesh.template add_property_map<Vertex_index, IEdge>
      ("v:iedge", Intersection_graph::null_iedge()).first;
    m_data->v_active_map = m_data->mesh.template add_property_map<Vertex_index, bool>
      ("v:active", true).first;
    m_data->e_iedge_map = m_data->mesh.template add_property_map<Edge_index, IEdge>
      ("e:iedge", Intersection_graph::null_iedge()).first;
    m_data->input_map = m_data->mesh.template add_property_map<Face_index, KSR::size_t>
      ("f:input", KSR::no_element()).first;
    m_data->k_map = m_data->mesh.template add_property_map<Face_index, unsigned int>
      ("f:k", 0).first;
  }

  const Plane_3& plane() const { return m_data->plane; }

  const Mesh& mesh() const { return m_data->mesh; }
  Mesh& mesh() { return m_data->mesh; }

  void set_point (const Vertex_index& vertex_index, const Point_2& point)
  {
    m_data->mesh.point(vertex_index) = point;
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

  std::array<Vertex_index, 4>
  add_bbox_polygon (const std::array<Point_2, 4>& points,
                    const std::array<IVertex, 4>& ivertices)
  {
    std::array<Vertex_index, 4> vertices;
    for (std::size_t i = 0; i < 4; ++ i)
    {
      Vertex_index vi = m_data->mesh.add_vertex(points[i]);
      m_data->v_ivertex_map[vi] = ivertices[i];
      vertices[i] = vi;
    }
    
    Face_index fi = m_data->mesh.add_face (vertices);
    m_data->input_map[fi] = KSR::no_element();

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
    m_data->input_map[fi] = input_idx;

    return KSR::size_t(fi);
  }

  Edge_index edge (const Vertex_index& v0, const Vertex_index& v1)
  {
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
