// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_POLYLINE_GRAPH_POINTS_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_POLYLINE_GRAPH_POINTS_H

#include <CGAL/license/Shape_detection.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/utils.h>
#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {

  template<
  typename GeomTraits,
  typename PolygonMesh,
  typename FaceRange = typename PolygonMesh::Face_range,
  typename VertexRange = typename PolygonMesh::Vertex_range,
  typename VertexToPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type>
  class Polyline_graph_points {

  private:
    struct PVertex {
      std::size_t index = std::size_t(-1);
      std::set<std::size_t> neighbors;
      std::set<std::size_t> regions;
    };

  public:
    using Traits = GeomTraits;
    using Face_graph = PolygonMesh;
    using Face_range = FaceRange;
    using Vertex_range = VertexRange;
    using Vertex_to_point_map = VertexToPointMap;

    using PVertex_range = std::vector<PVertex>;
    using Point_map = internal::Polyline_graph_point_map<PVertex, Vertex_range, Vertex_to_point_map>;

  private:
    using Face_to_region_map = internal::Item_to_region_index_map;
    using Face_to_index_map = internal::Item_to_index_property_map<Face_range>;
    using Vertex_to_index_map = internal::Item_to_index_property_map<Vertex_range>;

  public:
    Polyline_graph_points(
      const PolygonMesh& pmesh,
      const std::vector< std::vector<std::size_t> >& regions,
      const VertexToPointMap vertex_to_point_map) :
    m_face_graph(pmesh),
    m_regions(regions),
    m_face_range(faces(m_face_graph)),
    m_vertex_range(vertices(m_face_graph)),
    m_vertex_to_point_map(vertex_to_point_map),
    m_face_to_region_map(m_face_range, regions),
    m_face_to_index_map(m_face_range),
    m_vertex_to_index_map(m_vertex_range),
    m_point_map(m_vertex_range, m_vertex_to_point_map) {

      CGAL_precondition(m_face_range.size() > 0);
      CGAL_precondition(m_vertex_range.size() > 0);
      build_graph();
    }

    void build_graph() {

      clear();
      std::vector<std::size_t> vertex_map(
        m_vertex_range.size(), std::size_t(-1));
      for (const auto& edge : edges(m_face_graph)) {
        const auto hedge1 = halfedge(edge, m_face_graph);
        const auto hedge2 = opposite(hedge1, m_face_graph);

        const auto face1 = face(hedge1, m_face_graph);
        const auto face2 = face(hedge2, m_face_graph);

        const std::size_t fi1 = get(m_face_to_index_map, face1);
        const std::size_t fi2 = get(m_face_to_index_map, face2);
        CGAL_precondition(fi1 != fi2);

        int r1 = -1, r2 = -1;
        if (fi1 != std::size_t(-1))
          r1 = get(m_face_to_region_map, fi1);
        if (fi2 != std::size_t(-1))
          r2 = get(m_face_to_region_map, fi2);
        if (r1 == r2) continue;
        CGAL_precondition(r1 != r2);
        add_graph_edge(edge, r1, r2, vertex_map);
      }
    }

    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      neighbors.clear();
      CGAL_precondition(query_index < m_vertices.size());
      const auto& vertex = m_vertices[query_index];
      for (const std::size_t neighbor : vertex.neighbors)
        neighbors.push_back(neighbor);
    }

    const PVertex_range& pvertex_range() const {
      return m_vertices;
    }

    const Point_map& point_map() const {
      return m_point_map;
    }

    void clear() {
      m_vertices.clear();
    }

  private:
    const Face_graph& m_face_graph;
    const std::vector< std::vector<std::size_t> >& m_regions;
    const Face_range m_face_range;
    const Vertex_range m_vertex_range;
    const Vertex_to_point_map m_vertex_to_point_map;
    const Face_to_region_map m_face_to_region_map;
    const Face_to_index_map m_face_to_index_map;
    const Vertex_to_index_map m_vertex_to_index_map;
    const Point_map m_point_map;
    PVertex_range m_vertices;

    template<typename EdgeType>
    void add_graph_edge(
      const EdgeType& edge, const int region1, const int region2,
      std::vector<std::size_t>& vertex_map) {

      const auto s = source(edge, m_face_graph);
      const auto t = target(edge, m_face_graph);
      const std::size_t vsource = get(m_vertex_to_index_map, s);
      const std::size_t vtarget = get(m_vertex_to_index_map, t);
      CGAL_precondition(vsource != std::size_t(-1));
      CGAL_precondition(vtarget != std::size_t(-1));

      CGAL_precondition(vsource != vtarget);
      CGAL_precondition(region1 != region2);
      CGAL_precondition(
        vertex_map.size() == m_vertex_range.size());

      CGAL_precondition(vsource < vertex_map.size());
      if (vertex_map[vsource] == std::size_t(-1)) { // add new vertex
        vertex_map[vsource] = m_vertices.size();
        m_vertices.push_back(PVertex());
        m_vertices.back().index = vsource;
      }

      CGAL_precondition(vtarget < vertex_map.size());
      if (vertex_map[vtarget] == std::size_t(-1)) { // add new vertex
        vertex_map[vtarget] = m_vertices.size();
        m_vertices.push_back(PVertex());
        m_vertices.back().index = vtarget;
      }

      // Update vertex info.
      const std::size_t v1 = vertex_map[vsource];
      const std::size_t v2 = vertex_map[vtarget];
      CGAL_precondition(v1 != std::size_t(-1) && v1 < m_vertices.size());
      CGAL_precondition(v2 != std::size_t(-1) && v2 < m_vertices.size());

      auto& vertex1 = m_vertices[v1];
      auto& vertex2 = m_vertices[v2];
      vertex1.neighbors.insert(v2);
      vertex2.neighbors.insert(v1);
      vertex1.regions.insert(region1);
      vertex1.regions.insert(region2);
      vertex2.regions.insert(region1);
      vertex2.regions.insert(region2);
    }
  };

} // namespace internal
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_POLYLINE_GRAPH_POINTS_H
