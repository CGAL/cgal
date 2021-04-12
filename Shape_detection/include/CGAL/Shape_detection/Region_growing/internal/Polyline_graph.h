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
#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {

  template<
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
    using Face_graph = PolygonMesh;
    using Face_range = FaceRange;
    using Vertex_range = VertexRange;
    using Vertex_to_point_map = VertexToPointMap;

    using Point_range = std::vector<PVertex>;
    using Point_map   = internal::Polyline_graph_point_map<PVertex, Vertex_range, Vertex_to_point_map>;

  private:
    using Face_to_region_map = internal::Item_to_region_index_map;
    using Face_to_index_map = internal::Item_to_index_property_map<Face_range>;
    using Vertex_to_index_map = internal::Item_to_index_property_map<Vertex_range>;

  public:
    template<typename NamedParameters>
    Polyline_graph_points(
      const PolygonMesh& pmesh,
      const std::vector< std::vector<std::size_t> >& regions,
      const NamedParameters& np) :
    m_face_graph(pmesh),
    m_regions(regions),
    m_face_range(faces(m_face_graph)),
    m_vertex_range(vertices(m_face_graph)),
    m_vertex_to_point_map(parameters::choose_parameter(parameters::get_parameter(
      np, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, pmesh))),
    m_face_to_region_map(m_face_range, regions),
    m_face_to_index_map(m_face_range),
    m_vertex_to_index_map(m_vertex_range),
    m_point_map(m_vertex_range, m_vertex_to_point_map) {

      CGAL_precondition(m_face_range.size() > 0);
      CGAL_precondition(m_vertex_range.size() > 0);
      CGAL_precondition(regions.size() > 0);
      build_graph();
    }

    void build_graph() {

      clear();
      int r1 = -1, r2 = -1;
      std::vector<std::size_t> pvertex_map(
        m_vertex_range.size(), std::size_t(-1));
      for (const auto& edge : edges(m_face_graph)) {
        std::tie(r1, r2) = get_regions(edge);
        if (r1 == r2) continue;
        add_graph_edge(edge, r1, r2, pvertex_map);
      }
    }

    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      neighbors.clear();
      CGAL_precondition(query_index < m_pvertices.size());
      const auto& pvertex = m_pvertices[query_index];
      for (const std::size_t neighbor : pvertex.neighbors)
        neighbors.push_back(neighbor);
    }

    const Point_range& point_range() const {
      return m_pvertices;
    }

    const Point_map& point_map() const {
      return m_point_map;
    }

    void clear() {
      m_pvertices.clear();
    }

    void release_memory() {
      m_pvertices.shrink_to_fit();
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
    Point_range m_pvertices;

    template<typename EdgeType>
    std::pair<int, int> get_regions(const EdgeType& edge) const {

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
      return std::make_pair(r1, r2);
    }

    template<typename EdgeType>
    void add_graph_edge(
      const EdgeType& edge, const int region1, const int region2,
      std::vector<std::size_t>& pvertex_map) {

      const auto s = source(edge, m_face_graph);
      const auto t = target(edge, m_face_graph);
      const std::size_t vsource = get(m_vertex_to_index_map, s);
      const std::size_t vtarget = get(m_vertex_to_index_map, t);
      CGAL_precondition(vsource != std::size_t(-1));
      CGAL_precondition(vtarget != std::size_t(-1));

      CGAL_precondition(vsource != vtarget);
      CGAL_precondition(region1 != region2);
      CGAL_precondition(
        pvertex_map.size() == m_vertex_range.size());

      CGAL_precondition(vsource < pvertex_map.size());
      if (pvertex_map[vsource] == std::size_t(-1)) { // add new pvertex
        pvertex_map[vsource] = m_pvertices.size();
        m_pvertices.push_back(PVertex());
        m_pvertices.back().index = vsource;
      }

      CGAL_precondition(vtarget < pvertex_map.size());
      if (pvertex_map[vtarget] == std::size_t(-1)) { // add new pvertex
        pvertex_map[vtarget] = m_pvertices.size();
        m_pvertices.push_back(PVertex());
        m_pvertices.back().index = vtarget;
      }

      // Update pvertex info.
      const std::size_t pv1 = pvertex_map[vsource];
      const std::size_t pv2 = pvertex_map[vtarget];

      CGAL_precondition(pv1 != std::size_t(-1) && pv1 < m_pvertices.size());
      CGAL_precondition(pv2 != std::size_t(-1) && pv2 < m_pvertices.size());

      auto& pvertex1 = m_pvertices[pv1];
      auto& pvertex2 = m_pvertices[pv2];

      pvertex1.neighbors.insert(pv2);
      pvertex1.regions.insert(region1);
      pvertex1.regions.insert(region2);

      pvertex2.neighbors.insert(pv1);
      pvertex2.regions.insert(region1);
      pvertex2.regions.insert(region2);
    }
  };

} // namespace internal
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_POLYLINE_GRAPH_POINTS_H
