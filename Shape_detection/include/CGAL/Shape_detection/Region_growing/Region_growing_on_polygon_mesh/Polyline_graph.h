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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_POLYLINE_GRAPH_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_POLYLINE_GRAPH_H

#include <CGAL/license/Shape_detection.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>

namespace CGAL {
namespace Shape_detection {
namespace Polygon_mesh {

  /*!
    \ingroup PkgShapeDetectionRGOnMesh

    \brief Polygon mesh edges connected into a graph.

    This class returns all edges, which form polylines splitting the polygon mesh
    being a `PolygonMesh` into planar regions.

    \tparam GeomTraits
    a model of `Kernel`

    \tparam PolygonMesh
    a model of `FaceListGraph`

    \tparam FaceRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator` and
    value type is the face type of a polygon mesh

    \tparam EdgeRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator` and
    value type is the edge type of a polygon mesh

    \tparam VertexToPointMap
    a model of `ReadablePropertyMap` whose key type is the vertex type of a polygon mesh and
    value type is `Kernel::Point_3`

    \cgalModels `NeighborQuery`
  */
  template<
  typename GeomTraits,
  typename PolygonMesh,
  typename FaceRange = typename PolygonMesh::Face_range,
  typename EdgeRange = typename PolygonMesh::Edge_range,
  typename VertexToPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type>
  class Polyline_graph {

  private:
    struct PEdge {
      std::size_t index = std::size_t(-1);
      std::set<std::size_t> sneighbors;
      std::set<std::size_t> tneighbors;
      std::pair<int, int> regions;
    };

  public:
    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Face_graph = PolygonMesh;
    using Face_range = FaceRange;
    using Edge_range = EdgeRange;
    using Vertex_to_point_map = VertexToPointMap;
    /// \endcond

    /*!
      a model of `ConstRange` whose iterator type is `RandomAccessIterator` and
      value type is the internal edge type.
    */
    using Segment_range = std::vector<PEdge>;

    /*!
      a model of `ReadablePropertyMap` whose key type is the value type of `Segment_range`
      and value type is `Kernel::Segment_3`.
    */
    using Segment_map = internal::Polyline_graph_segment_map<PEdge, Face_graph, Edge_range, Vertex_to_point_map>;

    /// @}

  private:
    using Face_to_region_map = internal::Item_to_region_index_map;
    using Face_to_index_map = internal::Item_to_index_property_map<Face_range>;
    using Edge_to_index_map = internal::Item_to_index_property_map<Edge_range>;

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param pmesh
      an instance of a `PolygonMesh` that represents a polygon mesh

      \param regions
      a set of regions, where each region is an `std::vector<std::size_t>` and
      contains indices of all faces in a polygon mesh been identified as the same plane.

      \param vertex_to_point_map
      an instance of `VertexToPointMap` that maps a polygon mesh
      vertex to `Kernel::Point_3`

      \pre `faces(pmesh).size() > 0`
      \pre `edges(pmesh).size() > 0`
      \pre `regions.size() > 0`
    */
    Polyline_graph(
      const PolygonMesh& pmesh,
      const std::vector< std::vector<std::size_t> >& regions,
      const VertexToPointMap vertex_to_point_map) :
    m_face_graph(pmesh),
    m_regions(regions),
    m_face_range(faces(m_face_graph)),
    m_edge_range(edges(m_face_graph)),
    m_vertex_to_point_map(vertex_to_point_map),
    m_face_to_region_map(m_face_range, regions),
    m_face_to_index_map(m_face_range),
    m_edge_to_index_map(m_edge_range),
    m_segment_map(m_face_graph, m_edge_range, m_vertex_to_point_map) {

      CGAL_precondition(m_face_range.size() > 0);
      CGAL_precondition(m_edge_range.size() > 0);
      CGAL_precondition(regions.size() > 0);
      build_graph();
    }

    /// \cond SKIP_IN_MANUAL
    void build_graph() {

      clear();
      int r1 = -1, r2 = -1;
      std::vector<std::size_t> pedge_map(
        m_edge_range.size(), std::size_t(-1));
      for (const auto& edge : m_edge_range) {
        std::tie(r1, r2) = get_regions(edge);
        if (r1 == r2) continue;
        add_graph_edge(edge, r1, r2, pedge_map);
      }

      for (std::size_t i = 0; i < m_pedges.size(); ++i) {
        auto& pedge = m_pedges[i];
        CGAL_precondition(pedge.regions.first != pedge.regions.second);
        CGAL_precondition(pedge.index != std::size_t(-1));
        CGAL_precondition(pedge.index < m_edge_range.size());

        const auto& edge = *(m_edge_range.begin() + pedge.index);
        const auto s = source(edge, m_face_graph);
        const auto t = target(edge, m_face_graph);

        add_vertex_neighbors(s, i, pedge_map, pedge.sneighbors);
        add_vertex_neighbors(t, i, pedge_map, pedge.tneighbors);
      }
    }
    /// \endcond

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief implements `NeighborQuery::operator()()`.

      This operator retrieves indices of all edges from `segment_range()`,
      which are neighbors of the edge with the index `query_index`.
      These indices are returned in `neighbors`.

      \param query_index
      index of the query edge

      \param neighbors
      indices of edges, which are neighbors of the query edge

      \pre `query_index < segment_range().size()`
    */
    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      neighbors.clear();
      CGAL_precondition(query_index < segment_range().size());
      const auto& pedge = m_pedges[query_index];
      for (const std::size_t sneighbor : pedge.sneighbors)
        neighbors.push_back(sneighbor);
      for (const std::size_t tneighbor : pedge.tneighbors)
        neighbors.push_back(tneighbor);
    }

    /*!
      \brief returns an instance of `Segment_range` to access edges, which
      form polylines
    */
    const Segment_range& segment_range() const {
      return m_pedges;
    }

    /*!
      \brief returns an instance of `Segment_map` that maps an edge from `segment_range()`
      to `Kernel::Segment_3`.
    */
    const Segment_map& segment_map() const {
      return m_segment_map;
    }

    /*!
      \brief returns indices of all edges from `segment_range()`,
      which are connected to the source vertex of the edge
      with the index `query_index`.

      \pre `query_index < segment_range().size()`
    */
    const std::set<std::size_t>& source_neighbors(
      const std::size_t query_index) const {

      CGAL_precondition(query_index < segment_range().size());
      return m_pedges[query_index].sneighbors;
    }

    /*!
      \brief returns indices of all edges from `segment_range()`,
      which are connected to the target vertex of the edge
      with the index `query_index`.

      \pre `query_index < segment_range().size()`
    */
    const std::set<std::size_t>& target_neighbors(
      const std::size_t query_index) const {

      CGAL_precondition(query_index < segment_range().size());
      return m_pedges[query_index].tneighbors;
    }

    /*!
      \brief returns the edge index of the original polygon mesh edge, which
      corresponds to the edge from `segment_range()` with the index `query_index`.

      \pre `query_index < segment_range().size()`
    */
    std::size_t edge_index(const std::size_t query_index) const {

      CGAL_precondition(query_index < segment_range().size());
      return m_pedges[query_index].index;
    }

    /*!
      \brief returns indices of the two distinct regions, which are shared by the edge
      from `segment_range()` with the index `query_index`.

      \pre `query_index < segment_range().size()`
    */
    const std::pair<int, int>& edge_regions(
      const std::size_t query_index) const {

      CGAL_precondition(query_index < segment_range().size());
      return m_pedges[query_index].regions;
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    void clear() {
      m_pedges.clear();
    }

    void release_memory() {
      m_pedges.shrink_to_fit();
    }
    /// \endcond

  private:
    const Face_graph& m_face_graph;
    const std::vector< std::vector<std::size_t> >& m_regions;
    const Face_range m_face_range;
    const Edge_range m_edge_range;
    const Vertex_to_point_map m_vertex_to_point_map;
    const Face_to_region_map m_face_to_region_map;
    const Face_to_index_map m_face_to_index_map;
    const Edge_to_index_map m_edge_to_index_map;
    const Segment_map m_segment_map;
    Segment_range m_pedges;

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
      std::vector<std::size_t>& pedge_map) {

      PEdge pedge;
      CGAL_precondition(region1 != region2);
      const std::size_t ei = get(m_edge_to_index_map, edge);
      CGAL_precondition(ei != std::size_t(-1));

      pedge.index = ei;
      pedge.regions.first = region1;
      pedge.regions.second = region2;
      CGAL_precondition(pedge.index < pedge_map.size());
      pedge_map[pedge.index] = m_pedges.size();
      m_pedges.push_back(pedge);
    }

    template<typename VertexType>
    void add_vertex_neighbors(
      const VertexType& vertex, const std::size_t curr_pe,
      const std::vector<std::size_t>& pedge_map,
      std::set<std::size_t>& neighbors) const {

      const auto query_hedge = halfedge(vertex, m_face_graph);
      const auto hedges = halfedges_around_target(query_hedge, m_face_graph);
      CGAL_precondition(hedges.size() > 0);
      for (const auto& hedge : hedges) {
        const auto e = edge(hedge, m_face_graph);
        const std::size_t ei = get(m_edge_to_index_map, e);
        CGAL_precondition(ei < pedge_map.size());
        const std::size_t pe = pedge_map[ei];
        if (pe == std::size_t(-1)) continue;
        if (pe == curr_pe) continue;
        CGAL_precondition(pe < m_pedges.size());
        neighbors.insert(pe);
      }
    }
  };

} // namespace Polygon_mesh
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_POLYLINE_GRAPH_H
