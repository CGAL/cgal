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

    \tparam PolygonMesh
    a model of `FaceListGraph`

    \tparam VertexPointMap
    a model of `ReadablePropertyMap` whose key type is the vertex type of a polygon mesh and
    value type is `Kernel::Point_3`

    \cgalModels `NeighborQuery`
  */
  template<
    typename PolygonMesh,
    typename VertexPointMap = typename property_map_selector<PolygonMesh, CGAL::vertex_point_t>::const_type
    >
  class Polyline_graph {

    using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;
    using edge_descriptor = typename boost::graph_traits<PolygonMesh>::edge_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
    using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;

    struct PEdge {
      std::size_t index = std::size_t(-1);
      edge_descriptor ed;
      std::set<std::size_t> sneighbors;
      std::set<std::size_t> tneighbors;
      std::pair<long, long> regions;
    };

    using Face_to_index_map = internal::Item_to_index_property_map<face_descriptor>;
    using Edge_to_index_map = internal::Item_to_index_property_map<edge_descriptor>;

    struct Transform_pedge {
      const edge_descriptor operator()(const PEdge& pedge) const {
        return pedge.ed;
      }
    };

    using Pedge_iterator = typename std::vector<PEdge>::const_iterator;
    using Transform_iterator = boost::transform_iterator<Transform_pedge, Pedge_iterator>;

  public:
    /// \name Types
    /// @{

    #ifdef DOXYGEN_NS
      /*!
        a model of `ConstRange` whose iterator type is `RandomAccessIterator` and
        value type is `edge_descriptor` of the `PolygonMesh`.
      */
      typedef unspecified_type Segment_range;

      /*!
        a model of `ReadablePropertyMap` whose key type is the value type of `Segment_range`
        and value type is `Kernel::Segment_3`.
      */
      typedef unspecified_type Segment_map;
    #endif

    /// @}

    /// \cond SKIP_IN_MANUAL
    using Segment_range = Iterator_range<Transform_iterator>;
    using Segment_map = Segment_from_edge_descriptor_map<PolygonMesh, VertexPointMap>;
    /// \endcond

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.
      \tparam FaceToRegionMap
      a model of `ReadablePropertyMap` whose key type is `face_descriptor` of the `PolygonMesh`
      and value type is `std::size_t`

      \tparam NamedParameters
      a sequence of optional \ref bgl_namedparameters "Named Parameters"

      \param pmesh a polygon mesh

      \param face_to_region_map maps each face of `pmesh` to
             the corresponding planar regions id

      \param np a sequence of \ref bgl_namedparameters "Named Parameters"
        among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{vertex_point_map}
          \cgalParamDescription{an instance of `VertexPointMap` that maps a polygon mesh
          vertex to `Kernel::Point_3`}
          \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
        \cgalParamNEnd
        \cgalParamNBegin{edge_index_map}
          \cgalParamDescription{a property map associating to each edge of `pmesh` a unique index between `0` and `num_edges(pmesh) - 1`}
          \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%edge_descriptor`
                         as key type and `std::size_t` as value type}
          \cgalParamDefault{an automatically indexed internal map}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `faces(pmesh).size() > 0`
      \pre `edges(pmesh).size() > 0`
    */
    template<typename FaceToRegionMap,
             typename NamedParameters = parameters::Default_named_parameters>
    Polyline_graph(
      const PolygonMesh& pmesh,
      FaceToRegionMap face_to_region_map,
      const NamedParameters& np = parameters::default_values())
      :  m_vpm(parameters::choose_parameter(parameters::get_parameter(
          np, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, pmesh)))
      ,  m_segment_map(&pmesh, m_vpm)
    {
      clear();

      std::vector<std::size_t> pedge_map(
        edges(pmesh).size(), std::size_t(-1));

      typedef typename GetInitializedEdgeIndexMap<PolygonMesh, NamedParameters>::const_type EdgeIndexMap;
      EdgeIndexMap eimap = get_initialized_edge_index_map(pmesh, np);

      // collect edges either on the boundary or having two different incident regions
      for (const auto& edge : edges(pmesh))
      {
        halfedge_descriptor h1 = halfedge(edge, pmesh),
                            h2 = opposite(h1, pmesh);

        face_descriptor f1  = face(h1, pmesh),
                        f2  = face(h2, pmesh);

        std::size_t r1 = -1, r2 = -1;

        if (f1 != boost::graph_traits<PolygonMesh>::null_face())
          r1 = get(face_to_region_map, f1);
        if (f2 != boost::graph_traits<PolygonMesh>::null_face())
          r2 = get(face_to_region_map, f2);

        if (r1 == r2) continue;
        add_graph_edge(edge, r1, r2, pedge_map, eimap);
      }

      // build adjacency between edges
      for (std::size_t i = 0; i < m_pedges.size(); ++i) {
        auto& pedge = m_pedges[i];
        CGAL_precondition(pedge.regions.first != pedge.regions.second);
        CGAL_precondition(pedge.index != std::size_t(-1));

        const edge_descriptor edge = pedge.ed;
        const vertex_descriptor s = source(edge, pmesh);
        const vertex_descriptor t = target(edge, pmesh);

        add_vertex_neighbors(s, i, pedge_map, pedge.sneighbors, pmesh, eimap);
        add_vertex_neighbors(t, i, pedge_map, pedge.tneighbors, pmesh, eimap);
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
    const Segment_range segment_range() const {
      return CGAL::make_range(
        boost::make_transform_iterator(m_pedges.begin(), Transform_pedge()),
        boost::make_transform_iterator(m_pedges.end(), Transform_pedge()));
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
    const std::pair<long, long>& edge_regions(
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
    const VertexPointMap m_vpm;

    const Segment_map m_segment_map;
    std::vector<PEdge> m_pedges;

    template <class EdgeIndexMap>
    void add_graph_edge(
      edge_descriptor edge, const long region1, const long region2,
      std::vector<std::size_t>& pedge_map,
      EdgeIndexMap eimap) {

      PEdge pedge;
      CGAL_precondition(region1 != region2);
      const std::size_t ei = get(eimap, edge);
      CGAL_precondition(ei != std::size_t(-1));

      pedge.index = ei;
      pedge.ed = edge;
      pedge.regions.first = region1;
      pedge.regions.second = region2;
      CGAL_precondition(pedge.index < pedge_map.size());
      pedge_map[pedge.index] = m_pedges.size();
      m_pedges.push_back(pedge);
    }

    template<typename EdgeIndexMap>
    void add_vertex_neighbors(
      const vertex_descriptor vertex, const std::size_t curr_pe,
      const std::vector<std::size_t>& pedge_map,
      std::set<std::size_t>& neighbors,
      const PolygonMesh& pmesh,
      EdgeIndexMap eimap) const {

      const auto query_hedge = halfedge(vertex, pmesh);
      const auto hedges = halfedges_around_target(query_hedge, pmesh);
      CGAL_precondition(hedges.size() > 0);
      for (const auto& hedge : hedges) {
        const auto e = edge(hedge, pmesh);
        const std::size_t ei = get(eimap, e);
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
