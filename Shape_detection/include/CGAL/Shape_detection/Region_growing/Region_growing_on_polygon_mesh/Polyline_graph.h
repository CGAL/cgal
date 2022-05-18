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

    using edge_iterator = typename boost::graph_traits<PolygonMesh>::edge_iterator;
    using Edge_range = Iterator_range<typename edge_iterator>;

    struct PEdge {
      std::size_t index = std::size_t(-1);
      edge_descriptor ed;
      std::vector<std::size_t> neighbors;
      std::pair<std::size_t, std::size_t> regions;
    };

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
             a corresponding planar region id

      \param np a sequence of \ref bgl_namedparameters "Named Parameters"
        among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{vertex_point_map}
          \cgalParamDescription{an instance of `VertexPointMap` that maps a polygon mesh
          vertex to `Kernel::Point_3`}
          \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
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
      : Polyline_graph(pmesh, edges(pmesh), face_to_region_map, np) {}

    /*!
      \brief initializes all internal data structures.
      \tparam FaceToRegionMap
      a model of `ReadablePropertyMap` whose key type is `face_descriptor` of the `PolygonMesh`
      and value type is `std::size_t`

      \tparam NamedParameters
      a sequence of optional \ref bgl_namedparameters "Named Parameters"

      \param pmesh a polygon mesh

      \param edge_range contains all edges in `pmesh` to be considered in the graph

      \param face_to_region_map maps each face of `pmesh` to
             a corresponding planar region id

      \param np a sequence of \ref bgl_namedparameters "Named Parameters"
        among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{vertex_point_map}
          \cgalParamDescription{an instance of `VertexPointMap` that maps a polygon mesh
          vertex to `Kernel::Point_3`}
          \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `faces(pmesh).size() > 0`
      \pre `edges(pmesh).size() > 0`
    */
    template<typename FaceToRegionMap,
      typename NamedParameters = parameters::Default_named_parameters>
      Polyline_graph(
        const PolygonMesh& pmesh,
        const Edge_range edge_range,
        FaceToRegionMap face_to_region_map,
        const NamedParameters& np = parameters::default_values()) {

      clear();

      typedef typename boost::property_map<PolygonMesh, CGAL::dynamic_edge_property_t<std::size_t> >::const_type EdgeIndexMap;
      EdgeIndexMap eimap = get(CGAL::dynamic_edge_property_t<std::size_t>(), pmesh);

      // collect edges either on the boundary or having two different incident regions
      for (const auto& edge : edge_range)
      {
        halfedge_descriptor h1 = halfedge(edge, pmesh),
          h2 = opposite(h1, pmesh);

        face_descriptor f1 = face(h1, pmesh),
          f2 = face(h2, pmesh);

        std::size_t r1 = -1, r2 = -1;

        if (f1 != boost::graph_traits<PolygonMesh>::null_face())
          r1 = get(face_to_region_map, f1);
        if (f2 != boost::graph_traits<PolygonMesh>::null_face())
          r2 = get(face_to_region_map, f2);

        if (r1 == r2)
          put(eimap, edge, std::size_t(-1));
        else
          add_graph_edge(edge, r1, r2, eimap);
      }

      // build adjacency between edges
      typedef typename boost::property_map<PolygonMesh, CGAL::dynamic_vertex_property_t<bool> >::const_type VisitedVertexMap;
      VisitedVertexMap visited_vertices = get(CGAL::dynamic_vertex_property_t<bool>(), pmesh);
      for (std::size_t i = 0; i < m_pedges.size(); ++i)
      {
        put(visited_vertices, source(m_pedges[i].ed, pmesh), false);
        put(visited_vertices, target(m_pedges[i].ed, pmesh), false);
      }

      for (std::size_t i = 0; i < m_pedges.size(); ++i)
      {
        auto& pedge = m_pedges[i];
        CGAL_precondition(pedge.regions.first != pedge.regions.second);
        CGAL_precondition(pedge.index != std::size_t(-1));

        std::array<vertex_descriptor, 2> vrts = { source(pedge.ed, pmesh), target(pedge.ed, pmesh) };

        for (int k = 0; k < 2; ++k)
        {
          if (!get(visited_vertices, vrts[k]))
          {
            add_vertex_neighbors(vrts[k], pmesh, eimap);
            put(visited_vertices, vrts[k], true);
          }
        }
      }
    }
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
      std::vector<std::size_t>& neighbors) const
    {
      neighbors.clear();
      CGAL_precondition(query_index < segment_range().size());
      const auto& pedge = m_pedges[query_index];
      neighbors=pedge.neighbors;
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
      which are connected to a vertex of the edge
      with the index `query_index`.

      \pre `query_index < segment_range().size()`
    */
    const std::vector<std::size_t>& neighbors(
      const std::size_t query_index) const
    {
      CGAL_precondition(query_index < segment_range().size());
      return m_pedges[query_index].neighbors;
    }

    /*!
      \brief returns the edge descriptor which
      corresponds to the edge from `segment_range()` with the index `query_index`.

      \pre `query_index < segment_range().size()`
    */
    edge_descriptor descriptor(const std::size_t query_index) const {
      CGAL_precondition(query_index < segment_range().size());
      return m_pedges[query_index].ed;
    }

    /*!
      \brief returns indices of the two distinct regions, which are shared by the edge
      from `segment_range()` with the index `query_index`.

      \pre `query_index < segment_range().size()`
    */
    const std::pair<std::size_t, std::size_t>& edge_regions(
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
      edge_descriptor edge, const std::size_t region1, const std::size_t region2,
      EdgeIndexMap eimap)
    {
      PEdge pedge;
      CGAL_precondition(region1 != region2);
      const std::size_t ei = m_pedges.size();
      put(eimap, edge, ei);

      pedge.index = ei;
      pedge.ed = edge;
      pedge.regions = {region1, region2};
      m_pedges.push_back(pedge);
    }

    template<typename EdgeIndexMap>
    void add_vertex_neighbors(
      const vertex_descriptor vertex,
      const PolygonMesh& pmesh,
      EdgeIndexMap eimap)
    {
      std::vector<std::size_t> nei;
      for (const auto& hedge : halfedges_around_target(vertex, pmesh))
      {
        const auto e = edge(hedge, pmesh);
        const std::size_t ei = get(eimap, e);
        if (ei == std::size_t(-1)) continue;
        if (nei.size()==2) return;
        nei.push_back(ei);
      }
      if (nei.size()==2)
      {
        m_pedges[nei[0]].neighbors.push_back(nei[1]);
        m_pedges[nei[1]].neighbors.push_back(nei[0]);
      }
    }
  };

} // namespace Polygon_mesh
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_POLYLINE_GRAPH_H
