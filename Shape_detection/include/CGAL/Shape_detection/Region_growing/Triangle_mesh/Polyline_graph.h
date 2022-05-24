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
namespace Triangle_mesh {

  /*!
    \ingroup PkgShapeDetectionRGOnMesh

    \brief Triangle mesh edges connected into a graph.

    This class returns all edges, which form polylines splitting the triangle mesh
    being a `TriangleMesh` into planar regions.

    \tparam TriangleMesh
    a model of `FaceListGraph`

    \tparam VertexPointMap
    a model of `ReadablePropertyMap` whose key type is the vertex type of a triangle mesh and
    value type is `Kernel::Point_3`

    \cgalModels `NeighborQuery`
  */
  template<
    typename TriangleMesh,
    typename VertexPointMap = typename property_map_selector<TriangleMesh, CGAL::vertex_point_t>::const_type
    >
  class Polyline_graph {

    using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
    using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
    using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

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
        value type is `edge_descriptor` of the `TriangleMesh`.
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
    using Segment_map = Segment_from_edge_descriptor_map<TriangleMesh, VertexPointMap>;
    /// \endcond

  public:
    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.
      \tparam FaceToRegionMap
      a model of `ReadablePropertyMap` whose key type is `face_descriptor` of the `TriangleMesh`
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
          \cgalParamDescription{an instance of `VertexPointMap` that maps a triangle mesh
          vertex to `Kernel::Point_3`}
          \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `faces(pmesh).size() > 0`
      \pre `edges(pmesh).size() > 0`
    */
    template<typename FaceToRegionMap,
             typename CGAL_NP_TEMPLATE_PARAMETERS>
    Polyline_graph(
      const TriangleMesh& tmesh,
      FaceToRegionMap face_to_region_map,
      const CGAL_NP_CLASS& np = parameters::default_values())
      : Polyline_graph(tmesh, edges(tmesh), face_to_region_map, np)
    {}

    /*!
      \brief initializes all internal data structures.
      \tparam FaceToRegionMap
      a model of `ReadablePropertyMap` whose key type is `face_descriptor` of the `TriangleMesh`
      and value type is `std::size_t`

      \tparam NamedParameters
      a sequence of optional \ref bgl_namedparameters "Named Parameters"

      \tparam EdgeRange a model of `ConstRange` with `edge_descriptor` as iterator value type.

      \param tmesh a triangle mesh

      \param edge_range contains all edges in `tmesh` to be considered in the graph

      \param face_to_region_map maps each face of `tmesh` to
             a corresponding planar region id

      \param np a sequence of \ref bgl_namedparameters "Named Parameters"
        among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{vertex_point_map}
          \cgalParamDescription{an instance of `VertexPointMap` that maps a triangle mesh
          vertex to `Kernel::Point_3`}
          \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `faces(pmesh).size() > 0`
      \pre `edges(pmesh).size() > 0`
    */
    template<typename FaceToRegionMap,
             typename EdgeRange,
             typename NamedParameters = parameters::Default_named_parameters>
    Polyline_graph(
      const TriangleMesh& tmesh,
      const EdgeRange& edge_range,
      FaceToRegionMap face_to_region_map,
      const NamedParameters& np = parameters::default_values())
    :  m_vpm(parameters::choose_parameter(parameters::get_parameter(
              np, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, tmesh)))
    ,  m_segment_map(&tmesh, m_vpm)
    {
      clear();

      typedef typename boost::property_map<TriangleMesh, CGAL::dynamic_edge_property_t<std::size_t> >::const_type EdgeIndexMap;
      EdgeIndexMap eimap = get(CGAL::dynamic_edge_property_t<std::size_t>(), tmesh);

      // init map
      for(edge_descriptor e : edges(tmesh))
        put(eimap, e, std::size_t(-1));

      // collect edges either on the boundary or having two different incident regions
      for (const auto& edge : edge_range)
      {
        halfedge_descriptor h1 = halfedge(edge, tmesh),
          h2 = opposite(h1, tmesh);

        face_descriptor f1 = face(h1, tmesh),
          f2 = face(h2, tmesh);

        std::size_t r1 = -1, r2 = -1;

        if (f1 != boost::graph_traits<TriangleMesh>::null_face())
          r1 = get(face_to_region_map, f1);
        if (f2 != boost::graph_traits<TriangleMesh>::null_face())
          r2 = get(face_to_region_map, f2);

        if (r1 != r2)
          add_graph_edge(edge, r1, r2, eimap);
      }

      // build adjacency between edges
      typedef typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<bool> >::const_type VisitedVertexMap;
      VisitedVertexMap visited_vertices = get(CGAL::dynamic_vertex_property_t<bool>(), tmesh);
      for (std::size_t i = 0; i < m_pedges.size(); ++i)
      {
        put(visited_vertices, source(m_pedges[i].ed, tmesh), false);
        put(visited_vertices, target(m_pedges[i].ed, tmesh), false);
      }

      for (std::size_t i = 0; i < m_pedges.size(); ++i)
      {
        auto& pedge = m_pedges[i];
        CGAL_precondition(pedge.regions.first != pedge.regions.second);
        CGAL_precondition(pedge.index != std::size_t(-1));

        std::array<vertex_descriptor, 2> vrts = { source(pedge.ed, tmesh), target(pedge.ed, tmesh) };

        for (int k = 0; k < 2; ++k)
        {
          if (!get(visited_vertices, vrts[k]))
          {
            add_vertex_neighbors(vrts[k], tmesh, eimap);
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
      const TriangleMesh& pmesh,
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

} // namespace Triangle_mesh
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_TRIANGLE_MESH_POLYLINE_GRAPH_H
