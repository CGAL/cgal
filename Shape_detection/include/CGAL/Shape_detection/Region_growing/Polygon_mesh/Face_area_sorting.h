// Copyright (c) 2024-2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sébastien Loriot
//


#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_FACE_AREA_SORTING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_FACE_AREA_SORTING_H

#include <CGAL/license/Shape_detection.h>
#include <CGAL/Shape_detection/Region_growing/internal/utils.h>

namespace CGAL {
namespace Shape_detection {
namespace Polygon_mesh {

  /*!
    \ingroup PkgShapeDetectionRGOnMesh

    \brief Sorting of polygon mesh faces with respect to their area.

    `Items` of faces in a polygon mesh are sorted in decreasing area.

    \tparam GeomTraits
    a model of `Kernel`

    \tparam PolygonMesh
    a model of `FaceListGraph`

    \tparam VertexToPointMap
    a model of `ReadablePropertyMap` whose key type is the vertex type of a polygon mesh and
    value type is `Kernel::Point_3`
  */
  template<
  typename GeomTraits,
  typename PolygonMesh,
  typename VertexToPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type>
  class Face_area_sorting
  {
    using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
    using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
    using FT = typename GeomTraits::FT;

  public:
    /// \name Types
    /// @{

    /// Item type.
    using Item = typename boost::graph_traits<PolygonMesh>::face_descriptor;

    /// Seed range.
    using Seed_range = std::vector<Item>;

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam NamedParameters
      a sequence of \ref bgl_namedparameters "Named Parameters"

      \param pmesh
      an instance of `PolygonMesh` that represents a polygon mesh

      \param np
      a sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{vertex_point_map}
          \cgalParamDescription{an instance of `VertexToPointMap` that maps a polygon mesh
          vertex to `Kernel::Point_3`}
          \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
        \cgalParamNEnd
        \cgalParamNBegin{geom_traits}
          \cgalParamDescription{an instance of `GeomTraits`}
          \cgalParamDefault{`GeomTraits()`}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `faces(pmesh).size() > 0`
    */
    template<typename CGAL_NP_TEMPLATE_PARAMETERS>
    Face_area_sorting(
      const PolygonMesh& pmesh,
      const CGAL_NP_CLASS& np = parameters::default_values())
        : m_face_graph(pmesh)
        , m_vpm(parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                                             get_const_property_map(CGAL::vertex_point, pmesh)))
        , m_traits(parameters::choose_parameter<GeomTraits>(parameters::get_parameter(np, internal_np::geom_traits)))
    {
      CGAL_precondition(faces(pmesh).size() > 0);

      m_ordered.resize(faces(pmesh).size());

      std::size_t index = 0;
      for (Item item : faces(pmesh))
        m_ordered[index++] = item;
      m_scores.resize(m_ordered.size(), 0.);
    }

    /*!
      \brief initializes all internal data structures.
      3 Parameter constructor with dummy parameter provided for compatibility with other sorting types.

      \tparam Dummy
      Dummy parameter, not used.

      \tparam NamedParameters
      a sequence of \ref bgl_namedparameters "Named Parameters"

      \param pmesh
      an instance of `PolygonMesh` that represents a polygon mesh

      \param np
      a sequence of \ref bgl_namedparameters "Named Parameters"
      among the ones listed below

      \cgalNamedParamsBegin
        \cgalParamNBegin{vertex_point_map}
          \cgalParamDescription{an instance of `VertexToPointMap` that maps a polygon mesh
          vertex to `Kernel::Point_3`}
          \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
        \cgalParamNEnd
        \cgalParamNBegin{geom_traits}
          \cgalParamDescription{an instance of `GeomTraits`}
          \cgalParamDefault{`GeomTraits()`}
        \cgalParamNEnd
      \cgalNamedParamsEnd

      \pre `faces(pmesh).size() > 0`
    */
    template<typename Dummy, typename CGAL_NP_TEMPLATE_PARAMETERS>
      Face_area_sorting(
        const PolygonMesh& pmesh,
        const Dummy&,
        const CGAL_NP_CLASS& np = parameters::default_values())
      : m_face_graph(pmesh)
      , m_vpm(parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
        get_const_property_map(CGAL::vertex_point, pmesh)))
      , m_traits(parameters::choose_parameter<GeomTraits>(parameters::get_parameter(np, internal_np::geom_traits)))
    {
      CGAL_precondition(faces(pmesh).size() > 0);

      m_ordered.resize(faces(pmesh).size());

      std::size_t index = 0;
      for (Item item : faces(pmesh))
        m_ordered[index++] = item;
      m_scores.resize(m_ordered.size(), 0.);
    }

    /// @}

    /// \name Sorting
    /// @{

    /*!
      \brief sorts `Items` of input faces.
    */
    void sort() {
      compute_scores();
      CGAL_precondition(m_scores.size() > 0);

      auto cmp = [this](const std::size_t i, const std::size_t j)
      {
        CGAL_precondition(i < m_scores.size());
        CGAL_precondition(j < m_scores.size());
        return m_scores[i] > m_scores[j];
      };
      std::vector<std::size_t> order(m_ordered.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), cmp);

      std::vector<Item> tmp(m_ordered.size());
      for (std::size_t i = 0; i < m_ordered.size(); i++)
        tmp[i] = m_ordered[order[i]];

      m_ordered.swap(tmp);
    }
    /// @}

    /// \name Access
    /// @{

    /*!
      \brief returns an instance of `Seed_range` to access the ordered `Items`
      of input faces.
    */
    const Seed_range &ordered() {
      return m_ordered;
    }
    /// @}

  private:
    const PolygonMesh& m_face_graph;
    VertexToPointMap m_vpm;
    const GeomTraits m_traits;
    Seed_range m_ordered;
    std::vector<FT> m_scores;

    void compute_scores()
    {
      auto squared_area = m_traits.compute_squared_area_3_object();
      std::size_t idx = 0;
      for (Item item : m_ordered)
      {
        halfedge_descriptor hd = halfedge(item, m_face_graph);
        std::vector<typename GeomTraits::Point_3> pts;

        for (vertex_descriptor v : vertices_around_face(hd,m_face_graph))
          pts.push_back( get(m_vpm, v) );

        if (pts.size()==3)
          m_scores[idx++] = approximate_sqrt(squared_area(pts[0], pts[1], pts[2]));
        else
        {
          std::vector<typename GeomTraits::Triangle_3> triangles;
          internal::triangulate_face<GeomTraits>(pts, triangles);
          for (const typename GeomTraits::Triangle_3& tr : triangles)
            m_scores[idx] += approximate_sqrt(squared_area(tr[0], tr[1], tr[2]));
          ++idx;
        }
      }
    }
  };

} // namespace Polygon_mesh
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_FACE_AREA_SORTING_H

