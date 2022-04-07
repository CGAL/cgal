// Copyright (c) 2018 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Florent Lafarge, Simon Giraudot, Thien Hoang, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_LEAST_SQUARES_PLANE_FIT_SORTING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_LEAST_SQUARES_PLANE_FIT_SORTING_H

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <vector>

// Boost includes.
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>

// Face graph includes.
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Eigen_diagonalize_traits.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/utils.h>
#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>

namespace CGAL {
namespace Shape_detection {
namespace Polygon_mesh {

  /*!
    \ingroup PkgShapeDetectionRGOnMesh

    \brief Sorting of polygon mesh faces with respect to the local plane fit quality.

    Indices of faces in a polygon mesh are sorted with respect to the quality of the
    least squares plane fit applied to the vertices of neighboring faces of each face.

    \tparam GeomTraits
    must be a model of `Kernel`.

    \tparam PolygonMesh
    must be a model of `FaceListGraph`.

    \tparam NeighborQuery
    must be a model of `NeighborQuery`.

    \tparam FaceRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator` and
    value type is the face type of a polygon mesh.

    \tparam VertexToPointMap
    must be an `LvaluePropertyMap` whose key type is the vertex type of a polygon mesh and
    value type is `Kernel::Point_3`.
  */
  template<
  typename GeomTraits,
  typename PolygonMesh,
  typename NeighborQuery,
  typename FaceRange = typename PolygonMesh::Face_range,
  typename VertexToPointMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type>
  class Least_squares_plane_fit_sorting {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Face_graph = PolygonMesh;
    using Neighbor_query = NeighborQuery;
    using Face_range = FaceRange;
    using Vertex_to_point_map = VertexToPointMap;
    using Seed_map = internal::Seed_property_map;
    /// \endcond

    #ifdef DOXYGEN_RUNNING
      /*!
        an `LvaluePropertyMap` whose key and value type is `std::size_t`.
        This map provides an access to the ordered indices of polygon mesh faces.
      */
      typedef unspecified_type Seed_map;
    #endif

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param pmesh
      an instance of `PolygonMesh` that represents a polygon mesh

      \param neighbor_query
      an instance of `NeighborQuery` that is used internally to
      access face's neighbors

      \param vertex_to_point_map
      an instance of `VertexToPointMap` that maps a polygon mesh
      vertex to `Kernel::Point_3`

      \pre `faces(pmesh).size() > 0`
    */
    Least_squares_plane_fit_sorting(
      const PolygonMesh& pmesh,
      NeighborQuery& neighbor_query,
      const VertexToPointMap vertex_to_point_map = VertexToPointMap()) :
    m_face_graph(pmesh),
    m_neighbor_query(neighbor_query),
    m_face_range(faces(m_face_graph)),
    m_vertex_to_point_map(vertex_to_point_map),
    m_to_local_converter() {

      CGAL_precondition(m_face_range.size() > 0);

      m_order.resize(m_face_range.size());
      for (std::size_t i = 0; i < m_face_range.size(); ++i)
        m_order[i] = i;
      m_scores.resize(m_face_range.size());
    }

    /// @}

    /// \name Sorting
    /// @{

    /*!
      \brief sorts indices of polygon mesh faces.
    */
    void sort() {

      compute_scores();
      CGAL_postcondition(m_scores.size() > 0);

      Compare_scores cmp(m_scores);
      std::sort(m_order.begin(), m_order.end(), cmp);
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief returns an instance of `Seed_map` to access the ordered indices
      of polygon mesh faces.
    */
    Seed_map seed_map() {
      return Seed_map(m_order);
    }

    /// @}

  private:

    // Types.
    using Local_traits = Exact_predicates_inexact_constructions_kernel;
    using Local_FT = typename Local_traits::FT;
    using Local_point_3 = typename Local_traits::Point_3;
    using Local_plane_3 = typename Local_traits::Plane_3;
    using To_local_converter = Cartesian_converter<Traits, Local_traits>;
    using Compare_scores = internal::Compare_scores<Local_FT>;

    // Functions.
    void compute_scores() {

      std::vector<std::size_t> neighbors;
      std::vector<Local_point_3> points;

      for (std::size_t i = 0; i < m_face_range.size(); ++i) {

        neighbors.clear();
        m_neighbor_query(i, neighbors);
        neighbors.push_back(i);

        points.clear();
        for (std::size_t j = 0; j < neighbors.size(); ++j) {
          CGAL_precondition(neighbors[j] < m_face_range.size());

          const auto face = *(m_face_range.begin() + neighbors[j]);
          const auto hedge = halfedge(face, m_face_graph);

          const auto vertices = vertices_around_face(hedge, m_face_graph);
          for (const auto vertex : vertices) {

            const auto& tmp_point = get(m_vertex_to_point_map, vertex);
            points.push_back(m_to_local_converter(tmp_point));
          }
        }
        CGAL_postcondition(points.size() > 0);

        Local_plane_3 fitted_plane;
        Local_point_3 fitted_centroid;

        m_scores[i] = CGAL::linear_least_squares_fitting_3(
          points.begin(), points.end(),
          fitted_plane, fitted_centroid,
          CGAL::Dimension_tag<0>(),
          Local_traits(),
          CGAL::Eigen_diagonalize_traits<Local_FT, 3>());
      }
    }

    // Fields.
    const Face_graph& m_face_graph;
    Neighbor_query& m_neighbor_query;
    const Face_range m_face_range;
    const Vertex_to_point_map m_vertex_to_point_map;

    std::vector<std::size_t> m_order;
    std::vector<Local_FT> m_scores;

    const To_local_converter m_to_local_converter;
  };

} // namespace Polygon_mesh
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_LEAST_SQUARES_PLANE_FIT_SORTING_H
