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

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_ONE_RING_NEIGHBOR_QUERY_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_ONE_RING_NEIGHBOR_QUERY_H

#include <CGAL/license/Shape_detection.h>

// Internal includes.
#include <CGAL/Shape_detection/Region_growing/internal/property_map.h>

namespace CGAL {
namespace Shape_detection {
namespace Polygon_mesh {

  /*!
    \ingroup PkgShapeDetectionRGOnMesh

    \brief Edge-adjacent faces connectivity in a polygon mesh.

    This class returns all faces, which are edge-adjacent to a query face in a
    polygon mesh being a `PolygonMesh`.

    \tparam PolygonMesh
    a model of `FaceListGraph`

    \tparam FaceRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator` and
    value type is the face type of a polygon mesh

    \cgalModels{NeighborQuery}
  */
  template<typename PolygonMesh>
  class One_ring_neighbor_query
  {
    using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;
    using Face_graph = PolygonMesh;
  public:
    /// Item type.
    using Item = typename boost::graph_traits<PolygonMesh>::face_descriptor;
    using Region = std::vector<Item>;

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param pmesh
      an instance of a `PolygonMesh` that represents a polygon mesh

      \pre `faces(pmesh).size() > 0`
    */
    One_ring_neighbor_query(const PolygonMesh& pmesh)
      : m_face_graph(pmesh)
    {}

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief implements `NeighborQuery::operator()()`.

      This operator retrieves all faces,
      which are edge-adjacent to the face `query`.
      These `Items` are returned in `neighbors`.

      \param query
      `Item` of the query face

      \param neighbors
      `Items` of faces, which are neighbors of the query face

      \pre `query_index < faces(pmesh).size()`
    */
    void operator()(
      const Item query,
      std::vector<Item>& neighbors) const {

      neighbors.clear();
      const auto query_hedge = halfedge(query, m_face_graph);

      for (face_descriptor face : faces_around_face(query_hedge, m_face_graph))
      {
        if (face != boost::graph_traits<PolygonMesh>::null_face())
        {
          neighbors.push_back(face);
        }
      }
    }

    /// @}

  private:
    const Face_graph& m_face_graph;
  };

} // namespace Polygon_mesh
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_ONE_RING_NEIGHBOR_QUERY_H
