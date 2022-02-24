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

    \cgalModels `NeighborQuery`
  */
  template<typename PolygonMesh
#ifndef CGAL_NO_DEPRECATED_CODE
  , typename FaceRange = typename PolygonMesh::Face_range
#endif
  >
  class One_ring_neighbor_query
  {
    using Face_to_index_map = typename boost::property_map<PolygonMesh, CGAL::dynamic_face_property_t<std::size_t> >::const_type;
  public:
    /// \cond SKIP_IN_MANUAL
    using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;
    using Face_graph = PolygonMesh;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param pmesh
      an instance of a `PolygonMesh` that represents a polygon mesh

      \pre `faces(pmesh).size() > 0`
    */
    One_ring_neighbor_query(
      const PolygonMesh& pmesh) :
    m_face_graph(pmesh),
    m_face_to_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), pmesh))
    {
      m_face_range.reserve(num_faces(pmesh)); // a bit larger if has garbage
      for (face_descriptor f : faces(pmesh))
      {
        put(m_face_to_index_map, f, m_face_range.size());
        m_face_range.push_back(f);
      }
      CGAL_precondition(m_face_range.size() > 0);
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief implements `NeighborQuery::operator()()`.

      This operator retrieves indices of all faces,
      which are edge-adjacent to the face with the index `query_index`.
      These indices are returned in `neighbors`.

      \param query_index
      index of the query face

      \param neighbors
      indices of faces, which are neighbors of the query face

      \pre `query_index < faces(pmesh).size()`
    */
    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      neighbors.clear();
      CGAL_precondition(query_index < m_face_range.size());
      face_descriptor query_face = m_face_range[query_index];
      const auto query_hedge = halfedge(query_face, m_face_graph);

      for (face_descriptor face : faces_around_face(query_hedge, m_face_graph))
      {
        if (face != boost::graph_traits<PolygonMesh>::null_face())
        {
          const std::size_t face_index = get(m_face_to_index_map, face);
          neighbors.push_back(face_index);
        }
      }
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    // A property map that can be used to access indices of the input faces.
    const Face_to_index_map& face_to_index_map() const {
      return m_face_to_index_map;
    }
    /// \endcond

  private:
    const Face_graph& m_face_graph;
    std::vector<face_descriptor> m_face_range;
    Face_to_index_map m_face_to_index_map;
  };

} // namespace Polygon_mesh
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYGON_MESH_ONE_RING_NEIGHBOR_QUERY_H
