// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_POLYGON_MESH_PROCESSING_GET_BORDER_H
#define CGAL_POLYGON_MESH_PROCESSING_GET_BORDER_H

#include <CGAL/license/Polygon_mesh_processing/core.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/Polygon_mesh_processing/border.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/boost/graph/border.h>"
#include <CGAL/Installation/internal/deprecation_warning.h>

#include <CGAL/boost/graph/border.h>

namespace CGAL {
namespace Polygon_mesh_processing {

  /*!
  * \ingroup PMPDeprecated
  * \deprecated This function is deprecated since \cgal 6.2. Users should use instead  `CGAL::border_halfedges()`.
  *
  * \brief collects the border halfedges of a surface patch defined as a face range.
  *
  * For each returned halfedge `h`, `opposite(h, pmesh)` belongs to a face of the patch,
  * but `face(h, pmesh)` does not belong to the patch.
  *
  * @tparam PolygonMesh model of `HalfedgeGraph`
  * @tparam FaceRange a model of `Range` with value type `boost::graph_traits<PolygonMesh>::%face_descriptor`.
  * @tparam HalfedgeOutputIterator model of `OutputIterator`
     holding `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
     for patch border
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param pmesh the polygon mesh to which the faces in `face_range` belong
  * @param face_range the range of faces defining the patch whose border halfedges
  *                   are collected
  * @param out the output iterator that collects the border halfedges of the patch,
  *            seen from outside.
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{face_index_map}
  *     \cgalParamDescription{a property map associating to each face of `pmesh` a unique index between `0` and `num_faces(pmesh) - 1`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
  *                    as key type and `std::size_t` as value type}
  *     \cgalParamDefault{an automatically indexed internal map}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @returns `out`
  *
  * @see `extract_boundary_cycles()`
  */
  template<typename PolygonMesh
         , typename FaceRange
         , typename HalfedgeOutputIterator
         , typename NamedParameters = parameters::Default_named_parameters>
  CGAL_DEPRECATED HalfedgeOutputIterator border_halfedges(const FaceRange& face_range,
                                                          const PolygonMesh& pmesh,
                                                          HalfedgeOutputIterator out,
                                                          const NamedParameters& np = parameters::default_values())
  {
    return CGAL::border_halfedges(face_range, pmesh, out, np);
  }

  template<typename PolygonMesh
         , typename HalfedgeOutputIterator>
  CGAL_DEPRECATED HalfedgeOutputIterator border_halfedges(const PolygonMesh& pmesh,
                                                          HalfedgeOutputIterator out)
  {
    return CGAL::border_halfedges(pmesh, out);
  }

  // counts the number of connected components of the boundary of the mesh.
  // This function is deprecated since \cgal 6.2. Users should use instead  `CGAL::number_of_borders()`.
  //
  // @tparam PolygonMesh model of `HalfedgeGraph`.
  //
  // @param pmesh the polygon mesh to which `face_range` belong
  //
  template<typename PolygonMesh>
  CGAL_DEPRECATED unsigned int number_of_borders(const PolygonMesh& pmesh)
  {
    return CGAL::number_of_borders(pmesh);
  }

  /// @ingroup PMPDeprecated
  /// \deprecated This function is deprecated since \cgal 6.2. Users should use instead  `CGAL::extract_boundary_cycles()`.
  ///
  /// extracts boundary cycles as a list of halfedges, with one halfedge per border.
  ///
  /// @tparam PolygonMesh a model of `HalfedgeListGraph`
  /// @tparam OutputIterator a model of `OutputIterator` holding objects of type
  ///   `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
  ///
  /// @param pm a polygon mesh
  /// @param out an output iterator where the border halfedges will be put
  ///
  /// @see `border_halfedges()`
  template <typename PolygonMesh, typename OutputIterator>
  CGAL_DEPRECATED OutputIterator extract_boundary_cycles(const PolygonMesh& pm,
                                                         OutputIterator out)
  {
    return CGAL::extract_boundary_cycles(pm, out);
  }

} // end of namespace Polygon_mesh_processing
} // end of namespace CGAL


#endif //CGAL_POLYGON_MESH_PROCESSING_GET_BORDER_H
