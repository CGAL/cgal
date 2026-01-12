// Copyright (c) 2021 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Antonio Gomes, Dmitry Anisimov
//

#ifndef CGAL_BARYCENTRIC_BOUNDARY_COORDINATES_3_H
#define CGAL_BARYCENTRIC_BOUNDARY_COORDINATES_3_H

#include <CGAL/license/Barycentric_coordinates_3.h>

#include <CGAL/Barycentric_coordinates_3/internal/utils_3.h>
#include <CGAL/Barycentric_coordinates_3/barycentric_enum_3.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

namespace CGAL {
namespace Barycentric_coordinates {

/*!
  \ingroup PkgBarycentricCoordinates3RefFunctions

  \brief computes boundary barycentric coordinates with respect to a closed convex triangle mesh.

  This function computes boundary barycentric coordinates at a given `query` point
  with respect to the vertices of a simple `polyhedron`, that is one
  coordinate per vertex. The coordinates are stored in a destination range
  beginning at `oi`.

  If `query` is at the vertex, the corresponding coordinate is set to one, while
  all other coordinates are zero. If `query` is on the face, the three corresponding
  coordinates are triangle coordinates, while all other coordinates are set to zero.
  If `query` is not on the boundary, all the coordinates are set to zero.

  \tparam TriangleMesh
  must be a model of the concept `FaceListGraph`.

  \tparam Point
  A model of `GeomTraits::Point_3` with `GeomTraits` being the type of the named parameter `geom_traits`.

  \tparam OutputIterator
  must be an output iterator accepting `GeomTraits::FT` with `GeomTraits` being the type of the named parameter `geom_traits`.

  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param tmesh
  an instance of `TriangleMesh`, which must be a convex simplicial polyhedron

  \param query
  a query point

  \param oi
  the beginning of the destination range with the computed coordinates

  \param np
  an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \return an output iterator to the element in the destination range,
  one past the last coordinate stored + the flag indicating whether the
  query point belongs to the polyhedron boundary

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
                     as key type and `%Point` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `TriangleMesh`.}
    \cgalParamNEnd
    \cgalParamNBegin{geom_traits}
      \cgalParamDescription{an instance of a geometric traits class}
      \cgalParamType{a class model of `BarycentricTraits_3`}
      \cgalParamDefault{a \cgal Kernel deduced from the value type of the vertex-point map, using `CGAL::Kernel_traits`}
      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \pre boost::num_vertices(`tmesh`) >= 4.
  \pre CGAL::is_triangle_mesh(`tmesh`)
  \pre CGAL::is_closed(`tmesh`).
  \pre CGAL::is_strongly_convex_3(`tmesh`).
*/
template<typename TriangleMesh,
         typename Point,
         typename OutputIterator,
         typename NamedParameters = parameters::Default_named_parameters>
std::pair<OutputIterator, bool>
boundary_coordinates_3(const TriangleMesh& tmesh,
                       const Point& query,
                       OutputIterator oi,
                       const NamedParameters& np = parameters::default_values())
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;
  static_assert(std::is_same_v<Geom_traits, typename Kernel_traits<Point>::Kernel>);

  using VPM = typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type;
  VPM vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, tmesh));

  Barycentric_coordinates::internal::Edge_case edge_case;

  if constexpr (std::is_same_v<typename GetGeomTraits<TriangleMesh, NamedParameters>::GT_from_NP, internal_np::Param_not_found>)
    edge_case = internal::locate_wrt_polyhedron(vpm, tmesh, query, oi, Geom_traits());
  else
    edge_case = internal::locate_wrt_polyhedron(vpm, tmesh, query, oi, parameters::get_parameter(np, internal_np::geom_traits));

  if (edge_case == internal::Edge_case::BOUNDARY)
    return { oi, true };
  else {
    internal::get_default(num_vertices(tmesh), oi);
    return { oi, false };
  }
}

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif
