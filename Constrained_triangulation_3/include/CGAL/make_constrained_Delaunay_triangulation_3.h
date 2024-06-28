// Copyright (c) 2024  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MAKE_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H
#define CGAL_MAKE_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H

#include <CGAL/license/Constrained_triangulation_3.h>

#include <CGAL/Constrained_Delaunay_triangulation_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_vertex_base_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>

namespace CGAL {

/*!
 * \ingroup PkgCT_3Classes
 * \brief The default 3D constrained Delaunay triangulation type.
 *
 * `Default_constrained_Delaunay_triangulation_3` is a metafunction that returns the
 * default 3D constrained Delaunay triangulation type for a given geometric
 * traits class.
 *
 * \tparam Geom_traits a geometric traits class.
 *
 * \return `type` is the default 3D constrained Delaunay triangulation type.
 *
 * \sa default_constrained_Delaunay_triangulation_3_t
*/
template <typename Geom_traits,
          typename Vb = Constrained_Delaunay_triangulation_vertex_base_3<Geom_traits>,
          typename Cb = Constrained_Delaunay_triangulation_cell_base_3<Geom_traits>>
struct Default_constrained_Delaunay_triangulation_3
{
  using Tds = Triangulation_data_structure_3<Vb, Cb>;
  using Tr = Triangulation_3<Geom_traits, Tds>;
  using type = Constrained_Delaunay_triangulation_3<Geom_traits, Tr>;
};

/*!
 * \ingroup PkgCT_3Classes
 * \brief The default 3D constrained Delaunay triangulation type.
 * \tparam Geom_traits a geometric traits class.
 *
 * This alias template names the default 3D constrained Delaunay triangulation
 * type for a given geometric traits class.
 *
 * \sa Default_constrained_Delaunay_triangulation_3
 */
template <typename Geom_traits>
using default_constrained_Delaunay_triangulation_3_t = typename Default_constrained_Delaunay_triangulation_3<Geom_traits>::type;

/*!
 * \ingroup PkgCT_3Functions
 * \brief Create a 3D constrained Delaunay triangulation conforming to the faces of a polygon mesh.
 *
 * The polygon mesh represents the polygonal constraints that will be enforced during the triangulation process.
 *
 * By default, each face of the polygon mesh is considered as a polygonal constraint for the triangulation. The
 * named parameter `face_patch_map` can be used to describe bigger polygonal constraints, possibly with holes. If
 * used, the argument of that parameter must be a property map that maps each face of the polygon mesh to a patch
 * identifier. Faces with the same patch identifier are considered as part of the same surface patch. Each of those
 * surface patches (defined as the union of the mesh faces with a given patch id) is supposed to be a polygon or a
 * polygon with holes, with coplanar vertices (or almost coplanar up to the precision of the number type used).
 *
 * The generated triangulation will be constrained to conform to the faces of the polygon mesh, or to the surface patches
 * described by the `face_patch_map` property map if provided.
 *
 * \tparam Triangulation An instance of the `CGAL::Triangulation_3` class template (or `CGAL::Default`).
 *           - Its `Geom_traits` type must be a model of `ConstrainedDelaunayTriangulationTraits_3`,
 *           - Its point type must be constructible from the point type of the polygon mesh,
 *           - its `Vertex` type must be a model of `ConstrainedDelaunayTriangulationVertexBase_3`, and
 *           - its `Cell` type must be a model of `ConstrainedDelaunayTriangulationCellBase_3`.
 * \tparam PolygonMesh a model of `FaceListGraph`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * If `Triangulation` is `CGAL::Default`, the geometric traits `Geom_traits` is deduced from the polygon mesh type
 * `PolygonMesh` and the named parameters `NamedParameters`. And then the default constrained Delaunay triangulation is
 * `Default_constrained_Delaunay_triangulation_3<Geom_traits>::type`.
 *
 * \param mesh The polygon mesh representing the constraints.
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `mesh`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `%Traits::Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, mesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `PolygonMesh`.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{`Traits`}
 *     \cgalParamDefault{the default constructed traits object `Traits{}`}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{face_patch_map}
 *    \cgalParamDescription{a property map associating a patch identifier to each face of `mesh`}
 *    \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
 *                   as key type and with any value type that is a *regular* type (model of `Regular`)}
 *   \cgalParamExtra{If this parameter is omitted, each face of the mesh is considered as a separate patch.}
 *   \cgalParamExtra{Faces with the same patch identifier are considered as part of the same surface patch.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \pre `mesh` must not have self-intersections:
 *      \link CGAL::Polygon_mesh_processing::does_self_intersect
 *      `CGAL::Polygon_mesh_processing::does_self_intersect(mesh, np) == false`
 *      \endlink
 */
template <typename Triangulation = CGAL::Default,
          typename PolygonMesh,
          typename NamedParams = parameters::Default_named_parameters>
auto make_constrained_Delaunay_triangulation_3(const PolygonMesh& mesh,
                                               const NamedParams& np = parameters::default_values())
{
  using Mesh_geom_traits = typename GetGeomTraits<PolygonMesh, NamedParams>::type;
  using CDT = typename CGAL::Default::Get<Triangulation, default_constrained_Delaunay_triangulation_3_t<Mesh_geom_traits>>::type;
  CDT cdt(mesh, np);
  auto remeshing_cdt{std::move(cdt).convert_for_remeshing()};
  static_assert(std::is_same_v<decltype(remeshing_cdt), CDT>);
  return remeshing_cdt;
}

} // end namespace CGAL

#endif // CGAL_MAKE_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H
