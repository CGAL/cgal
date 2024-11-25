// Copyright (c) 2024  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MAKE_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H
#define CGAL_MAKE_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H

#include <CGAL/license/Constrained_triangulation_3.h>

#include <CGAL/Conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_vertex_base_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>

namespace CGAL {

/*!
 * \ingroup PkgCT_3Functions
 * \brief creates a 3D constrained Delaunay triangulation conforming to the faces of a polygon mesh.
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
 * \tparam Triangulation An instance of the `CGAL::Conforming_constrained_Delaunay_triangulation_3` class template
 *   (or `CGAL::Default`).
 *           - Its `Traits` type must be a model of `ConformingConstrainedDelaunayTriangulationTraits_3`,
 *           - Its point type must be constructible from the point type of the polygon mesh,
 *           - its `Vertex` type must be a model of `ConformingConstrainedDelaunayTriangulationVertexBase_3`, and
 *           - its `Cell` type must be a model of `ConformingConstrainedDelaunayTriangulationCellBase_3`.
 * \tparam PolygonMesh a model of `FaceListGraph`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * If `Triangulation` is `CGAL::Default`, the geometric traits `Traits` is deduced from the polygon mesh type
 * `PolygonMesh` and the named parameters `NamedParameters`. Then, the default conforming constrained Delaunay
 * triangulation is `CGAL::Conforming_constrained_Delaunay_triangulation_3<Traits>`.
 *
 * \param mesh The polygon mesh representing the constraints.
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 *
 * \cgalNamedParamsBegin
 *
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `mesh`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `%Traits::Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, mesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *
 *                     must be available in `PolygonMesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *
 *     \cgalParamType{`Traits`}
 *     \cgalParamDefault{the default constructed traits object `Traits{}`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{face_patch_map}
 *    \cgalParamDescription{a property map associating a patch identifier to each face of `mesh`}
 *    \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
 *
 *                   as key type and with any value type that is a *regular* type (model of `Regular`)}
 *   \cgalParamExtra{If this parameter is omitted, each face of the mesh is considered as a separate patch.}
 *   \cgalParamExtra{Faces with the same patch identifier are considered as part of the same surface patch.}
 *   \cgalParamNEnd
 *
 * \cgalNamedParamsEnd
 *
 * \pre `mesh` must not have self-intersections:
 *      \link CGAL::Polygon_mesh_processing::does_self_intersect
 *      `CGAL::Polygon_mesh_processing::does_self_intersect(mesh, np) == false`
 *      \endlink
 */
template <typename Triangulation = CGAL::Default,
          typename PolygonMesh,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
auto make_conforming_constrained_Delaunay_triangulation_3(const PolygonMesh& mesh,
                                               const CGAL_NP_CLASS& np = parameters::default_values())
{
  using Mesh_geom_traits = typename GetGeomTraits<PolygonMesh, CGAL_NP_CLASS>::type;
  using CDT = typename CGAL::Default::Get<Triangulation,
                                          Conforming_constrained_Delaunay_triangulation_3<Mesh_geom_traits>>::type;
  CDT cdt(mesh, np);
  auto remeshing_cdt{std::move(cdt).convert_for_remeshing()};
  static_assert(std::is_same_v<decltype(remeshing_cdt), CDT>);
  return remeshing_cdt;
}

/*!
 * \ingroup PkgCT_3Functions
 * \brief creates a 3D constrained Delaunay triangulation conforming to the faces of a polygon soup.
 *
 * The polygon soup represents the polygonal constraints that will be enforced during the triangulation process.
 *
    * By default, each face of the polygon soup is considered as a polygonal constraint for the triangulation. The
    * named parameter `face_patch_map` can be used to describe bigger polygonal constraints, possibly with holes. If
    * used, the argument of that parameter must be a property map that maps each face of the polygon soup to a patch
    * identifier. Faces with the same patch identifier are considered as part of the same surface patch. Each of those
    * surface patches (defined as the union of the faces with a given patch id) is supposed to be a polygon or a
    * polygon with holes, with coplanar vertices (or almost coplanar up to the precision of the number type used).
 *
 * The generated triangulation will be constrained to conform to the faces of the polygon soup, or to the surface patches
 * described by the `face_patch_map` property map if provided.
 *
 * \tparam Triangulation An instance of the `CGAL::Conforming_constrained_Delaunay_triangulation_3` class template
 *   (or `CGAL::Default`).
 *           - Its `Traits` type must be a model of `ConformingConstrainedDelaunayTriangulationTraits_3`,
 *           - Its point type must be constructible from the point type of the polygon soup,
 *           - its `Vertex` type must be a model of `ConformingConstrainedDelaunayTriangulationVertexBase_3`, and
 *           - its `Cell` type must be a model of `ConformingConstrainedDelaunayTriangulationCellBase_3`.
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type
 * \tparam PolygonRange a model of the concept `RandomAccessContainer` whose value type is a model of the concept
 *    `RandomAccessContainer` whose value type is `std::size_t`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * If `Triangulation` is `CGAL::Default`, the geometric traits `Traits` is deduced from the point type
 * in `PointRange` and the named parameters `NamedParameters`. Then, the default conforming constrained
 * Delaunay triangulation is `CGAL::Conforming_constrained_Delaunay_triangulation_3<Traits>`.
 *
 * \param points a range of points representing the vertices of the polygon soup
 * \param polygons a range of ranges of indices representing the faces of the polygon soup
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *
 *   \cgalParamNBegin{point_map}
 *     \cgalParamDescription{a property map associating points to the elements of the range `points`}
 *     \cgalParamType{a model of `ReadablePropertyMap` whose value type is a point type convertible to the point type}
 *     \cgalParamDefault{`CGAL::Identity_property_map`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{`Traits`}
 *     \cgalParamDefault{the default constructed traits object `Traits{}`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{face_patch_map}
 *    \cgalParamDescription{a property map associating a patch identifier to each face of `soup`}
 *    \cgalParamType{a class model of `ReadWritePropertyMap` with `std::size_t`
 *                   as key type and with any value type that is a *regular* type (model of `Regular`)}
 *   \cgalParamExtra{If this parameter is omitted, each face of the polygon soup is considered as a separate patch.}
 *   \cgalParamExtra{Faces with the same patch identifier are considered as part of the same surface patch.}
 *   \cgalParamNEnd
 *
 * \cgalNamedParamsEnd
 *
 * \pre the polygon soup must not have self-intersections.
 */
template <typename Triangulation = CGAL::Default,
          typename PointRange,
          typename PolygonRange,
          typename NamedParameters = parameters::Default_named_parameters>
auto make_conforming_constrained_Delaunay_triangulation_3(const PointRange& points,
                                               const PolygonRange& polygons,
                                               const NamedParameters& np = parameters::default_values())
{
  using PointRange_value_type = CGAL::cpp20::remove_cvref_t<decltype(*points.begin())>;
  auto point_map = parameters::choose_parameter(parameters::get_parameter(np, internal_np::point_map),
                                               CGAL::Identity_property_map<PointRange_value_type>{});
  auto get_geom_traits_type = [&]() {
    auto geom_traits_np = parameters::get_parameter(np, internal_np::geom_traits);
    using Geom_traits_np_type = decltype(geom_traits_np);
    if constexpr (!std::is_same_v<Geom_traits_np_type, internal_np::Param_not_found>) {
      return Geom_traits_np_type{};
    } else {
      using Point = CGAL::cpp20::remove_cvref_t<decltype(get(point_map, *points.begin()))>;
      using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
      return Kernel{};
    }
  };

  using Geom_traits = decltype(get_geom_traits_type());
  using CDT =
      typename CGAL::Default::Get<Triangulation, Conforming_constrained_Delaunay_triangulation_3<Geom_traits>>::type;
  CDT cdt(points, polygons, np);
  auto remeshing_cdt{std::move(cdt).convert_for_remeshing()};
  static_assert(std::is_same_v<decltype(remeshing_cdt), CDT>);
  return remeshing_cdt;
}

} // end namespace CGAL

#endif // CGAL_MAKE_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H
