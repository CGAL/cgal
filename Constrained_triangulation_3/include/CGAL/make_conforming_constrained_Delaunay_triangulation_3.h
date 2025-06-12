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

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

namespace CGAL {

/*!
  \addtogroup PkgConstrainedTriangulation3FunctionsPolygonSoupOrMesh
  @{

  Free Functions Template for Creating Conforming Constrained Delaunay Triangulations {#PkgConstrainedTriangulation3Functions}
  ==========================================================================

  The following functions create a 3D conforming constrained Delaunay triangulation
  from either a polygon soup or a polygon mesh.

  Input Data        {#make_conforming_constrained_Delaunay_triangulation_3_input_data}
  ----------

  The input data (polygon mesh or polygon soup) represents the polygonal constraints enforced
  during the triangulation process.

  By default, each face of the input is considered a polygonal constraint for the triangulation. The
  named parameter `plc_face_id` can be used to describe larger polygonal constraints, possibly with holes. If
  used, this parameter must be a property map that associates each face of the input with a patch
  identifier. Faces with the same patch identifier are considered part of the same surface patch. Each of these
  surface patches (defined as the union of the input faces with a given patch identifier) is expected to be a polygon or
  a polygon with holes, with coplanar vertices (or nearly coplanar up to the precision of the number type used).

  The generated triangulation will conform to the faces of the input, or to the surface patches
  described by the `plc_face_id` property map if provided.


  \pre The input data must not be coplanar.

  Template Parameters        {#make_conforming_constrained_Delaunay_triangulation_3_template_parameters}
  -------------------

  For both function templates, the template arguments can be deduced from the function arguments, or defaulted.

  - The first template argument `Triangulation` defaults to `CGAL::Default`, and in that case the
    triangulation type is deduced from the input type and the named parameters (see below).
  - The following one or two template arguments are deduced from the input data: either a polygon mesh type,
    or a polygon soup defined by two types (a sequence of points and a sequence of sequences of indices).
  - The last template argument is the named parameters class, deduced from the function arguments.


  Returned Triangulation Type      {#make_conforming_constrained_Delaunay_triangulation_3_returned_type}
  ---------------------------

  For both functions, the template parameter `Triangulation` defines the type of the triangulation that is created
  and returned by the function.

    - If `Triangulation` is `CGAL::Default`, the geometric traits class type is deduced from the input data and
      the named parameters. If the named parameter `geom_traits` is provided, the traits class is deduced from it.
      Otherwise, the point type of the input data is used to deduce the traits class. Let's call it `Traits`.
      The returned triangulation type is then `CGAL::Conforming_constrained_Delaunay_triangulation_3<Traits>`.
    - Otherwise, `Triangulation` must be a specialization of the `CGAL::Conforming_constrained_Delaunay_triangulation_3`
      class template, with the following requirements:
        - its `Vertex` type must be a model of `ConformingConstrainedDelaunayTriangulationVertexBase_3`, and
        - its `Cell` type must be a model of `ConformingConstrainedDelaunayTriangulationCellBase_3`.

  In both cases, the traits class `Traits` must fulfill the following requirements:
        - It must be a model of the concept `ConformingConstrainedDelaunayTriangulationTraits_3`.
        - It must have a `Point_3` type that is constructible from the point type of the input data.
  @}
*/

/*!
 * \ingroup PkgConstrainedTriangulation3FunctionsPolygonSoupOrMesh
 *
 * \brief creates a 3D constrained Delaunay triangulation conforming to the faces of a polygon mesh.
 *
 * \tparam PolygonMesh a model of `FaceListGraph`
 * \tparam Triangulation an instance of the `CGAL::Conforming_constrained_Delaunay_triangulation_3` class template,
 *                       or `CGAL::Default`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param mesh the polygon mesh representing the constraints
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `mesh`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `%Traits::Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, mesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `PolygonMesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{plc_face_id}
 *    \cgalParamDescription{a property map associating a patch identifier to each facet of `mesh`.
 *                   Each identifier corresponds to a planar surface patch. Each surface
 *                   patch can be composed of several facets of `mesh`, forming a planar polygon.}
 *    \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
 *                   as key type and with any value type that is a *regular* type (model of `Regular`)}
 *   \cgalParamExtra{If this parameter is omitted, each facet of `mesh` is considered a separate patch.
 *                   Faces with the same patch identifier are considered part of the same surface patch.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{return_empty_on_invalid_input}
 *     \cgalParamDescription{a Boolean activating the check of preconditions on the input
 *                           before starting to construct the triangulation.
 *                           If the check is activated and the preconditions are not satisfied,
 *                           an empty triangulation is immediately returned.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *     \cgalParamExtra{If this parameter is set to `false` and the preconditions are not satisfied, the algorithm may throw an exception or crash.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{`Traits` as defined above in the section \ref make_conforming_constrained_Delaunay_triangulation_3_returned_type}
 *     \cgalParamDefault{the default constructed traits object `Traits{}`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return a 3D constrained Delaunay triangulation conforming to the faces of the polygon mesh, of a type
 *   described in the section \ref make_conforming_constrained_Delaunay_triangulation_3_returned_type above.
 *
 * \pre `mesh` must not have self-intersections.
 *      For triangulated surfaces, it can be checked using the function
 *      \link CGAL::Polygon_mesh_processing::does_self_intersect
 *      `CGAL::Polygon_mesh_processing::does_self_intersect(mesh, np) == false`
 *      \endlink
 */
template <typename Triangulation = CGAL::Default,
          typename PolygonMesh,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
auto make_conforming_constrained_Delaunay_triangulation_3(const PolygonMesh &mesh,
                                                          const CGAL_NP_CLASS &np = parameters::default_values())
{
  using Mesh_geom_traits = typename GetGeomTraits<PolygonMesh, CGAL_NP_CLASS>::type;
  using CDT = typename CGAL::Default::Get<Triangulation,
                                          Conforming_constrained_Delaunay_triangulation_3<Mesh_geom_traits>>::type;
  CDT cdt(mesh, np);
  auto remeshing_cdt{std::move(cdt)};
  static_assert(std::is_same_v<decltype(remeshing_cdt), CDT>);
  return remeshing_cdt;
}

/*!
 * \ingroup PkgConstrainedTriangulation3FunctionsPolygonSoupOrMesh
 *
 * \brief creates a 3D constrained Delaunay triangulation conforming to the faces of a polygon soup.
 *
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type
 *         of the polygon soup
 * \tparam PolygonRange a model of the concept `RandomAccessContainer` whose value type is a model of the concept
 *    `RandomAccessContainer` whose value type is `std::size_t`
 * \tparam Triangulation an instance of the `CGAL::Conforming_constrained_Delaunay_triangulation_3` class template,
 *         or `CGAL::Default`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param points a range of points representing the vertices of the polygon soup
 * \param polygons a range of ranges of indices representing the faces of the polygon soup
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *
 *   \cgalParamNBegin{point_map}
 *     \cgalParamDescription{a property map associating points to the elements of the range `points`}
 *     \cgalParamType{a model of `ReadablePropertyMap` whose
 *                    key type is the value type of `PointRange`, and
 *                    value type is convertible to the point type of the geometric traits class}
 *     \cgalParamDefault{`CGAL::Identity_property_map`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{plc_face_id}
 *    \cgalParamDescription{a property map associating a patch identifier to each facet of `soup`.
 *                   Each identifier corresponds to a planar surface patch. Each surface
 *                   patch can be composed of several faces of `soup`, forming a planar polygon.}
 *    \cgalParamType{a class model of `ReadablePropertyMap` with `std::size_t`
 *                   as key type and with any value type that is a *regular* type (model of `Regular`)}
 *   \cgalParamExtra{If this parameter is omitted, each facet of the polygon soup is considered a separate patch.}
 *   \cgalParamExtra{Otherwise facets with the same patch identifier are considered part of the same surface patch.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{return_empty_on_invalid_input}
 *     \cgalParamDescription{a Boolean activating the check of preconditions on the input
 *                           before starting to construct the triangulation.
 *                           If the check is activated and the preconditions are not satisfied,
 *                           an empty triangulation is immediately returned.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *     \cgalParamExtra{If this parameter is set to `false` and the preconditions are not satisfied, the algorithm may throw an exception or crash.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{`Traits` as defined above in the section \ref make_conforming_constrained_Delaunay_triangulation_3_returned_type}
 *     \cgalParamDefault{the default constructed traits object `Traits{}`}
 *   \cgalParamNEnd
 *
 * \cgalNamedParamsEnd
 *
 * \return a 3D constrained Delaunay triangulation conforming to the faces of the polygon soup, of a type
 *   described in the section \ref make_conforming_constrained_Delaunay_triangulation_3_returned_type above.
 *
 * \pre The polygon soup must be free of self-intersections. If the polygon soup is a triangle soup, this is equivalent to:
 *      \link CGAL::Polygon_mesh_processing::does_triangle_soup_self_intersect
 *     `CGAL::Polygon_mesh_processing::does_triangle_soup_self_intersect(points, polygons, np) == false`
 *     \endlink.
 */
template <typename Triangulation = CGAL::Default,
          typename PointRange,
          typename PolygonRange,
          typename NamedParameters = parameters::Default_named_parameters>
auto make_conforming_constrained_Delaunay_triangulation_3(const PointRange &points,
                                                          const PolygonRange &polygons,
                                                          const NamedParameters &np = parameters::default_values())
{
  using Geom_traits = typename GetPolygonGeomTraits<PointRange, PolygonRange, NamedParameters>::type;

  using Default_CDT = Conforming_constrained_Delaunay_triangulation_3<Geom_traits>;
  using CDT = typename CGAL::Default::Get<Triangulation, Default_CDT>::type;
  CDT cdt(points, polygons, np);
  auto remeshing_cdt{std::move(cdt)};
  static_assert(std::is_same_v<decltype(remeshing_cdt), CDT>);
  return remeshing_cdt;
}

} // end namespace CGAL

#endif // CGAL_MAKE_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H
