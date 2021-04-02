// Copyright (c) 2005-2017 GeometryFactory (France).  All Rights Reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Le-Jeng Shiue <Andy.Shiue@gmail.com>
//

#ifndef CGAL_SUBDIVISION_METHODS_3_H
#define CGAL_SUBDIVISION_METHODS_3_H

#include <CGAL/basic.h>

#include <CGAL/circulator.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Subdivision_method_3/subdivision_hosts_3.h>
#include <CGAL/Subdivision_method_3/subdivision_masks_3.h>

namespace CGAL {

/// The namespace containing the subdivision methods.
namespace Subdivision_method_3 {

/*!
\addtogroup PkgSurfaceSubdivisionMethod3Functions

A subdivision method recursively refines a coarse mesh and
generates an ever closer approximation to a smooth surface.
`Subdivision_method_3` consists of four subdivision methods
and their refinement hosts. Each refinement host is a template
function of a polygon mesh class and a
geometry policy class. It refines the connectivity of the
control mesh and computes the geometry of the refined mesh.
The geometry computation is dedicated to the custom
geometry policy. A geometry policy consists of functions
that compute the new point based on the subdivision stencil.
A stencil defines the footprint (a submesh of the control mesh)
of a new point.

The four supported refinement hosts are the
primal quadrilateral quadrisection (PQQ),
the primal triangle quadrisection (PTQ),
the dual quadrilateral quadrisection (DQQ),
and the \f$ \sqrt{3}\f$ triangulation.
These refinements are respectively used in
Catmull-Clark, Loop, Doo-Sabin and \f$ \sqrt{3}\f$ subdivisions.

\cgalHeading{Refinement Host}

A refinement host is a template function of
a polygon mesh class and a geometry mask class. It refines
the input polygon mesh, and computes new points through
the geometry masks.
`Subdivision_method_3` supports four refinement hosts:
`PQQ`, `PTQ`, `DQQ` and `Sqrt3`.

\image html RefSchemes.svg

\cgalHeading{Example}

This example program subdivides a polygonal mesh with
Catmull-Clark subdivision.

\cgalExample{Subdivision_method_3/CatmullClark_subdivision.cpp}

\sa `CGAL::CatmullClark_mask_3<PolygonMesh>`
\sa `CGAL::DooSabin_mask_3<PolygonMesh`
\sa `CGAL::Loop_mask_3<PolygonMesh`
\sa `CGAL::Sqrt3_mask_3<PolygonMesh>`
*/
/// @{

namespace parameters = CGAL::parameters;

// -----------------------------------------------------------------------------

#ifndef DOXYGEN_RUNNING
// Backward compatibility
#ifndef CGAL_NO_DEPRECATED_CODE
template <class PolygonMesh>
CGAL_DEPRECATED_MSG("you are using the deprecated API of CatmullClark_subdivision(), please update your code")
void CatmullClark_subdivision(PolygonMesh& pmesh, int step) {
  PQQ(pmesh, CatmullClark_mask_3<PolygonMesh>(&pmesh, get(vertex_point,pmesh)), step);
}
#endif
#endif

/*!
 *
 * applies Catmull-Clark subdivision several times on the control mesh `pmesh`.
 * The geometry of the refined mesh is computed by the geometry policy mask `CatmullClark_mask_3`.
 * This function overwrites the control mesh `pmesh` with the subdivided mesh.
 *
 * @tparam PolygonMesh a model of `MutableFaceGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param pmesh a polygon mesh
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `pmesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{number_of_iterations}
 *     \cgalParamDescription{the number of subdivision steps}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{`1`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \pre `pmesh` must be a triangle mesh.
 **/
template <class PolygonMesh, class NamedParameters>
void CatmullClark_subdivision(PolygonMesh& pmesh, const NamedParameters& np) {
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, pmesh));

  unsigned int step = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);
  CatmullClark_mask_3<PolygonMesh,Vpm> mask(&pmesh, vpm);

  for(unsigned int i = 0; i < step; i++)
    internal::PQQ_1step(pmesh, vpm, mask);
}

template <class PolygonMesh>
void CatmullClark_subdivision(PolygonMesh& pmesh)
{
  CatmullClark_subdivision(pmesh, CGAL::parameters::all_default());
}
// -----------------------------------------------------------------------------

#ifndef DOXYGEN_RUNNING
// backward compatibility
#ifndef CGAL_NO_DEPRECATED_CODE
template <class PolygonMesh>
CGAL_DEPRECATED_MSG("you are using the deprecated API of Loop_subdivision(), please update your code")
void Loop_subdivision(PolygonMesh& pmesh, int step) {
  PTQ(pmesh, Loop_mask_3<PolygonMesh>(&pmesh, get(vertex_point,pmesh)) , step);
}
#endif
#endif

/*!
 *
 * applies Loop subdivision several times on the control mesh `pmesh`.
 * The geometry of the refined mesh is computed by the geometry policy mask `Loop_mask_3`.
 * This function overwrites the control mesh `pmesh` with the subdivided mesh.

 * @tparam PolygonMesh a model of `MutableFaceGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param pmesh a polygon mesh
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `pmesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{number_of_iterations}
 *     \cgalParamDescription{the number of subdivision steps}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{`1`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 **/
template <class PolygonMesh, class NamedParameters>
void Loop_subdivision(PolygonMesh& pmesh, const NamedParameters& np) {
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, pmesh));

  unsigned int step = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);
  Loop_mask_3<PolygonMesh,Vpm> mask(&pmesh, vpm);

  for(unsigned int i = 0; i < step; i++)
    internal::PTQ_1step(pmesh, vpm, mask);
}

template <class PolygonMesh>
void Loop_subdivision(PolygonMesh& pmesh)
{
  Loop_subdivision(pmesh, CGAL::parameters::all_default());
}
// -----------------------------------------------------------------------------

#ifndef DOXYGEN_RUNNING
// backward compatibility
#ifndef CGAL_NO_DEPRECATED_CODE
template <class PolygonMesh>
CGAL_DEPRECATED_MSG("you are using the deprecated API of DooSabin_subdivision(), please update your code")
void DooSabin_subdivision(PolygonMesh& pmesh, int step) {
  DQQ(pmesh, DooSabin_mask_3<PolygonMesh>(&pmesh, get(vertex_point, pmesh)), step);
}
#endif
#endif

/*!
 *
 * applies DooSabin subdivision several times on the control mesh `pmesh`.
 * The geometry of the refined mesh is computed by the geometry policy mask `DooSabin_mask_3`.
 * This function overwrites the control mesh `pmesh` with the subdivided mesh.

 * @tparam PolygonMesh a model of `MutableFaceGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param pmesh a polygon mesh
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `pmesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{number_of_iterations}
 *     \cgalParamDescription{the number of subdivision steps}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{`1`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 **/
template <class PolygonMesh, class NamedParameters>
void DooSabin_subdivision(PolygonMesh& pmesh, const NamedParameters& np) {
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, pmesh));

  unsigned int step = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);
  DooSabin_mask_3<PolygonMesh,Vpm> mask(&pmesh, vpm);

  for(unsigned int i = 0; i < step; i++)
    internal::DQQ_1step(pmesh, vpm, mask);
}

template <class PolygonMesh>
void DooSabin_subdivision(PolygonMesh& pmesh)
{
  DooSabin_subdivision(pmesh, CGAL::parameters::all_default());
}
// -----------------------------------------------------------------------------

#ifndef DOXYGEN_RUNNING
// backward compatibility
#ifndef CGAL_NO_DEPRECATED_CODE
template <class PolygonMesh>
CGAL_DEPRECATED_MSG("you are using the deprecated API of Sqrt3_subdivision(), please update your code")
void Sqrt3_subdivision(PolygonMesh& pmesh, int step) {
  Sqrt3(pmesh, Sqrt3_mask_3<PolygonMesh>(&pmesh, get(vertex_point,pmesh)), step);
}
#endif
#endif

/*!
 *
 * applies \f$ \sqrt{3}\f$-subdivision several times on the control mesh `pmesh`.
 * The geometry of the refined mesh is computed by the geometry policy mask `Sqrt3_mask_3`.
 * This function overwrites the control mesh `pmesh` with the subdivided mesh.
 *
 * \attention The border subdivision only happens every second subdivision step
 *            during a <em>single</em> call of this function.
 *
 * @tparam PolygonMesh a model of `MutableFaceGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param pmesh a polygon mesh
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `pmesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{number_of_iterations}
 *     \cgalParamDescription{the number of subdivision steps}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{`1`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \pre `pmesh` must be a triangle mesh.
 **/
template <class PolygonMesh, class NamedParameters>
void Sqrt3_subdivision(PolygonMesh& pmesh, const NamedParameters& np) {
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, pmesh));

  unsigned int step = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);
  Sqrt3_mask_3<PolygonMesh,Vpm> mask(&pmesh, vpm);

  for(unsigned int i = 0; i < step; i++)
    internal::Sqrt3_1step(pmesh, vpm, mask, (i%2==1));
}

template <class PolygonMesh>
void Sqrt3_subdivision(PolygonMesh& pmesh)
{
  Sqrt3_subdivision(pmesh, CGAL::parameters::all_default());
}
/// @}

} // namespace Subdivision_method_3

} // namespace CGAL

#endif // CGAL_SUBDIVISION_METHODS_3_H
