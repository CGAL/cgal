// ======================================================================
//
// Copyright (c) 2005-2017 GeometryFactory (France).  All Rights Reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s): Le-Jeng Shiue <Andy.Shiue@gmail.com>
//
// ======================================================================

#ifndef CGAL_SUBDIVISION_METHODS_3_H
#define CGAL_SUBDIVISION_METHODS_3_H

#include <CGAL/basic.h>

#include <CGAL/circulator.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/Subdivision_method_3/subdivision_hosts_3.h>
#include <CGAL/Subdivision_method_3/subdivision_masks_3.h>

namespace CGAL {

/// The namespace containing the subdivision methods.
namespace Subdivision_method_3 {

/*!
\addtogroup PkgSurfaceSubdivisionMethods3Functions

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

namespace parameters = Polygon_mesh_processing::parameters;

// -----------------------------------------------------------------------------

#ifndef DOXYGEN_RUNNING
// backward compatibility
template <class PolygonMesh>
void CatmullClark_subdivision(PolygonMesh& pmesh, int step = 1) {
  PQQ(pmesh, CatmullClark_mask_3<PolygonMesh>(&pmesh, get(vertex_point,pmesh)), step);
}
#endif

/*!
 *
 * applies Catmull-Clark subdivision several times on the control mesh `pmesh`.
 * The geometry of the refined mesh is computed by the geometry policy mask `CatmullClark_mask_3`.
 * This function overwrites the control mesh `pmesh` with the subdivided mesh.
 *
 * @tparam PolygonMesh a model of `MutableFaceGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh a polygon mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *       If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the number of subdivision steps, by default 1.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * \pre `pmesh` must be a triangle mesh.
 **/
template <class PolygonMesh, class NamedParameters>
void CatmullClark_subdivision(PolygonMesh& pmesh, const NamedParameters& np) {
  using boost::choose_param;
  using boost::get_param;
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
  Vpm vpm = choose_param(get_param(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, pmesh));

  unsigned int step = choose_param(get_param(np, internal_np::number_of_iterations), 1);
  CatmullClark_mask_3<PolygonMesh,Vpm> mask(&pmesh, vpm);

  for(unsigned int i = 0; i < step; i++)
    internal::PQQ_1step(pmesh, vpm, mask);
}

// -----------------------------------------------------------------------------

#ifndef DOXYGEN_RUNNING
// backward compatibility
template <class PolygonMesh>
void Loop_subdivision(PolygonMesh& pmesh, int step = 1) {
  PTQ(pmesh, Loop_mask_3<PolygonMesh>(&pmesh, get(vertex_point,pmesh)) , step);
}
#endif

/*!
 *
 * applies Loop subdivision several times on the control mesh `pmesh`.
 * The geometry of the refined mesh is computed by the geometry policy mask `Loop_mask_3`.
 * This function overwrites the control mesh `pmesh` with the subdivided mesh.

 * @tparam PolygonMesh a model of `MutableFaceGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh a polygon mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *       If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the number of subdivision steps, by default 1.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 **/
template <class PolygonMesh, class NamedParameters>
void Loop_subdivision(PolygonMesh& pmesh, const NamedParameters& np) {
  using boost::choose_param;
  using boost::get_param;
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
  Vpm vpm = choose_param(get_param(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, pmesh));

  unsigned int step = choose_param(get_param(np, internal_np::number_of_iterations), 1);
  Loop_mask_3<PolygonMesh,Vpm> mask(&pmesh, vpm);

  for(unsigned int i = 0; i < step; i++)
    internal::PTQ_1step(pmesh, vpm, mask);
}

// -----------------------------------------------------------------------------

#ifndef DOXYGEN_RUNNING
// backward compatibility
template <class PolygonMesh>
void DooSabin_subdivision(PolygonMesh& pmesh, int step = 1) {
  DQQ(pmesh, DooSabin_mask_3<PolygonMesh>(&pmesh, get(vertex_point, pmesh)), step);
}
#endif

/*!
 *
 * applies DooSabin subdivision several times on the control mesh `pmesh`.
 * The geometry of the refined mesh is computed by the geometry policy mask `DooSabin_mask_3`.
 * This function overwrites the control mesh `pmesh` with the subdivided mesh.

 * @tparam PolygonMesh a model of `MutableFaceGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh a polygon mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *       If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the number of subdivision steps, by default 1.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 **/
template <class PolygonMesh, class NamedParameters>
void DooSabin_subdivision(PolygonMesh& pmesh, const NamedParameters& np) {
  using boost::choose_param;
  using boost::get_param;
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
  Vpm vpm = choose_param(get_param(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, pmesh));

  unsigned int step = choose_param(get_param(np, internal_np::number_of_iterations), 1);
  DooSabin_mask_3<PolygonMesh,Vpm> mask(&pmesh, vpm);

  for(unsigned int i = 0; i < step; i++)
    internal::DQQ_1step(pmesh, vpm, mask);
}

// -----------------------------------------------------------------------------

#ifndef DOXYGEN_RUNNING
// backward compatibility
template <class PolygonMesh>
void Sqrt3_subdivision(PolygonMesh& pmesh, int step = 1) {
  Sqrt3(pmesh, Sqrt3_mask_3<PolygonMesh>(&pmesh, get(vertex_point,pmesh)), step);
}
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
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh a polygon mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *       If this parameter is omitted, an internal property map for
 *       `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} the number of subdivision steps, by default 1.
 *     \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * \pre `pmesh` must be a triangle mesh.
 **/
template <class PolygonMesh, class NamedParameters>
void Sqrt3_subdivision(PolygonMesh& pmesh, const NamedParameters& np) {
  using boost::choose_param;
  using boost::get_param;
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
  Vpm vpm = choose_param(get_param(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, pmesh));

  unsigned int step = choose_param(get_param(np, internal_np::number_of_iterations), 1);
  Sqrt3_mask_3<PolygonMesh,Vpm> mask(&pmesh, vpm);

  for(unsigned int i = 0; i < step; i++)
    internal::Sqrt3_1step(pmesh, vpm, mask, (i%2==1));
}

/// @}

} // namespace Subdivision_method_3

} // namespace CGAL

#endif // CGAL_SUBDIVISION_METHODS_3_H
