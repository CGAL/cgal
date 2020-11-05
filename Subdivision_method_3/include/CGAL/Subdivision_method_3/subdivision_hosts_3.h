// Copyright (c) 2017 GeometryFactory (France).  All Rights Reserved.
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

#ifndef CGAL_SUBDIVISION_HOSTS_3_H
#define CGAL_SUBDIVISION_HOSTS_3_H

#include <CGAL/basic.h>

#include <vector>

#include <CGAL/circulator.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Subdivision_method_3/internal/subdivision_hosts_impl_3.h>

namespace CGAL {

/// The namespace containing the subdivision methods.
namespace Subdivision_method_3 {

/*!
\addtogroup PkgSurfaceSubdivisionMethod3Functions
*/
/// @{

namespace parameters = CGAL::parameters;

// -----------------------------------------------------------------------------

#ifndef DOXYGEN_RUNNING
// backward compatibility
template <class PolygonMesh, class Mask>
void PQQ(PolygonMesh& pmesh, Mask mask, int step = 1) {
  // todo:  static assert that PolygonMesh == Mask::PolygonMesh
  for (int i = 0; i < step; i++)
    internal::PQQ_1step(pmesh, get(vertex_point, pmesh), mask);
}
#endif

/*!
 * applies the PQQ refinement several times on the control mesh `pmesh`.
 * The geometry of the refined mesh is computed by the geometry policy `mask`.
 * This function overwrites the control mesh `pmesh` with the refined mesh.
 *
 * @tparam PolygonMesh a model of `MutableFaceGraph`
 * @tparam Mask a model of `PQQMask_3`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param pmesh a polygon mesh
 * @param mask a geometry policy mask
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
template <class PolygonMesh, class Mask, class NamedParameters>
void PQQ(PolygonMesh& pmesh, Mask mask, const NamedParameters& np) {
  // todo:  static assert that PolygonMesh == Mask::PolygonMesh
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, pmesh));

  unsigned int step = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);
  for (unsigned int i = 0; i < step; i++)
    internal::PQQ_1step(pmesh, vpm, mask);
}

// -----------------------------------------------------------------------------

#ifndef DOXYGEN_RUNNING
// backward compatibility
template <class PolygonMesh,  class Mask>
void PTQ(PolygonMesh& pmesh, Mask mask, int step = 1) {
  for (int i = 0; i < step; i++)
    internal::PTQ_1step(pmesh, get(vertex_point,pmesh), mask);
}
#endif

/*!
 * applies the PTQ refinement several times on the control mesh `pmesh`.
 * The geometry of the refined mesh is computed by the geometry policy `mask`.
 * This function overwrites the control mesh `pmesh` with the refined mesh.
 *
 * @tparam PolygonMesh a model of `MutableFaceGraph`
 * @tparam Mask a model of `PTQMask_3`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param pmesh a polygon mesh
 * @param mask a geometry policy mask
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
template <class PolygonMesh, class Mask, class NamedParameters>
void PTQ(PolygonMesh& pmesh, Mask mask, const NamedParameters& np) {
  // todo:  static assert that PolygonMesh == Mask::PolygonMesh
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, pmesh));

  unsigned int step = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);
  for (unsigned int i = 0; i < step; i++)
    internal::PTQ_1step(pmesh, vpm, mask);
}

// -----------------------------------------------------------------------------

#ifndef DOXYGEN_RUNNING
// backward compatibility
template <class PolygonMesh, class Mask>
void DQQ(PolygonMesh& pmesh, Mask mask, int step = 1) {
  for (int i = 0; i < step; ++i) {
    internal::DQQ_1step(pmesh, get(vertex_point, pmesh), mask);
  }
}
#endif

/*!
 * applies the DQQ refinement several times on the control mesh `pmesh`.
 * The geometry of the refined mesh is computed by the geometry policy `mask`.
 * This function overwrites the control mesh `pmesh` with the refined mesh.
 *
 * @tparam PolygonMesh a model of `MutableFaceGraph`
 * @tparam Mask a model of `DQQMask_3`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param pmesh a polygon mesh
 * @param mask a geometry policy mask
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
template <class PolygonMesh, class Mask, class NamedParameters>
void DQQ(PolygonMesh& pmesh, Mask mask, const NamedParameters& np) {
  // todo:  static assert that PolygonMesh == Mask::PolygonMesh
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, pmesh));

  unsigned int step = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);
  for (unsigned int i = 0; i < step; i++)
    internal::DQQ_1step(pmesh, vpm, mask);
}

// -----------------------------------------------------------------------------

#ifndef DOXYGEN_RUNNING
// backward compatibility
template <class PolygonMesh, class Mask>
void Sqrt3(PolygonMesh& pmesh, Mask mask, int step = 1) {
  for (int i = 0; i < step; i++)
    internal::Sqrt3_1step(pmesh, get(vertex_point, pmesh), mask, (i%2==1));
}
#endif

/*!
 * applies the \f$ \sqrt{3}\f$ refinement several times on the control mesh `pmesh`.
 * The geometry of the refined mesh is computed by the geometry policy `mask`.
 * This function overwrites the control mesh `pmesh` with the refined mesh.
 *
 * \attention The border subdivision only happens every second subdivision step
 *            during a <em>single</em> call of this function.
 *
 * @tparam PolygonMesh a model of `MutableFaceGraph`
 * @tparam Mask a model of `Sqrt3Mask_3`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param pmesh a polygon mesh
 * @param mask a geometry policy mask
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
template <class PolygonMesh, class Mask, class NamedParameters>
void Sqrt3(PolygonMesh& pmesh, Mask mask, const NamedParameters& np) {
  // todo:  static assert that PolygonMesh == Mask::PolygonMesh
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type Vpm;
  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, pmesh));

  unsigned int step = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);
  for (unsigned int i = 0; i < step; i++)
    internal::Sqrt3_1step(pmesh, vpm, mask, (i%2==1));
}

/// @}

} // namespace Subdivision_method_3

} // namespace CGAL

#endif // CGAL_SUBDIVISION_HOSTS_3_H
