// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : exude_mesh_3 function definition.
//******************************************************************************

#ifndef CGAL_EXUDE_MESH_3_H
#define CGAL_EXUDE_MESH_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Mesh_3/sliver_criteria.h>
#include <CGAL/Mesh_3/Slivers_exuder.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/Mesh_3/parameters_defaults.h>
#include <CGAL/Named_function_parameters.h>

namespace CGAL {
/*!
 * @ingroup PkgMesh3Functions
 *
 * The function `exude_mesh_3()` performs a sliver exudation process on a Delaunay mesh.
 *
 * The sliver exudation process consists in optimizing the weights of vertices
 * of the weighted Delaunay triangulation in such a way that slivers disappear and
 * the quality of the mesh improves.
 *
 * @warning This optimizer modifies the weight of vertices of the triangulation and,
 * if called, must be the last optimizer to be called. If the mesh is refined after
 * this optimization has been performed, all improvements will be lost.
 *
 * @tparam C3T3 a model of the concept `MeshComplex_3InTriangulation_3`.
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param c3t3 the initial mesh that will be modified by the algorithm to represent the final optimized mesh.
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{time_limit}
 *     \cgalParamDescription{is used to set up, in seconds, a CPU time limit after which the optimization process
 *                           is stopped. This time is measured using the `Real_timer` class. The default value is
 *                           0 and means that there is no time limit.}
 *     \cgalParamType{`double`}
 *     \cgalParamPrecondition{`time_limit >= 0`}
 *     \cgalParamDefault{0}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{sliver_bound}
 *     \cgalParamDescription{is designed to give, in degrees, a targeted lower bound on dihedral angles of mesh cells.
 *                           The exudation process considers in turn all the mesh cells that have a smallest dihedral
 *                           angle less than `sliver_bound` and tries to make them disappear by weighting their vertices.
 *                           The optimization process stops when every cell in the mesh achieves this quality. The
 *                           default value is 0 and means that there is no targeted bound: the exuder then runs as long
 *                           as it can improve the smallest dihedral angles of the set of cells incident to some vertices.}
 *     \cgalParamType{`double`}
 *     \cgalParamPrecondition{`0<= sliver_bound <= 180`}
 *     \cgalParamDefault{0}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return a value of type `CGAL::Mesh_optimization_return_code` which is:
 * <UL>
 * <LI>`CGAL::BOUND_REACHED` when the targeted bound for the smallest dihedral angle in the mesh is reached.
 * <LI>`CGAL::TIME_LIMIT_REACHED` when the time limit is reached.
 * <LI>`CGAL::CANT_IMPROVE_ANYMORE` when exudation process stops because it can no longer improve
 * the smallest dihedral angle of the set of cells incident to some vertex in the mesh.
 * </UL>
 *
 * \cgalHeading{Example}
 *
 * \code{.cpp}
 * // Exude without sliver_bound, using at most 10s CPU time
 * exude_mesh_3(c3t3,
 *              parameters::time_limit(10));
 * \endcode
 *
 * \sa `CGAL::Mesh_optimization_return_code`
 * \sa `CGAL::make_mesh_3()`
 * \sa `CGAL::refine_mesh_3()`
 * \sa `CGAL::perturb_mesh_3()`
 * \sa `CGAL::lloyd_optimize_mesh_3()`
 * \sa `CGAL::odt_optimize_mesh_3()`
 *
 */
template<typename C3T3, typename CGAL_NP_TEMPLATE_PARAMETERS>
Mesh_optimization_return_code exude_mesh_3(C3T3& c3t3,const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  double time_limit=choose_parameter(get_parameter(np,internal_np::maximum_running_time),parameters::default_values_for_mesh_3::time_limit);
  double sliver_bound= choose_parameter(get_parameter(np,internal_np::lower_sliver_bound),parameters::default_values_for_mesh_3::exude_sliver_bound);
  return exude_mesh_3_impl(c3t3,time_limit,sliver_bound);
}

#ifndef CGAL_NO_DEPRECATED_CODE
template<typename C3T3>
CGAL_DEPRECATED
Mesh_optimization_return_code exude_mesh_3(C3T3& c3t3, double time_limit, double sliver_bound = 0)
{
  return exude_mesh_3(c3t3, CGAL::parameters::time_limit(time_limit).sliver_bound(sliver_bound));
}
#endif
#ifndef DOXYGEN_RUNNING
// Overload handling parameters passed with operator=
template<typename C3T3,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
Mesh_optimization_return_code exude_mesh_3(C3T3& c3t3,
                                           const CGAL_NP_CLASS_1&  np1,
                                           const CGAL_NP_CLASS_2&  np2,
                                           const NP& ... nps)
{
  return exude_mesh_3(c3t3,internal_np::combine_named_parameters(np1, np2, nps...));
}

template <typename C3T3>
Mesh_optimization_return_code
exude_mesh_3_impl(C3T3& c3t3,
                  const double time_limit,
                  const double sliver_bound)
{
  typedef typename C3T3::Triangulation Tr;
  typedef Mesh_3::Min_dihedral_angle_criterion<Tr> Sc;
  //typedef Mesh_3::Radius_radio_criterion<Tr> Sc;
  typedef typename Mesh_3::Slivers_exuder<C3T3, Sc> Exuder;

  // Create exuder
  Sc criterion(sliver_bound, c3t3.triangulation());
  Exuder exuder(c3t3, criterion);

  // Set time_limit
  exuder.set_time_limit(time_limit);

  // Launch exudation
  return exuder();
}
#endif //DOXYGEN_RUNNING

} //namespace CGAL

#endif // CGAL_EXUDE_MESH_3_H
