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
// File Description : lloyd_optimize_mesh_3 function definition.
//******************************************************************************

#ifndef CGAL_LLOYD_OPTIMIZE_MESH_3_H
#define CGAL_LLOYD_OPTIMIZE_MESH_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/Mesh_3/Mesh_global_optimizer.h>
#include <CGAL/Mesh_3/Lloyd_move.h>
#include <CGAL/Mesh_3/Mesh_sizing_field.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/Mesh_3/parameters_defaults.h>
#include <CGAL/Mesh_3/internal/check_weights.h>


namespace CGAL {
/*!
 * \ingroup PkgMesh3Functions
 *
 * The function `lloyd_optimize_mesh_3()` is a mesh optimization process
 * based on the minimization of a global energy function.
 *
 * In `lloyd_optimize_mesh_3()`, the minimized global energy may be interpreted
 * as the \f$ L^1\f$-norm of the error achieved
 * when the function \f$ x^2\f$ is interpolated on the mesh domain
 * using a piecewise linear function which is linear
 * in each cell of the Voronoi diagram of the mesh vertices.
 *
 * The optimizer `lloyd_optimize_mesh_3()` works in iterative steps.
 * At each iteration, mesh vertices are moved into
 * positions that bring to zero the energy gradient
 * and the Delaunay triangulation is updated.
 * Vertices on the mesh boundaries are handled
 * in a special way so as to preserve an accurate
 * representation of the domain boundaries.
 *
 * \tparam C3T3 a model of the concept `MeshComplex_3InTriangulation_3`.
 * \tparam MD a model of the concept `MeshDomain_3`.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param c3t3 the initial mesh that will be modified by the algorithm to represent the final optimized mesh.
 * @param domain the domain used to create the `c3t3` parameter
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{time_limit}
 *     \cgalParamDescription{to set up, in seconds, a CPU time limit after which the optimization process is stopped.
 *                           This time is measured using `CGAL::Real_timer`. 0 means that there is no time limit.}
 *     \cgalParamType{`double`}
 *     \cgalParamPrecondition{`time_limit >= 0`}
 *     \cgalParamDefault{0}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{max_iteration_number}
 *     \cgalParamDescription{limit on the number of performed iterations. 0 means that there is
 *                           no limit on the number of performed iterations.}
 *     \cgalParamPrecondition{`max_iteration_number >=0`}
 *     \cgalParamType{`int`}
 *     \cgalParamDefault{0}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{freeze_bound}
 *     \cgalParamDescription{designed to reduce running time of each optimization iteration.
 *                           Any vertex that has a displacement less than a given fraction of the length
 *                           of its shortest incident edge, is frozen (i.e.\ is not relocated).
 *                           The parameter `freeze_bound` gives the threshold ratio.
 *                           If it is set to 0, freezing of vertices is disabled.}
 *     \cgalParamPrecondition{`0<= freeze_bound <=1`}
 *     \cgalParamType{`double`}
 *     \cgalParamDefault{0.01}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{convergence}
 *     \cgalParamDescription{threshold ratio of stopping criterion based on convergence: the optimization process is stopped
 *                           when at the last iteration the displacement of any vertex is less than
 *                           a given fraction of the length of the shortest edge incident to that vertex.}
 *     \cgalParamPrecondition{`0 <=convergence <= 1`}
 *     \cgalParamType{`double`}
 *     \cgalParamDefault{0.02}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{do_freeze}
 *     \cgalParamDescription{completes the `freeze_bound` parameter. If it is set to `true` (default value),
 *                           frozen vertices will not move anymore in next iterations. Otherwise, at each iteration, any vertex that
 *                           moves, unfreezes all its incident vertices.}
 *     \cgalParamType{`bool`}
 *     \cgalParamDefault{true}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return  a value of type `CGAL::Mesh_optimization_return_code` which is:
 * <UL>
 * <LI>`CGAL::TIME_LIMIT_REACHED` when the time limit is reached.
 * <LI>`CGAL::MAX_ITERATION_NUMBER_REACHED` when `lloyd_optimize_mesh_3()` stops because it has performed `max_iteration_number` iterations.
 * <LI>`CGAL::CONVERGENCE_REACHED` when `lloyd_optimize_mesh_3()` stops because the convergence criterion
 * is achieved.
 * <LI>`CGAL::ALL_VERTICES_FROZEN` when all vertices have been frozen, when the
 * `do_freeze` parameter is set to true.
 * <LI>`CGAL::CANT_IMPROVE_ANYMORE` when `lloyd_optimize_mesh_3()` stops because
 * most vertices have been frozen, and no better convergence can be reached.
 * </UL>
 *
 * \cgalHeading{Example}
 *
 *
 * \code{.cpp}
 * // Lloyd-smoothing until convergence reaches 0.01, freezing vertices which
 * // move less than 0.001*shortest_incident_edge_length
 * lloyd_optimize_mesh_3(c3t3,
 *                       domain,
 *                       parameters::convergence(0.01).
 *                       parameters::freeze_bound(0.001).
 *                       parameters::do_freeze(true));
 *
 * \endcode
 *
 * \sa `CGAL::Mesh_optimization_return_code`
 * \sa `CGAL::make_mesh_3()`
 * \sa `CGAL::refine_mesh_3()`
 * \sa `CGAL::exude_mesh_3()`
 * \sa `CGAL::perturb_mesh_3()`
 * \sa `CGAL::odt_optimize_mesh_3()`
 *
 * \note This function requires the \ref thirdpartyEigen library.
 */
template<typename C3T3, typename MeshDomain, typename CGAL_NP_TEMPLATE_PARAMETERS>
Mesh_optimization_return_code lloyd_optimize_mesh_3(C3T3& c3t3, const MeshDomain& domain,const CGAL_NP_CLASS& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;
    std::size_t max_iterations = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 0);
    const double convergence_ratio = choose_parameter(get_parameter(np, internal_np::convergence_ratio), parameters::default_values_for_mesh_3::lloyd_convergence_ratio);
    const double freeze_bound = choose_parameter(get_parameter(np, internal_np::vertex_freeze_bound), parameters::default_values_for_mesh_3::lloyd_freeze_ratio);
    const double time_limit = choose_parameter(get_parameter(np, internal_np::maximum_running_time), parameters::default_values_for_mesh_3::time_limit);
    bool do_freeze = choose_parameter(get_parameter(np,internal_np::freeze),true);
    return lloyd_optimize_mesh_3_impl(c3t3, domain, time_limit, max_iterations, convergence_ratio, freeze_bound, do_freeze);
}

#ifndef DOXYGEN_RUNNING
// Overload handling parameters passed with operator=
template<typename C3T3, typename MeshDomain,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
Mesh_optimization_return_code lloyd_optimize_mesh_3(C3T3& c3t3, const MeshDomain& domain, const CGAL_NP_CLASS_1&  np1, const CGAL_NP_CLASS_2&  np2, const NP& ... nps)
{
    return lloyd_optimize_mesh_3(c3t3,domain, internal_np::combine_named_parameters(np1, np2, nps...));
}

template <typename C3T3, typename MeshDomain>
Mesh_optimization_return_code
lloyd_optimize_mesh_3_impl(C3T3& c3t3,
                           const MeshDomain& domain,
                           const double time_limit,
                           std::size_t max_iteration_number,
                           const double convergence,
                           const double freeze_bound
                           , const bool do_freeze)
{
  CGAL_precondition(
    !Mesh_3::internal::has_non_protecting_weights(c3t3.triangulation(), domain));

  typedef typename C3T3::Triangulation  Tr;

  typedef Mesh_3::Mesh_sizing_field<Tr>               Sizing;
  typedef typename Mesh_3::Lloyd_move<C3T3,Sizing>    Move;

  typedef typename
    Mesh_3::Mesh_global_optimizer<C3T3,MeshDomain,Move> Lloyd_optimizer;

  // Create optimizer
  Lloyd_optimizer opt (c3t3,
                       domain,
                       freeze_bound,
                       do_freeze,
                       convergence);

  // Set max time
  opt.set_time_limit(time_limit);

  // 1000 iteration max to avoid infinite loops
  if ( 0 == max_iteration_number )
    max_iteration_number = 1000;

  // Launch optimization
  return opt(static_cast<int>(max_iteration_number));
}
#endif //DOXYGEN_RUNNING

}  // end namespace CGAL

#endif // CGAL_LLOYD_OPTIMIZE_MESH_3_H
