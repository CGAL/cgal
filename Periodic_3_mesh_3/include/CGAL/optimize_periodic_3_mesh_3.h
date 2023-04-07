// Copyright (c) 2017 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     :  Stephane Tayeb,
//                  Mael Rouxel-Labbé
//
//******************************************************************************
// File Description : Free functions for P3M3 optimizers, which are simply
//                    identical to Mesh_3's optimizer free functions.
//******************************************************************************

#ifndef CGAL_OPTIMIZE_PERIODIC_3_MESH_3_H
#define CGAL_OPTIMIZE_PERIODIC_3_MESH_3_H

#include <CGAL/license/Periodic_3_mesh_3.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/optimize_mesh_3.h>

namespace CGAL {
// ---------------------------------- pertuber ---------------------------------
/*!
 * \ingroup PkgPeriodic3Mesh3Functions
 *
 * The function `perturb_periodic_3_mesh_3()` is a mesh optimizer that
 * improves the quality of a Delaunay mesh by changing the positions of some vertices of the mesh.
 *
 * This function directly calls `perturb_mesh_3()`, but is provided for convenience.
 * Further information can be found in the documentation of the function `perturb_mesh_3()`.
 */
template<typename C3T3, typename MeshDomain, typename CGAL_NP_TEMPLATE_PARAMETERS>
Mesh_optimization_return_code perturb_periodic_3_mesh_3(C3T3& c3t3, const MeshDomain& domain, const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  double time_limit = choose_parameter(get_parameter(np,internal_np::maximum_running_time),parameters::default_values_for_mesh_3::time_limit);
  auto sliver_bound = choose_parameter(get_parameter(np,internal_np::lower_sliver_bound), parameters::default_values_for_mesh_3::perturb_sliver_bound);
  auto sliver_criterion = choose_parameter(get_parameter(np, internal_np::sliver_criteria), parameters::default_values_for_mesh_3::default_sliver_criterion(c3t3,sliver_bound));
  auto perturbation_vector = choose_parameter(get_parameter(np,internal_np::perturb_vector), default_perturbation_vector(c3t3,domain,sliver_criterion));
  return perturb_mesh_3_impl(c3t3, domain, time_limit, sliver_criterion, perturbation_vector);
}

// Overload handling parameters passed with operator=
template<typename C3T3, typename MeshDomain,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
Mesh_optimization_return_code perturb_periodic_3_mesh_3(C3T3& c3t3, const MeshDomain& domain,
                                                        const CGAL_NP_CLASS_1&  np1,
                                                        const CGAL_NP_CLASS_2&  np2,
                                                        const NP& ... nps)
{
  return perturb_periodic_3_mesh_3(c3t3,domain, internal_np::combine_named_parameters(np1, np2, nps...));
}
// Overload handling parameters passed with operator=
template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
Mesh_optimization_return_code perturb_periodic_3_mesh_3(const CGAL_NP_CLASS_1&  np1,
                                                        const CGAL_NP_CLASS_2&  np2,
                                                        const NP& ... nps)
{
  return perturb_periodic_3_mesh_3(internal_np::combine_named_parameters(np1, np2, nps...));
}

// ---------------------------------- exuder -----------------------------------
/*!
 * \ingroup PkgPeriodic3Mesh3Functions
 *
 * The function `exude_periodic_3_mesh_3()` performs a sliver exudation process
 * on a periodic Delaunay mesh.
 *
 * The sliver exudation process consists in optimizing the weights of vertices
 * of the periodic weighted Delaunay triangulation in such a way that slivers disappear
 * and the quality of the mesh improves.
 *
 * \warning This optimizer modifies the weight of vertices of the periodic triangulation
 * and, if called, must be the last optimizer to be called. If the mesh is refined after
 * this optimization has been performed, all improvements will be lost.
 *
 * This function directly calls `exude_mesh_3()`, but is provided for convenience.
 * Further information can be found in the documentation of the function `exude_mesh_3()`.
 */
template<typename C3T3, typename CGAL_NP_TEMPLATE_PARAMETERS>
Mesh_optimization_return_code exude_periodic_3_mesh_3(C3T3& c3t3,const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  double time_limit=choose_parameter(get_parameter(np,internal_np::maximum_running_time),0);
  double sliver_bound= choose_parameter(get_parameter(np,internal_np::lower_sliver_bound),parameters::default_values_for_mesh_3::exude_sliver_bound);
  return exude_mesh_3_impl(c3t3,time_limit,sliver_bound);
}
// Overload handling parameters passed with operator=
template<typename C3T3,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
Mesh_optimization_return_code exude_periodic_3_mesh_3(C3T3& c3t3,
                                                      const CGAL_NP_CLASS_1&  np1,
                                                      const CGAL_NP_CLASS_2&  np2,
                                                      const NP& ... nps)
{
  return exude_periodic_3_mesh_3(c3t3,internal_np::combine_named_parameters(np1, np2, nps...));
}
// Overload handling parameters passed with operator=
template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
Mesh_optimization_return_code exude_periodic_3_mesh_3(const CGAL_NP_CLASS_1&  np1,
                                                      const CGAL_NP_CLASS_2&  np2,
                                                      const NP& ... nps)
{
  return exude_periodic_3_mesh_3(internal_np::combine_named_parameters(np1, np2, nps...));
}



// ------------------------------ odt optimizer --------------------------------
/*!
 * \ingroup PkgPeriodic3Mesh3Functions
 *
 * The function `odt_optimize_periodic_3_mesh_3()` is a periodic mesh optimization
 * process based on the minimization of a global energy function.
 *
 * This function directly calls `odt_optimize_mesh_3()`, but is provided for convenience.
 * Further information can be found in the documentation of the function `odt_optimize_mesh_3()`.
 */
template<typename C3T3,typename MeshDomain,typename CGAL_NP_TEMPLATE_PARAMETERS>
Mesh_optimization_return_code odt_optimize_periodic_3_mesh_3(C3T3& c3t3, const MeshDomain& domain, const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  double time_limit=choose_parameter(get_parameter(np,internal_np::maximum_running_time),0);
  std::size_t max_iteration_number=choose_parameter(get_parameter(np,internal_np::number_of_iterations),0);
  double convergence=choose_parameter(get_parameter(np,internal_np::convergence_ratio),0.02);
  double freeze_bound=choose_parameter(get_parameter(np,internal_np::vertex_freeze_bound),0.01);
  bool do_freeze=choose_parameter(get_parameter(np,internal_np::freeze),true);
  return odt_optimize_mesh_3_impl(c3t3, domain, time_limit, max_iteration_number, convergence, freeze_bound, do_freeze);
}

// Overload handling parameters passed with operator=
template<typename C3T3, typename MeshDomain,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
Mesh_optimization_return_code odt_optimize_periodic_3_mesh_3(C3T3& c3t3, const MeshDomain& domain,
                                                             const CGAL_NP_CLASS_1&  np1,
                                                             const CGAL_NP_CLASS_2&  np2,
                                                             const NP& ... nps)
{
  return odt_optimize_periodic_3_mesh_3(c3t3, domain, internal_np::combine_named_parameters(np1, np2, nps...));
}
// Overload handling parameters passed with operator=
template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
Mesh_optimization_return_code odt_optimize_periodic_3_mesh_3(const CGAL_NP_CLASS_1&  np1,
                                                             const CGAL_NP_CLASS_2&  np2,
                                                             const NP& ... nps)
{
  return odt_optimize_periodic_3_mesh_3(internal_np::combine_named_parameters(np1, np2, nps...));
}


// ------------------------------- lloyd optimizer -----------------------------
/*!
 * \ingroup PkgPeriodic3Mesh3Functions
 *
 * The function `lloyd_optimize_periodic_3_mesh_3()` is a periodic mesh optimization
 * process based on the minimization of a global energy function.
 *
 * This function directly calls `lloyd_optimize_mesh_3()`, but is provided for convenience.
 * Further information can be found in the documentation of the function `lloyd_optimize_mesh_3()`.
 *
 * \note This function requires the \ref thirdpartyEigen library.
*/
template<typename C3T3, typename MeshDomain, typename CGAL_NP_TEMPLATE_PARAMETERS>
Mesh_optimization_return_code lloyd_optimize_periodic_3_mesh_3(C3T3& c3t3, const MeshDomain& domain,const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  std::size_t max_iterations = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 0);
  const double convergence_ratio = choose_parameter(get_parameter(np, internal_np::convergence_ratio), 0.02);
  const double freeze_bound = choose_parameter(get_parameter(np, internal_np::vertex_freeze_bound), 0.001);
  const double time_limit = choose_parameter(get_parameter(np, internal_np::maximum_running_time), 0.);
  bool do_freeze = choose_parameter(get_parameter(np,internal_np::freeze),true);
  return lloyd_optimize_mesh_3_impl(c3t3, domain, time_limit, max_iterations, convergence_ratio, freeze_bound, do_freeze);
}

// Overload handling parameters passed with operator=
template<typename C3T3, typename MeshDomain,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
Mesh_optimization_return_code lloyd_optimize_periodic_3_mesh_3(C3T3& c3t3, const MeshDomain& domain,
                                                               const CGAL_NP_CLASS_1&  np1,
                                                               const CGAL_NP_CLASS_2&  np2,
                                                               const NP& ... nps)
{
  return lloyd_optimize_periodic_3_mesh_3(c3t3,domain, internal_np::combine_named_parameters(np1, np2, nps...));
}

} // namespace CGAL

#endif // CGAL_OPTIMIZE_PERIODIC_3_MESH_3_H
