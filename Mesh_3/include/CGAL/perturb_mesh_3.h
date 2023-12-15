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
// File Description : perturb_mesh_3 function definition.
//******************************************************************************

#ifndef CGAL_PERTURB_MESH_3_H
#define CGAL_PERTURB_MESH_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_3/sliver_criteria.h>
#include <CGAL/Mesh_3/Sliver_perturber.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/Mesh_3/parameters_defaults.h>
#include <CGAL/Mesh_3/internal/check_weights.h>
#include <CGAL/use.h>
#include <CGAL/Named_function_parameters.h>

#include <vector>

namespace CGAL {
/*!
 * \ingroup PkgMesh3Functions
 *
 * The function `perturb_mesh_3()` is a mesh optimizer that
 * improves the quality of a Delaunay mesh
 * by changing the positions of some vertices of the mesh.
 *
 * The perturber tries to improve the dihedral angles of the worst cells in the mesh
 * degree by degree: the
 * step number `n` is considered as successful
 * if after this step the worst tetrahedron of the mesh has a minimal dihedral
 * angle larger than `n` degrees.
 * The perturber exits if this is not the case.
 *
 * \tparam C3T3 a model of the concept `MeshComplex_3InTriangulation_3`
 * \tparam MD a model of the concept `MeshDomain_3`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 *  @param c3t3 the initial mesh and is modified by the algorithm to represent the final optimized mesh
 *  @param domain the domain used to create the `c3t3` parameter
 *  @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{time_limit}
 *     \cgalParamDescription{is used to set up, in seconds, a CPU time limit after which the optimization process
 *                           is stopped. This time is measured using the `Real_timer` class. The default value is
 *                           0 and means that there is no time limit.}
 *     \cgalParamType{`double`}
 *     \cgalParamPrecondition{`0 <= sliver_bound <= 180`}
 *   \cgalParamDefault{0}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{sliver_bound}
 *     \cgalParamDescription{is designed to give, in degrees, a targeted lower bound on dihedral angles of mesh cells.
 *                          The exudation process considers in turn all the mesh cells that have a smallest dihedral
 *                          angle less than sliver_bound and tries to make them disappear by weighting their vertices.
 *                          The optimization process stops when every cell in the mesh achieves this quality. The
 *                          default value is 0 and means that there is no targeted bound: the exuder then runs as long
 *                          as it can improve the smallest dihedral angles of the set of cells incident to some vertices.}
 *     \cgalParamType{`double`}
 *     \cgalParamPrecondition{`time_limit >= 0`}
 *     \cgalParamDefault{0}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return a value of type `CGAL::Mesh_optimization_return_code` which is:
 * <UL>
 *   <LI>`CGAL::BOUND_REACHED` when the targeted bound for the smallest dihedral angle in the mesh is reached.
 *   <LI>`CGAL::TIME_LIMIT_REACHED` when the time limit is reached.
 *   <LI>`CGAL::CANT_IMPROVE_ANYMORE` when the perturbation process stops because the last step is unsuccessful.
 * </UL>
 *
 * \cgalHeading{Example}
 *
 * \code{.cpp}
 * // Perturb until every dihedral angle of the mesh is >= 10 degrees
 * // No time bound is set
 * perturb_mesh_3(c3t3,
 *                domain,
 *                parameters::sliver_bound = 10);
 * \endcode
 *
 * \sa `CGAL::Mesh_optimization_return_code`
 * \sa `CGAL::make_mesh_3()`
 * \sa `CGAL::refine_mesh_3()`
 * \sa `CGAL::exude_mesh_3()`
 * \sa `CGAL::lloyd_optimize_mesh_3()`
 * \sa `CGAL::odt_optimize_mesh_3()`
 *
 */
template<typename C3T3, typename MeshDomain, typename CGAL_NP_TEMPLATE_PARAMETERS>
Mesh_optimization_return_code perturb_mesh_3(C3T3& c3t3, const MeshDomain& domain, const CGAL_NP_CLASS& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;
    double time_limit = choose_parameter(get_parameter(np,internal_np::maximum_running_time),parameters::default_values_for_mesh_3::time_limit);
    auto sliver_bound = choose_parameter(get_parameter(np,internal_np::lower_sliver_bound), parameters::default_values_for_mesh_3::perturb_sliver_bound);
    auto sliver_criterion = choose_parameter(get_parameter(np, internal_np::sliver_criteria), parameters::default_values_for_mesh_3::default_sliver_criterion(c3t3,sliver_bound));
    auto perturbation_vector = choose_parameter(get_parameter(np,internal_np::perturb_vector), default_perturbation_vector(c3t3,domain,sliver_criterion));
    return perturb_mesh_3_impl(c3t3, domain, time_limit, sliver_criterion, perturbation_vector);
}


#ifndef DOXYGEN_RUNNING
// Overload handling parameters passed with operator=
template<typename C3T3, typename MeshDomain,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
Mesh_optimization_return_code perturb_mesh_3(C3T3& c3t3, const MeshDomain& domain,
                                             const CGAL_NP_CLASS_1&  np1,
                                             const CGAL_NP_CLASS_2&  np2,
                                             const NP& ... nps)
{
  return perturb_mesh_3(c3t3,domain, internal_np::combine_named_parameters(np1, np2, nps...));
}

template <typename C3T3,
          typename MeshDomain,
          typename SliverCriterion>
std::vector<typename Mesh_3::Sliver_perturber<C3T3,MeshDomain,SliverCriterion>::Perturbation*>
default_perturbation_vector(const C3T3&,
                            const MeshDomain&,
                            const SliverCriterion&)
{
  typedef MeshDomain Md;
  typedef SliverCriterion Sc;
  typedef Mesh_3::Sliver_perturber<C3T3,Md,Sc>            Perturber;
  typedef typename Perturber::Perturbation                Perturbation;

  typedef Mesh_3::Sq_radius_perturbation<C3T3,Md,Sc>      Sq_radius;
  typedef Mesh_3::Volume_perturbation<C3T3,Md,Sc>         Volume;
  typedef Mesh_3::Dihedral_angle_perturbation<C3T3,Md,Sc> Dihedral_angle;
  typedef Mesh_3::Li_random_perturbation<C3T3,Md,Sc>      Li_random;

  std::vector<Perturbation*> perturbation_vect;
  perturbation_vect.push_back(new Sq_radius(40,0.05));
  perturbation_vect.push_back(new Volume(40,0.05));
  perturbation_vect.push_back(new Dihedral_angle(40,0.05));
  perturbation_vect.push_back(new Li_random(100,0.15));

  return perturbation_vect;
}


template <typename C3T3,
          typename MeshDomain,
          typename SliverCriterion,
          typename PPerturbationVector>
Mesh_optimization_return_code
perturb_mesh_3_impl(C3T3& c3t3,
                    const MeshDomain& domain,
                    const double time_limit,
                    const SliverCriterion& sliver_criterion,
                    const PPerturbationVector& perturbation_vector)
{
  CGAL_precondition(
    !Mesh_3::internal::has_non_protecting_weights(c3t3.triangulation(), domain));

  typedef MeshDomain Md;
  typedef SliverCriterion Sc;

  typedef Mesh_3::Sliver_perturber<C3T3,Md,Sc> Perturber;

  // Build perturber
  Perturber perturber(c3t3, domain, sliver_criterion);

  // Add perturbations
  for(std::size_t i = 0; i < perturbation_vector.size(); ++i)
    perturber.add_perturbation( perturbation_vector[i] );

  // Set max time
  perturber.set_time_limit(time_limit);

  // Launch perturber
  return perturber();
}
#endif //DOXYGEN_RUNNING

} //namespace CGAL


#include <CGAL/enable_warnings.h>

#endif // CGAL_PERTURB_MESH_3_H
