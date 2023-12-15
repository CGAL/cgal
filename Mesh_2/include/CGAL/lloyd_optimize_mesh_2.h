// Copyright (c) 2014-2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Jane Tournois
//

#ifndef CGAL_LLOYD_OPTIMIZE_MESH_2_H
#define CGAL_LLOYD_OPTIMIZE_MESH_2_H

#include <CGAL/license/Mesh_2.h>

#include <CGAL/Mesh_2/Mesh_global_optimizer_2.h>
#include <CGAL/Mesh_2/Lloyd_move_2.h>
#include <CGAL/Mesh_2/Mesh_sizing_field.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/iterator.h>

#include <CGAL/Named_function_parameters.h>

#include <fstream>

namespace CGAL
{

/*!
 * @ingroup PkgMesh2Functions
 *
 * The function `lloyd_optimize_mesh_2()` is a mesh optimization process
 * based on the minimization of a global energy function.
 *
 * In `lloyd_optimize_mesh_2()`, the minimized global energy may be interpreted
 * as the \f$ L^1\f$-norm of the error achieved
 * when the function \f$ x^2\f$ is interpolated on the mesh domain
 * using a piecewise linear function which is linear
 * in each cell of the Voronoi diagram of the mesh vertices.
 *
 * The optimizer `lloyd_optimize_mesh_2()` works in iterative steps.
 * At each iteration, mesh vertices are moved into
 * positions that bring to zero the energy gradient
 * and the Delaunay triangulation is updated.
 * Vertices on the mesh boundaries are not moved.
 *
 * @tparam CDT is required to be or derive from `CGAL::Constrained_Delaunay_triangulation_2`,
 * with vertex base and face base of its underlying `TriangulationDataStructure_2`
 * being models of `DelaunayMeshVertexBase_2` and `DelaunayMeshFaceBase_2`, respectively.
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param cdt the initial mesh that will be modified by the algorithm to represent the final optimized mesh.
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{seeds}
 *     \cgalParamDescription{2D points used to define the domain to mesh.
 *                           If `seeds_are_in_domain==true`, the mesh domain is the union of
 *                           the bounded connected components including at least
 *                           one seed. If `seeds_are_in_domain==false`, the domain is the
 *                           union of the bounded components including no seed.}
 *     \cgalParamType{a class model of `ConstRange` whose iterator is a model of `InputIterator` with `CDT::Point_2` as value type.}
 *     \cgalParamDefault{No seed.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{seeds_are_in_domain}
 *     \cgalParamDescription{specified if seeds indicate bounded connected components inside or outside of the domain.}
 *     \cgalParamType{`bool`}
 *     \cgalParamDefault{false}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{max_iteration_number}
 *     \cgalParamDescription{limit on the number of performed iterations. 0 means that there is
 *                           no limit on the number of performed iterations.}
 *     \cgalParamPrecondition{`max_iteration_number >=0`}
 *     \cgalParamType{`int`}
 *     \cgalParamDefault{0}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{time_limit}
 *     \cgalParamDescription{CPU time limit (in seconds) after which the optimization process is stopped.
 *                           This time is measured using `CGAL::Real_timer`. 0 means that there is no time limit.}
 *     \cgalParamType{`double`}
 *     \cgalParamPrecondition{`time_limit` \f$ \geq\f$ 0}
 *     \cgalParamDefault{0}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{freeze_bound}
 *     \cgalParamDescription{designed to reduce running time of each optimization iteration.
 *                           Any vertex that has a displacement less than a given fraction of the length
 *                           of its shortest incident edge, is frozen (i.e.\ is not relocated).
 *                           The parameter `freeze_bound` gives the threshold ratio.
 *                           If it is set to 0, freezing of vertices is disabled.}
 *     \cgalParamPrecondition{`0<= freeze_bound <=1`}
 *     \cgalParamType{`double`}
 *     \cgalParamDefault{0.001}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{convergence}
 *     \cgalParamDescription{threshold ratio of stopping criterion based on convergence: the optimization process is stopped
 *                           when at the last iteration the displacement of any vertex is less than
 *                           a given fraction of the length of the shortest edge incident to that vertex.}
 *     \cgalParamPrecondition{`0 <=convergence <= 1`}
 *     \cgalParamType{`double`}
 *     \cgalParamDefault{0.001}
 *   \cgalParamNEnd
 *
 * \cgalNamedParamsEnd
 *
 * @returns an enum value providing some information about the outcome of the algorithm.
 *
 * \sa `CGAL::Mesh_optimization_return_code`
 * \sa `CGAL::refine_Delaunay_mesh_2()`
 */
template<typename CDT, typename CGAL_NP_TEMPLATE_PARAMETERS>
Mesh_optimization_return_code
lloyd_optimize_mesh_2(CDT& cdt, const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::get_parameter_reference;
  using parameters::is_default_parameter;

  std::size_t max_iterations = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 0);
  const double convergence_ratio = choose_parameter(get_parameter(np, internal_np::convergence_ratio), 0.001);
  const double freeze_bound = choose_parameter(get_parameter(np, internal_np::vertex_freeze_bound), 0.001);
  const double time_limit = choose_parameter(get_parameter(np, internal_np::maximum_running_time), 0.);
  // Seeds
  // {
  typedef std::vector<typename CDT::Point_2> Default_seeds;
  typedef typename internal_np::Lookup_named_param_def<internal_np::seeds_t,
                                                       CGAL_NP_CLASS,
                                                       Default_seeds>::reference Seeds;
  Default_seeds ds;
  Seeds seeds = choose_parameter(get_parameter_reference(np, internal_np::seeds), ds);
  // }
  const bool mark =  choose_parameter(get_parameter(np, internal_np::seeds_are_in_domain), false);

  if (is_default_parameter<CGAL_NP_CLASS,internal_np::i_seed_begin_iterator_t>::value ||
      is_default_parameter<CGAL_NP_CLASS,internal_np::i_seed_end_iterator_t>::value)
  {
    return lloyd_optimize_mesh_2_impl(cdt,
                                      max_iterations,
                                      convergence_ratio,
                                      freeze_bound,
                                      time_limit,
                                      seeds.begin(),
                                      seeds.end(),
                                      mark);
  }
  else
  {
    return lloyd_optimize_mesh_2_impl(cdt,
                                      max_iterations,
                                      convergence_ratio,
                                      freeze_bound,
                                      time_limit,
                                      choose_parameter(get_parameter(np, internal_np::i_seed_begin_iterator), CGAL::Emptyset_iterator()),
                                      choose_parameter(get_parameter(np, internal_np::i_seed_end_iterator), CGAL::Emptyset_iterator()),
                                      mark);
  }
}

#ifndef DOXYGEN_RUNNING
  // Overload handling parameters passed with operator=
  template<typename CDT,
           typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
           typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
           typename ... NP>
  Mesh_optimization_return_code
  lloyd_optimize_mesh_2(CDT& cdt,
                        const CGAL_NP_CLASS_1&  np1,
                        const CGAL_NP_CLASS_2&  np2,
                        const NP& ... nps)
  {
    return lloyd_optimize_mesh_2(cdt, internal_np::combine_named_parameters(np1, np2, nps...));
  }

  /**
  * this partial specialization is a workaround
  * to avoid compilation errors when seeds_begin and seeds_end are
  * not initialized. Indeed, there is no way to have a
  * "default empty iterator" for these named parameters.
  * Emptyset_iterator implements OutputIterator,
  * but stands here for "any empty input iterator"
  * (and any other type could).
  */
  template<typename CDT>
  Mesh_optimization_return_code
  lloyd_optimize_mesh_2_impl(CDT& cdt,
                             const std::size_t max_iterations,
                             const double convergence_ratio,
                             const double freeze_bound,
                             const double time_limit,
                             CGAL::Emptyset_iterator,
                             CGAL::Emptyset_iterator,
                             const bool mark)
  {
    std::list<typename CDT::Point> seeds;
    return lloyd_optimize_mesh_2_impl(cdt, max_iterations, convergence_ratio,
      freeze_bound, time_limit, seeds.begin(), seeds.end(), mark);
  }

  template<typename CDT, typename InputIterator>
  Mesh_optimization_return_code
  lloyd_optimize_mesh_2_impl(CDT& cdt,
                             const std::size_t max_iterations,
                             const double convergence_ratio,
                             const double freeze_bound,
                             const double time_limit,
                             InputIterator seeds_begin,
                             InputIterator seeds_end,
                             const bool mark)
  {
    typedef Mesh_2::Mesh_sizing_field<CDT>           Sizing;
    typedef Mesh_2::Lloyd_move_2<CDT, Sizing>        Mv;
    typedef Mesh_2::Mesh_global_optimizer_2<CDT, Mv> Optimizer;

    Optimizer lloyd(cdt,
                    convergence_ratio,
                    freeze_bound);
    lloyd.set_time_limit(time_limit);
    lloyd.set_seeds(seeds_begin, seeds_end, mark);

#ifdef CGAL_MESH_2_OPTIMIZERS_DEBUG
    std::ofstream os("before_lloyd.angles.txt");
    lloyd.output_angles_histogram(os);
    os.close();
#endif

    // 1000 iteration max to avoid infinite loop
    std::size_t nb_iterations = (0 == max_iterations)
      ? 1000
      : max_iterations;

    //run optimization
    Mesh_optimization_return_code rc = lloyd(nb_iterations);

#ifdef CGAL_MESH_2_OPTIMIZERS_DEBUG
    std::ofstream os2("after_lloyd.angles.txt");
    lloyd.output_angles_histogram(os2);
    os2.close();
#endif

    return rc;
  }
#endif // DOXYGEN_RUNNING

} //end namespace CGAL

#endif
