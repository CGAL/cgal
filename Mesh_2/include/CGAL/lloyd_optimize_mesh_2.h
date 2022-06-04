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

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_2/Mesh_global_optimizer_2.h>
#include <CGAL/Mesh_2/Lloyd_move_2.h>
#include <CGAL/Mesh_2/Mesh_sizing_field.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/iterator.h>
#include <CGAL/boost/parameter.h>
#include <boost/parameter/preprocessor.hpp>

#include <CGAL/Named_function_parameters.h>

#include <fstream>

#ifndef DOXYGEN_RUNNING
// see <CGAL/config.h>
CGAL_PRAGMA_DIAG_PUSH
// see <CGAL/boost/parameter.h>
CGAL_IGNORE_BOOST_PARAMETER_NAME_WARNINGS
#endif

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
 *     \cgalParamType{a class model of `ConstRange` with iterator being a model of `InputIterator` with `CDT::Point_2` as value type.}
 *     \cgalParamDefault{No seed.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{seeds_are_in_domain}
 *     \cgalParamDescription{specified if seeds indicates bounded connected components inside or outside of the domain.}
 *     \cgalParamType{`bool`}
 *     \cgalParamDefault{false}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{number_of_iterations}
 *     \cgalParamDescription{limit on the number of performed iterations. 0 means that there is
 *                           no limit on the number of performed iterations.}
 *     \cgalParamExtra{\pre `number_of_iterations >=0`}
 *     \cgalParamType{`int`}
 *     \cgalParamDefault{0}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{maximum_running_time}
 *     \cgalParamDescription{to set up, in seconds, a CPU time limit after which the optimization process is stopped.
 *                           This time is measured using `CGAL::Timer`. 0 means that there is no time limit.}
 *     \cgalParamType{`double`}
 *     \cgalParamExtra{\pre `maximum_running_time` \f$ \geq\f$ 0}
 *     \cgalParamDefault{0}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{freeze_bound}
 *     \cgalParamDescription{designed to reduce running time of each optimization iteration.
 *                           Any vertex that has a displacement less than a given fraction of the length
 *                           of its shortest incident edge, is frozen (i.e.\ is not relocated).
 *                           The parameter `freeze_bound` gives the threshold ratio.
 *                           If it is set to 0, freezing of vertices is disabled.}
 *     \cgalParamExtra{\pre `0<= freeze_bound <=1}
 *     \cgalParamType{`double`}
 *     \cgalParamDefault{0.001}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{convergence_ratio}
 *     \cgalParamDescription{threshold ratio of stopping criterion based on convergence: the optimization process is stopped
 *                           when at the last iteration the displacement of any vertex is less than
 *                           a given fraction of the length of the shortest edge incident to that vertex.}
 *     \cgalParamExtra{\pre `0 <=convergence <= 1`}
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

  int max_iterations = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 0);
  const double convergence_ratio = choose_parameter(get_parameter(np, internal_np::convergence_ratio), 0.001);
  const double freeze_bound = 0.001; /* choose_parameter(get_parameter(np, internal_np::freeze_bound), 0.001); */
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

  return lloyd_optimize_mesh_2_impl(cdt,
                                    max_iterations,
                                    convergence_ratio,
                                    freeze_bound,
                                    time_limit,
                                    seeds.begin(),
                                    seeds.end(),
                                    mark);
}

#ifndef DOXYGEN_RUNNING
#ifndef CGAL_NO_DEPRECATED_CODE
#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4003) // not enough actual parameters for macro
#endif
  // TODO: check how to use CGAL_DEPRECATED here
  BOOST_PARAMETER_FUNCTION(
  (Mesh_optimization_return_code),
  lloyd_optimize_mesh_2,
  parameters::tag,
  (required (in_out(cdt),*))
  (optional
    (max_iteration_number_, *, 0 )
    (convergence_, *, 0.001 )
    (time_limit_, *, 0. )
    (freeze_bound_, *, 0.001 )
    (seeds_begin_, *, CGAL::Emptyset_iterator())//see comments below
    (seeds_end_, *, CGAL::Emptyset_iterator())//see comments below
    (mark_, *, false) //if "false", seeds indicate "outside" regions
  )
  )
  {
    return lloyd_optimize_mesh_2_impl(cdt,
                                      max_iteration_number_,
                                      convergence_,
                                      freeze_bound_,
                                      time_limit_,
                                      seeds_begin_,
                                      seeds_end_,
                                      mark_);
  }

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_NO_DEPRECATED_CODE

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
                             const int max_iterations,
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
                             const int max_iterations,
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
    int nb_iterations = (0 == max_iterations)
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

CGAL_PRAGMA_DIAG_POP

#include <CGAL/enable_warnings.h>
#else

namespace CGAL {

/*!
\ingroup PkgMesh2Functions

\deprecated This function is deprecated since \cgal 5.5, the overload using `NamedParameters` must be used instead.

The function `lloyd_optimize_mesh_2()` is a mesh optimization process
based on the minimization of a global energy function.

In `lloyd_optimize_mesh_2()`, the minimized global energy may be interpreted
as the \f$ L^1\f$-norm of the error achieved
when the function \f$ x^2\f$ is interpolated on the mesh domain
using a piecewise linear function which is linear
in each cell of the Voronoi diagram of the mesh vertices.

The optimizer `lloyd_optimize_mesh_2()` works in iterative steps.
At each iteration, mesh vertices are moved into
positions that bring to zero the energy gradient
and the Delaunay triangulation is updated.
Vertices on the mesh boundaries are not moved.

\tparam CDT is required to be or derive from `CGAL::Constrained_Delaunay_triangulation_2`,
with vertex base and face base of its underlying `TriangulationDataStructure_2`
being models of `DelaunayMeshVertexBase_2` and `DelaunayMeshFaceBase_2`, respectively.
The argument `cdt`, passed by reference, provides the initial mesh
and is modified by the algorithm to represent the final optimized mesh.

\tparam PointIterator must be an iterator with value type `Kernel::Point_2`

The function has several optional parameters which are named parameters
(we use the Boost.Parameter library).
Therefore, when calling the function, the parameters can be provided in any order
provided that the names of the parameters are used
(see example at the bottom of this page).

\cgalHeading{Named Parameters}

- <b>`parameters::time_limit`</b>
is used to set up, in seconds,
a CPU time limit after which the optimization process is stopped. This time is
measured using `CGAL::Timer`.
The default value is 0 and means that there is no time limit.
\pre `time_limit` \f$ \geq\f$ 0

- <b>`parameters::%max_iteration_number`</b> sets a limit on the
number of performed iterations. The default value of 0 means that there is
no limit on the number of performed iterations.
\pre `max_iteration_number`\f$ \geq\f$ 0

- <b>`parameters::%convergence`</b> is a stopping criterion based on convergence:
the optimization process is stopped, when at the last iteration,
the displacement of any vertex is less than a given fraction of the
length of the shortest edge incident to that vertex.
The parameter `convergence` gives the threshold ratio.
\pre 0 \f$ \leq\f$ `convergence` \f$ \leq\f$ 1

- <b>`parameters::freeze_bound`</b> is designed to reduce running time of each
optimization iteration. Any vertex that has a displacement less than a given
fraction of the length of its shortest incident edge, is frozen (i.e.\ is
not relocated). The parameter `freeze_bound` gives the threshold ratio.
The default value is 0.001. If it is set to 0, freezing of vertices is disabled.
\pre 0 \f$ \leq\f$ `freeze_bound` \f$ \leq\f$ 1

- <b>`parameters::seeds_begin`</b> and <b>`parameters::seeds_end`</b>
are begin and end input iterators to iterate on seed points.
The sequence [`parameters::seeds_begin`, `parameters::seeds_end`)
defines the domain in which the mesh was generated, and should be optimized.

- <b>`parameters::mark`</b>. If `mark` is set to true, the mesh domain
is the union of the bounded connected components including at least one seed.
If `mark` is false, the domain is the union of the bounded components including
no seed. Note that the unbounded component of the plane is never optimized.
The default value is false.




\return
The function `lloyd_optimize_mesh_2()` returns a value of type `CGAL::Mesh_optimization_return_code`
which is:
<UL>
<LI>`CGAL::TIME_LIMIT_REACHED` when the time limit is reached.
<LI>`CGAL::MAX_ITERATION_NUMBER_REACHED` when `lloyd_optimize_mesh_2()` stops because it has performed `max_iteration_number` iterations.
<LI>`CGAL::CONVERGENCE_REACHED` when `lloyd_optimize_mesh_2()` stops because the convergence criterion
is met.
<LI>`CGAL::ALL_VERTICES_FROZEN` when all vertices have been frozen, when the
`freeze_bound` parameter is set to a positive value.
<LI>`CGAL::CANT_IMPROVE_ANYMORE` when `lloyd_optimize_mesh_2()` stops because
most vertices have been frozen, and no better convergence can be reached.
</UL>

\cgalHeading{Example}


\code{.cpp}
// Lloyd-smoothing until convergence reaches 0.01, freezing vertices which
// move less than 0.001*shortest_incident_edge_length
lloyd_optimize_mesh_2(cdt,
                      parameters::convergence=0.01,
                      parameters::freeze_bound=0.001);

\endcode

\sa `CGAL::Mesh_optimization_return_code`
\sa `CGAL::refine_Delaunay_mesh_2()`

*/

template<typename CDT, typename PointIterator>
CGAL::Mesh_optimization_return_code
lloyd_optimize_mesh_2(CDT& cdt,
  double parameters::time_limit=0,
  std::size_t parameters::max_iteration_number=0,
  double parameters::convergence=0.001,
  double parameters::freeze_bound = 0.001,
  PointIterator parameters::seeds_begin = PointIterator(),
  PointIterator parameters::seeds_end = PointIterator(),
  bool parameters::mark = false);

} /* namespace CGAL */

#endif // DOXYGEN_RUNNING

} //end namespace CGAL

#endif
