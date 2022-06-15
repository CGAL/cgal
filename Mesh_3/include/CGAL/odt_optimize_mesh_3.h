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
// File Description : odt_optimize_mesh_3 function definition.
//******************************************************************************

#ifndef CGAL_ODT_OPTIMIZE_MESH_3_H
#define CGAL_ODT_OPTIMIZE_MESH_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/Mesh_3/Mesh_global_optimizer.h>
#include <CGAL/Mesh_3/Odt_move.h>
#include <CGAL/Mesh_3/Mesh_sizing_field.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/Mesh_3/parameters_defaults.h>
#include <CGAL/Mesh_3/internal/check_weights.h>


namespace CGAL {
/*!
@ingroup PkgMesh3Functions

The function `odt_optimize_mesh_3()` is a mesh optimization process
based on the minimization of a global energy function.

In `odt_optimize_mesh_3()`, the minimized global energy may be interpreted
as the \f$ L^1\f$-norm of the error achieved
when the function \f$ x^2\f$ is interpolated on the mesh domain
using a piecewise linear function which is linear in each mesh cell.

The optimizer `odt_optimize_mesh_3()` works in iterative steps.
At each iteration, mesh vertices are moved into
positions that bring to zero the energy gradient
and the Delaunay triangulation is updated.
Vertices on the mesh boundaries are handled
in a special way so as to preserve an accurate
representation of the domain boundaries.

@pre `time_limit` \f$ \geq\f$ 0 and 0 \f$ \leq\f$ `convergence` \f$ \leq\f$ 1 and 0 \f$ \leq\f$ `freeze_bound` \f$ \leq\f$ 1

@tparam C3T3 iis required to be a model of the concept
`MeshComplex_3InTriangulation_3`.
The argument `c3t3`, passed by
reference, provides the initial mesh
and is modified by the algorithm
to represent the final optimized mesh.

 @tparam MeshDomain is required to be a model of the concept
`MeshDomain_3`. The argument `domain` must be the `MD`
object used to create the `c3t3` parameter.

 @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

 @param c3t3 the initial mesh that will be modified by the algorithm to represent the final optimized mesh.
 @param domain ...
 @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:

 \cgalNamedParamsBegin
   \cgalParamNBegin{time_limit_new}
     \cgalParamDescription{is used to set up, in seconds,
                           a CPU time limit after which the optimization process is stopped. This time is
                           measured using `Real_timer`.
                           The default value is 0 and means that there is no time limit.}
     \cgalParamType{`double`}
     \cgalParamDefault{0}

   \cgalParamNBegin{max_iteration_number_new}
     \cgalParamDescription{sets a limit on the number of performed iterations.
                           The default value of 0 means that there is
                           no limit on the number of performed iterations.}
     \cgalParamType{`std::size_t`}
     \cgalParamDefault{0}

   \cgalParamNBegin{convergence_new}
     \cgalParamDescription{is a stopping criterion based on convergence:
                           the optimization process is stopped, when at the last iteration,
                           the displacement of any vertex is less than a given percentage of the length
                           the shortest edge incident to that vertex.
                           The parameter `convergence` gives the threshold ratio.}
     \cgalParamType{`double`}
     \cgalParamDefault{0.02}

   \cgalParamNBegin{freeze_bound_new}
     \cgalParamDescription{is designed to reduce running time of each optimization iteration. Any vertex
                           that has a displacement less than a given percentage of the length of its shortest incident edge, is frozen (i.e.\ is
                           not relocated). The parameter `freeze_bound` gives the threshold ratio.}
     \cgalParamType{`double`}
     \cgalParamDefault{0.01}

   \cgalParamNBegin{do_freeze_new}
     \cgalParamDescription{completes the `freeze_bound` parameter. If it is set to `true` (default value),
                           frozen vertices will not move anymore in next iterations. Otherwise, at each iteration, any vertex that
                           moves, unfreezes all its incident vertices.}
     \cgalParamType{`bool`}
     \cgalParamDefault{true}
 \cgalNamedParamsEnd
 \return
The function `odt_optimize_mesh_3()` returns a value of type `CGAL::Mesh_optimization_return_code`
which is:
<UL>
<LI>`CGAL::TIME_LIMIT_REACHED` when the time limit is reached.
<LI>`CGAL::MAX_ITERATION_NUMBER_REACHED` when `odt_optimize_mesh_3()` stops because it has performed `max_iteration_number` iterations.
<LI>`CGAL::CONVERGENCE_REACHED` when `odt_optimize_mesh_3()` stops because the convergence criterion
is achieved.
<LI>`CGAL::ALL_VERTICES_FROZEN` when all vertices have been frozen, when the
`do_freeze` parameter is set to true.
<LI>`CGAL::CANT_IMPROVE_ANYMORE` when `odt_optimize_mesh_3()` stops because
most vertices have been frozen, and no better convergence can be reached.
</UL>

\cgalHeading{Example}

\code{.cpp}
// 100 iterations of ODT-smoothing
odt_optimize_mesh_3(c3t3,
                    domain,
                    parameters::max_iteration_number = 100,
                    parameters::convergence = 0);
\endcode

\sa `CGAL::Mesh_optimization_return_code`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`
\sa `CGAL::exude_mesh_3()`
\sa `CGAL::perturb_mesh_3()`
\sa `CGAL::lloyd_optimize_mesh_3()`
 */
template<typename C3T3,typename MeshDomain,typename CGAL_NP_TEMPLATE_PARAMETERS>
Mesh_optimization_return_code odt_optimize_mesh_3(C3T3& c3t3, MeshDomain& domain, const CGAL_NP_CLASS& np = parameters::default_values())
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
#ifndef DOXYGEN_RUNNING
   #ifndef CGAL_NO_DEPRECATED_CODE
   template<typename C3T3, typename MeshDomain, typename... NP_Pack>
   Mesh_optimization_return_code odt_optimize_mesh_3(C3T3& c3t3, MeshDomain& domain, const NP_Pack& ...nps)
   {
       return odt_optimize_mesh_3(c3t3, domain, internal_np::combine_named_parameters(nps...));
   }
   #endif //CGAL_NO_DEPRECATED_CODE

template <typename C3T3, typename MeshDomain>
Mesh_optimization_return_code
odt_optimize_mesh_3_impl(C3T3& c3t3,
                         const MeshDomain& domain,
                         const double time_limit,
                         std::size_t max_iteration_number,
                         const double convergence,
                         const double freeze_ratio,
                         const bool do_freeze )
{
  CGAL_precondition(
    !Mesh_3::internal::has_non_protecting_weights(c3t3.triangulation(), domain));

  typedef typename C3T3::Triangulation  Tr;

  typedef Mesh_3::Mesh_sizing_field<Tr>             Sizing;
  typedef typename Mesh_3::Odt_move<C3T3,Sizing>    Move;

  typedef typename
    Mesh_3::Mesh_global_optimizer<C3T3,MeshDomain,Move> Odt_optimizer;

  // Create optimizer
  Odt_optimizer opt(c3t3,
                    domain,
                    freeze_ratio,
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

#else
namespace CGAL {

/*!
\ingroup PkgMesh3Functions
\deprecated This function is deprecated since \cgal 5.5, the overload using `NamedParameters` must be used instead.

The function `odt_optimize_mesh_3()` is a mesh optimization process
based on the minimization of a global energy function.

In `odt_optimize_mesh_3()`, the minimized global energy may be interpreted
as the \f$ L^1\f$-norm of the error achieved
when the function \f$ x^2\f$ is interpolated on the mesh domain
using a piecewise linear function which is linear in each mesh cell.

The optimizer `odt_optimize_mesh_3()` works in iterative steps.
At each iteration, mesh vertices are moved into
positions that bring to zero the energy gradient
and the Delaunay triangulation is updated.
Vertices on the mesh boundaries are handled
in a special way so as to preserve an accurate
representation of the domain boundaries.

\pre `time_limit` \f$ \geq\f$ 0 and 0 \f$ \leq\f$ `convergence` \f$ \leq\f$ 1 and 0 \f$ \leq\f$ `freeze_bound` \f$ \leq\f$ 1

\tparam C3T3 is required to be a model of the concept
`MeshComplex_3InTriangulation_3`.
The argument `c3t3`, passed by
reference, provides the initial mesh
and is modified by the algorithm
to represent the final optimized mesh.

\tparam MD is required to be a model of the concept
`MeshDomain_3`. The argument `domain` must be the `MD`
object used to create the `c3t3` parameter.

The function has four optional parameters which are named parameters
(we use the Boost.Parameter library).
Therefore, when calling the function, the parameters can be provided in any order
provided that the names of the parameters are used
(see example at the bottom of this page).

\cgalHeading{Named Parameters}

- <b>`parameters::time_limit`</b>
is used to set up, in seconds,
a CPU time limit after which the optimization process is stopped. This time is
measured using `Real_timer`.
The default value is 0 and means that there is no time limit.

- <b>`parameters::max_iteration_number`</b> sets a limit on the
number of performed iterations. The default value of 0 means that there is
no limit on the number of performed iterations.

- <b>`parameters::convergence`</b> is a stopping criterion based on convergence:
the optimization process is stopped, when at the last iteration,
the displacement of any vertex is less than a given percentage of the length
the shortest edge incident to that vertex.
The parameter `convergence` gives the threshold ratio.

- <b>`parameters::freeze_bound`</b> is designed to reduce running time of each optimization iteration. Any vertex
that has a displacement less than a given percentage of the length (the  of its shortest incident edge, is frozen (i.e.\ is
not relocated). The parameter `freeze_bound` gives the threshold ratio.

- <b>`parameters::do_freeze`</b> completes the `freeze_bound` parameter. If it is set to `true` (default value),
frozen vertices will not move anymore in next iterations. Otherwise, at each iteration, any vertex that
moves, unfreezes all its incident vertices.


\return
The function `odt_optimize_mesh_3()` returns a value of type `CGAL::Mesh_optimization_return_code`
which is:
<UL>
<LI>`CGAL::TIME_LIMIT_REACHED` when the time limit is reached.
<LI>`CGAL::MAX_ITERATION_NUMBER_REACHED` when `odt_optimize_mesh_3()` stops because it has performed `max_iteration_number` iterations.
<LI>`CGAL::CONVERGENCE_REACHED` when `odt_optimize_mesh_3()` stops because the convergence criterion
is achieved.
<LI>`CGAL::ALL_VERTICES_FROZEN` when all vertices have been frozen, when the
`do_freeze` parameter is set to true.
<LI>`CGAL::CANT_IMPROVE_ANYMORE` when `odt_optimize_mesh_3()` stops because
most vertices have been frozen, and no better convergence can be reached.
</UL>

\cgalHeading{Example}

\code{.cpp}
// 100 iterations of ODT-smoothing
odt_optimize_mesh_3(c3t3,
                    domain,
                    parameters::max_iteration_number = 100,
                    parameters::convergence = 0);
\endcode

\sa `CGAL::Mesh_optimization_return_code`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`
\sa `CGAL::exude_mesh_3()`
\sa `CGAL::perturb_mesh_3()`
\sa `CGAL::lloyd_optimize_mesh_3()`

*/

template<typename C3T3, typename MD>
Mesh_optimization_return_code
odt_optimize_mesh_3(C3T3& c3t3,
                    const MD& domain,
                    double parameters::time_limit=0,
                    std::size_t parameters::max_iteration_number=0,
                    double parameters::convergence=0.02,
                    double parameters::freeze_bound = 0.01,
                    bool parameters::do_freeze=true);

} /* namespace CGAL */

#endif //DOXYGEN_RUNNING

}  // end namespace CGAL

#endif // CGAL_ODT_OPTIMIZE_MESH_3_H
