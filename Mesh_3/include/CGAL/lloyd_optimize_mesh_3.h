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

#include <CGAL/disable_warnings.h>

#include <CGAL/boost/parameter.h>
#include <CGAL/Mesh_3/Mesh_global_optimizer.h>
#include <CGAL/Mesh_3/Lloyd_move.h>
#include <CGAL/Mesh_3/Mesh_sizing_field.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/Mesh_3/parameters_defaults.h>
#include <CGAL/internal/Mesh_3/check_weights.h>

#include <boost/parameter/preprocessor.hpp>

namespace CGAL {

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4003) // not enough actual parameters for macro
#endif

// see <CGAL/config.h>
CGAL_PRAGMA_DIAG_PUSH
// see <CGAL/boost/parameter.h>
CGAL_IGNORE_BOOST_PARAMETER_NAME_WARNINGS

BOOST_PARAMETER_FUNCTION(
  (Mesh_optimization_return_code),
  lloyd_optimize_mesh_3,
  parameters::tag,
  (required (in_out(c3t3),*) (domain,*) )
  (optional
    (time_limit_, *, 0 )
    (max_iteration_number_, *, 0 )
    (convergence_, *, parameters::default_values::lloyd_convergence_ratio )
    (freeze_bound_, *, parameters::default_values::lloyd_freeze_ratio )
    (do_freeze_, *, parameters::default_values::do_freeze ))
)
{
  return lloyd_optimize_mesh_3_impl(c3t3, domain,
                                    time_limit_, max_iteration_number_,
                                    convergence_, freeze_bound_
                                    , do_freeze_);
}
CGAL_PRAGMA_DIAG_POP

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif


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


}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_LLOYD_OPTIMIZE_MESH_3_H
