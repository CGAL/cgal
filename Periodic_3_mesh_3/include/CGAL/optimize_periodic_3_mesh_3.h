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

#include <CGAL/optimize_mesh_3.h>

namespace CGAL {

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4003) // not enough actual parameters for macro
#endif

// see <CGAL/config.h>
CGAL_PRAGMA_DIAG_PUSH
// see <CGAL/boost/parameter.h>
CGAL_IGNORE_BOOST_PARAMETER_NAME_WARNINGS

// ---------------------------------- pertuber ---------------------------------

BOOST_PARAMETER_FUNCTION(
  (Mesh_optimization_return_code),
  perturb_periodic_3_mesh_3,
  parameters::tag,
  (required (in_out(c3t3),*) (domain,*))
  (optional
    (time_limit_, *, 0)
    (sliver_bound_, *, parameters::default_values::perturb_sliver_bound)
    (sliver_criterion_, *, parameters::default_values::default_sliver_criterion(c3t3, sliver_bound_))
    (perturbation_vector_, *, default_perturbation_vector(c3t3,domain,sliver_criterion_))
  )
)
{
  CGAL_USE(sliver_bound_);
  return perturb_mesh_3_impl(c3t3, domain, time_limit_, sliver_criterion_,
                             perturbation_vector_);
}

// ---------------------------------- exuder -----------------------------------

BOOST_PARAMETER_FUNCTION(
  (Mesh_optimization_return_code),
  exude_periodic_3_mesh_3,
  parameters::tag,
  (required (in_out(c3t3),*))
  (optional
    (time_limit_, *, 0)
    (sliver_bound_, *, parameters::default_values::exude_sliver_bound)
  )
)
{
  return exude_mesh_3_impl(c3t3, time_limit_, sliver_bound_);
}


// ------------------------------ odt optimizer --------------------------------

BOOST_PARAMETER_FUNCTION(
  (Mesh_optimization_return_code),
  odt_optimize_periodic_3_mesh_3,
  parameters::tag,
  (required (in_out(c3t3),*) (domain,*))
  (optional
    (time_limit_, *, 0)
    (max_iteration_number_, *, 0)
    (convergence_, *, parameters::default_values::odt_convergence_ratio)
    (freeze_bound_, *, parameters::default_values::odt_freeze_ratio)
    (do_freeze_, *, parameters::default_values::do_freeze)
   )
)
{
  return odt_optimize_mesh_3_impl(c3t3, domain,
                                  time_limit_, max_iteration_number_,
                                  convergence_, freeze_bound_, do_freeze_);
}


// ------------------------------- lloyd optimizer -----------------------------

BOOST_PARAMETER_FUNCTION(
  (Mesh_optimization_return_code),
  lloyd_optimize_periodic_3_mesh_3,
  parameters::tag,
  (required (in_out(c3t3),*) (domain,*))
  (optional
    (time_limit_, *, 0)
    (max_iteration_number_, *, 0)
    (convergence_, *, parameters::default_values::lloyd_convergence_ratio)
    (freeze_bound_, *, parameters::default_values::lloyd_freeze_ratio)
    (do_freeze_, *, parameters::default_values::do_freeze)
   )
)
{
  return lloyd_optimize_mesh_3_impl(c3t3, domain,
                                    time_limit_, max_iteration_number_,
                                    convergence_, freeze_bound_, do_freeze_);
}

CGAL_PRAGMA_DIAG_POP

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

} // namespace CGAL

#endif // CGAL_OPTIMIZE_PERIODIC_3_MESH_3_H
