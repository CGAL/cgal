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

#include <CGAL/boost/parameter.h>
#include <CGAL/Mesh_3/sliver_criteria.h>
#include <CGAL/Mesh_3/Sliver_perturber.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/Mesh_3/parameters_defaults.h>
#include <CGAL/internal/Mesh_3/check_weights.h>
#include <CGAL/use.h>

#include <boost/parameter/preprocessor.hpp>

#include <vector>

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
  perturb_mesh_3,
  parameters::tag,
  (required (in_out(c3t3),*) (domain,*) )
  (optional
    (time_limit_, *, 0 )
    (sliver_bound_, *, parameters::default_values::perturb_sliver_bound )
    (sliver_criterion_, *,
       parameters::default_values::default_sliver_criterion(c3t3,sliver_bound_))
    (perturbation_vector_, *,
       default_perturbation_vector(c3t3,domain,sliver_criterion_))
  )
)
{
  CGAL_USE(sliver_bound_);
  return perturb_mesh_3_impl(c3t3, domain, time_limit_, sliver_criterion_,
                             perturbation_vector_);
}
CGAL_PRAGMA_DIAG_POP

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

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



} //namespace CGAL


#include <CGAL/enable_warnings.h>

#endif // CGAL_PERTURB_MESH_3_H
