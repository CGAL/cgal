// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : perturb_mesh_3 function definition.
//******************************************************************************

#ifndef CGAL_PERTURB_MESH_3_H
#define CGAL_PERTURB_MESH_3_H


#include <CGAL/Mesh_3/global_parameters.h>
#include <CGAL/Mesh_3/sliver_criteria.h>
#include <CGAL/Mesh_3/Sliver_perturber.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/Mesh_3/parameters_defaults.h>

namespace CGAL {

BOOST_PARAMETER_FUNCTION(
  (Mesh_optimization_return_code),
  perturb_mesh_3,
  parameters::tag,
  (required (in_out(c3t3),*) (domain,*) )
  (optional
    (time_limit_, *, 0 )
    (sliver_bound_, *, parameters::default_values::perturb_sliver_bound )
  )
)
{
  return perturb_mesh_3_impl(c3t3, domain, time_limit_, sliver_bound_);
}



template <typename C3T3, typename MeshDomain> 
Mesh_optimization_return_code
perturb_mesh_3_impl(C3T3& c3t3,
                    const MeshDomain& domain,
                    const double time_limit,
                    const double sliver_bound)
{
  typedef MeshDomain Md;
  typedef typename C3T3::Triangulation::Geom_traits Gt;
  typedef Mesh_3::Min_dihedral_angle_criterion<Gt> Sc;
  //typedef Mesh_3::Radius_radio_criterion<Gt> Sc;
  
  typedef Mesh_3::Sliver_perturber<C3T3,Md,Sc>            Perturber;
  typedef Mesh_3::Sq_radius_perturbation<C3T3,Md,Sc>      Sq_radius;
  typedef Mesh_3::Volume_perturbation<C3T3,Md,Sc>         Volume;
  typedef Mesh_3::Dihedral_angle_perturbation<C3T3,Md,Sc> Dihedral_angle;
  typedef Mesh_3::Li_random_perturbation<C3T3,Md,Sc>      Li_random;
  
  // Build perturber
  Perturber perturber(c3t3,domain);
  
  perturber.add_perturbation(new Sq_radius(40,0.05));
  perturber.add_perturbation(new Volume(40,0.05));
  perturber.add_perturbation(new Dihedral_angle(40,0.05));
  perturber.add_perturbation(new Li_random(100,0.15));

  // Set max time
  perturber.set_time_limit(time_limit);
  
  // Launch perturber
  if ( sliver_bound != 0 )
    return perturber(sliver_bound);
  else
    return perturber();
}
  
  
  
} //namespace CGAL


#endif // CGAL_PERTURB_MESH_3_H
