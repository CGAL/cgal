// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
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
// File Description : defines constants (default values) for parameters of
// Mesh_3 global functions
//******************************************************************************

#ifndef CGAL_MESH_3_PARAMETERS_DEFAULTS_H
#define CGAL_MESH_3_PARAMETERS_DEFAULTS_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_3/sliver_criteria.h>

namespace CGAL {
namespace parameters { namespace default_values {

// exude_mesh_3  
const double exude_sliver_bound = 0.;

// perturb_mesh_3
const double perturb_sliver_bound = 0.;
template<typename C3T3>
CGAL::Mesh_3::Min_dihedral_angle_criterion
  <typename C3T3::Triangulation> 
  default_sliver_criterion(const C3T3& c3t3, const double& bound)
{
  typedef typename C3T3::Triangulation Tr;
  return CGAL::Mesh_3::Min_dihedral_angle_criterion<Tr>(bound, c3t3.triangulation());
}

// global optimizers
const bool do_freeze = true;

// lloyd_optimize_mesh_3
const double lloyd_freeze_ratio = 0.01;
const double lloyd_convergence_ratio = 0.02;

// odt_optimize_mesh_3
const double odt_freeze_ratio = 0.01;
const double odt_convergence_ratio = 0.02;

} } // end namespace parameters::default_values
} // end namespace CGAL

#endif // CGAL_MESH_3_PARAMETERS_DEFAULTS_H
