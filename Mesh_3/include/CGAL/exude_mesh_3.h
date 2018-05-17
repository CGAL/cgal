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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : exude_mesh_3 function definition.
//******************************************************************************

#ifndef CGAL_EXUDE_MESH_3_H
#define CGAL_EXUDE_MESH_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_3/sliver_criteria.h>
#include <CGAL/Mesh_3/Slivers_exuder.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/Mesh_3/parameters_defaults.h>
#include <CGAL/boost/parameter.h>

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
  exude_mesh_3,
  parameters::tag,
  (required (in_out(c3t3),*) )
  (optional
    (time_limit_, *, 0 )
    (sliver_bound_, *, parameters::default_values::exude_sliver_bound )
  )
)
{
  return exude_mesh_3_impl(c3t3, time_limit_, sliver_bound_);
}
CGAL_PRAGMA_DIAG_POP

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif



template <typename C3T3>
Mesh_optimization_return_code
exude_mesh_3_impl(C3T3& c3t3,
                  const double time_limit,
                  const double sliver_bound)
{
  typedef typename C3T3::Triangulation Tr;
  typedef Mesh_3::Min_dihedral_angle_criterion<Tr> Sc;
  //typedef Mesh_3::Radius_radio_criterion<Tr> Sc;
  typedef typename Mesh_3::Slivers_exuder<C3T3, Sc> Exuder;
  
  // Create exuder
  Sc criterion(sliver_bound, c3t3.triangulation());
  Exuder exuder(c3t3, criterion);

  // Set time_limit
  exuder.set_time_limit(time_limit);
  
  // Launch exudation
  return exuder();
}


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_EXUDE_MESH_3_H
