// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Soichiro Yamazaki <soichiro19998@gmail.com>
//
#ifndef CGAL_HEXMESHING_GENERATE_TWO_REFINEMENT_MESH_H
#define CGAL_HEXMESHING_GENERATE_TWO_REFINEMENT_MESH_H 1

#include <CGAL/hexmeshing/Hexmeshing_for_linear_cell_complex_sequential.h>
#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <string>

namespace CGAL
{
  template<typename LCC=Default, typename TriangleMesh,
           typename NamedParameters=parameters::Default_named_parameters>
  auto generate_hexahedral_mesh_using_two_refinement
  (const TriangleMesh& tmesh, int cube_cells_per_dim, int nb_levels,
   const NamedParameters& np=parameters::default_values())
  {
    bool trim=parameters::choose_parameter
      (parameters::get_parameter(np, internal_np::use_triming), true);
    bool smooth=parameters::choose_parameter
      (parameters::get_parameter(np, internal_np::use_smoothing), true);

    internal::Hexmeshing_for_linear_cell_complex<TriangleMesh> hdata(tmesh, cube_cells_per_dim);
    hdata.two_refinement(nb_levels, trim, smooth);
    if constexpr (std::is_same_v<LCC, Default>)
    { return hdata.lcc; }
    else
    { return LCC(hdata.lcc); }
  }
}

#endif // CGAL_HEXMESHING_GENERATE_TWO_REFINEMENT_MESH_H //
// EOF //
