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

#ifndef CGAL_GENERATE_HEXAHEDRAL_MESH_USING_TWO_REFINEMENT_H
#define CGAL_GENERATE_HEXAHEDRAL_MESH_USING_TWO_REFINEMENT_H 1

#include <CGAL/license/LCC_processing.h>

#include <CGAL/Hexmeshing/Hexmeshing_for_linear_cell_complex.h>
#include <CGAL/Hexmeshing/LCC_items.h>
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
    using parameters::choose_parameter;
    using parameters::get_parameter;

    bool trim=parameters::choose_parameter
      (parameters::get_parameter(np, internal_np::use_trimming), true);
    bool smooth=parameters::choose_parameter
      (parameters::get_parameter(np, internal_np::use_smoothing), true);

    using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
    using VPM = typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type;
    VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(CGAL::vertex_point, tmesh));

    internal::Hexmeshing_for_linear_cell_complex<GT, TriangleMesh, VPM> hdata(tmesh, vpm, cube_cells_per_dim);
    hdata.two_refinement(nb_levels, trim, smooth);
    if constexpr (std::is_same_v<LCC, Default>)
    { return hdata.lcc; }
    else
    { return LCC(hdata.lcc); }
  }
}

#endif // CGAL_GENERATE_HEXAHEDRAL_MESH_USING_TWO_REFINEMENT_H //
// EOF //
