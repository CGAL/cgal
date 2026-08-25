// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Soichiro Yamazaki <soichiro19998@gmail.com>
//
#ifndef CGAL_HEXMESHING_GENERATE_TWO_REFINEMENT_MESH_H
#define CGAL_HEXMESHING_GENERATE_TWO_REFINEMENT_MESH_H 1

#include <CGAL/hexmeshing/Hexmeshing_for_linear_cell_complex_sequential.h>
#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <string>

namespace CGAL
{
  template <typename LCC, typename TriangleMesh>
  LCC generate_hexahedral_mesh_using_two_refinement
  (const TriangleMesh& tmesh, int cube_cells_per_dim, int nb_levels, bool trim=true, bool smooth=true)
  {
    internal::Hexmeshing_for_linear_cell_complex<TriangleMesh> hdata(tmesh, cube_cells_per_dim);
    hdata.two_refinement(nb_levels, trim, smooth);
    LCC lcc(hdata.lcc);
    return lcc;
  }
}

#endif // CGAL_HEXMESHING_GENERATE_TWO_REFINEMENT_MESH_H //
// EOF //
