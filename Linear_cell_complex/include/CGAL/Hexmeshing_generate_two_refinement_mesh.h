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
#include <CGAL/hexmeshing/Hexmeshing_for_linear_cell_complex_sequential.h>
#include <CGAL/hexmeshing/Hexmeshing_outer_alias.h>
#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <string>

namespace CGAL {
  internal::Hexmeshing::LCC generate_two_refinement_mesh(const std::string& file, int cube_cells_per_dim, int nb_levels, bool trim=false) {
    internal::Hexmeshing_for_linear_cell_complex hdata(file, cube_cells_per_dim);
    hdata.two_refinement(nb_levels, trim);
    return hdata.lcc;
  }

  template <typename TriangleMesh=internal::Hexmeshing::Polyhedron>
  internal::Hexmeshing::LCC generate_two_refinement_mesh(TriangleMesh& poly, int cube_cells_per_dim, int nb_levels, bool trim=false) {
    internal::Hexmeshing_for_linear_cell_complex<TriangleMesh> hdata(poly, cube_cells_per_dim);
    hdata.two_refinement(nb_levels, trim);
    return hdata.lcc;
  }
}