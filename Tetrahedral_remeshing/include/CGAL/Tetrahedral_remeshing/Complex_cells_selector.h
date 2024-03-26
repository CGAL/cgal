// Copyright (c) 2024 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj, Jean-Marc Thiery, Tamy Boubekeur

#ifndef CGAL_INTERNAL_C3T3_CELLS_SELECTOR_H
#define CGAL_INTERNAL_C3T3_CELLS_SELECTOR_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/property_map.h>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
  template<typename Tr>
  struct Complex_cells_selector
  {
    using key_type = typename Tr::Cell_handle;
    using value_type = bool;
    using reference = bool;
    using category = boost::read_write_property_map_tag;

    friend value_type get(const Complex_cells_selector&, const key_type& c)
    {
      using SI = typename Tr::Cell::Subdomain_index;
      return c->subdomain_index() != SI();
    }
    friend void put(Complex_cells_selector&, const key_type&, const value_type)
    {} //nothing to do : subdomain indices are updated in remeshing};
  };


} // end namespace Tetrahedral_remeshing
} // end namespace CGAL

#endif // CGAL_INTERNAL_C3T3_CELLS_SELECTOR_H
