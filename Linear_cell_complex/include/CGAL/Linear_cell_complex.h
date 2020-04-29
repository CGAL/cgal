// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_LINEAR_CELL_COMPLEX_H
#define CGAL_LINEAR_CELL_COMPLEX_H 1

#include <CGAL/Linear_cell_complex_fwd.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_traits.h>
#include <CGAL/Linear_cell_complex_min_items.h>
#include <CGAL/Combinatorial_map.h>
#include <CGAL/CMap_linear_cell_complex_storages.h>

namespace CGAL {

  /** @file Linear_cell_complex.h
   * Definition of a linear cell complex, i.e. a combinatorial map with points
   * associated to all vertices. Deprecated class.
   */

#if !defined(CGAL_NO_DEPRECATED_CODE)
  template < unsigned int d_, unsigned int ambient_dim,
             class Traits_, class Items_, class Alloc_,
             template<unsigned int,class,class,class,class> class CMap,
             class Storage_ >
  class CGAL_DEPRECATED Linear_cell_complex:
    public Linear_cell_complex_for_combinatorial_map<d_, ambient_dim, Traits_, Items_,
                                                     Alloc_, CMap, Storage_>
  {};
#endif

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_H //
// EOF //
