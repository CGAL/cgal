// Copyright (c) 2017 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_LINEAR_CELL_COMPLEX_FOR_BGL_COMBINATORIAL_MAP_HELPER_H
#define CGAL_LINEAR_CELL_COMPLEX_FOR_BGL_COMBINATORIAL_MAP_HELPER_H 1

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_bgl_min_items.h>
#include <CGAL/Linear_cell_complex_traits.h>
#include <CGAL/CMap_linear_cell_complex_storages.h>

namespace CGAL {

  /** @file Linear_cell_complex_for_bgl_combinatorial_map_helper.h
   * Definition of a linear cell complex based on combinatorial map, to use with BGL.
   * This cmap has points associated to all vertices and faces enables. Moreover
   * these cells have id.
   */

  // Linear_cell_complex_for_bgl_combinatorial_map class.
  template < unsigned int d_, unsigned int ambient_dim = d_,
             class Traits_ = Linear_cell_complex_traits<ambient_dim>,
             class Alloc_ = CGAL_ALLOCATOR(int),
             template<unsigned int,class,class,class,class>
             class CMap = Combinatorial_map_base,
             class Storage_ = CMap_linear_cell_complex_storage_1
             <d_, ambient_dim, Traits_,
              CGAL::Linear_cell_complex_bgl_min_items, Alloc_> >  
  struct Linear_cell_complex_for_bgl_combinatorial_map_helper
  {
  public:
    /// Type of the Linear_cell_complex_for_combinatorial_map.
    typedef Linear_cell_complex_for_combinatorial_map
            <d_, ambient_dim, Traits_, CGAL::Linear_cell_complex_bgl_min_items,
             Alloc_, CMap, Storage_> type;
  }; 

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_FOR_BGL_COMBINATORIAL_MAP_HELPER_H //
// EOF //
