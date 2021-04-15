// Copyright (c) 2016 CNRS and LIRIS' Establishments (France).
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
#ifndef LINEAR_CELL_COMPLEX_FWD_H
#define LINEAR_CELL_COMPLEX_FWD_H 1

#include <CGAL/memory.h>
#include <CGAL/Combinatorial_map_fwd.h>
#include <CGAL/Generalized_map_fwd.h>

namespace CGAL {

template<unsigned int d_, unsigned int ambient_dim,
         class Traits_, class Items_, class Alloc_ >
class CMap_linear_cell_complex_storage_1;

template<unsigned int d_, unsigned int ambient_dim,
         class Traits_, class Items_, class Alloc_ >
class GMap_linear_cell_complex_storage_1;

template <unsigned int d>
struct LCC_default_kernel;

template <unsigned int d_,
          class Kernel=typename LCC_default_kernel<d_>::type >
struct Linear_cell_complex_traits;

#if defined(CGAL_CMAP_DART_DEPRECATED) && !defined(CGAL_NO_DEPRECATED_CODE)
  template <unsigned int d>
  struct CGAL_DEPRECATED Linear_cell_complex_min_items;
#else
  struct Linear_cell_complex_min_items;
#endif

template < unsigned int d_, unsigned int ambient_dim,
           class Traits_,
           class Items_,
           class Alloc_,
           template<unsigned int,class,class,class,class>
           class Map,
           class Refs_,
           class Storage_>
class Linear_cell_complex_base;

template < unsigned int d_, unsigned int ambient_dim = d_,
           class Traits_ = Linear_cell_complex_traits<ambient_dim>,
#if defined(CGAL_CMAP_DART_DEPRECATED) && !defined(CGAL_NO_DEPRECATED_CODE)
           class Items_ = Linear_cell_complex_min_items<d_>,
#else
           class Items_ = Linear_cell_complex_min_items,
#endif
           class Alloc_ = CGAL_ALLOCATOR(int),
           template<unsigned int,class,class,class,class>
           class CMap = Combinatorial_map_base,
           class Storage_ = CMap_linear_cell_complex_storage_1<d_, ambient_dim,
                                                               Traits_, Items_,
                                                               Alloc_> >
  class Linear_cell_complex_for_combinatorial_map;

  template < unsigned int d_, unsigned int ambient_dim = d_,
             class Traits_ = Linear_cell_complex_traits<ambient_dim>,
             class Items_ = Linear_cell_complex_min_items,
             class Alloc_ = CGAL_ALLOCATOR(int),
             template<unsigned int,class,class,class,class>
             class CMap = Generalized_map_base,
             class Storage_ = GMap_linear_cell_complex_storage_1<d_, ambient_dim,
                                                                 Traits_, Items_,
                                                                 Alloc_> >
    class Linear_cell_complex_for_generalized_map;

#if !defined(CGAL_NO_DEPRECATED_CODE)
  template < unsigned int d_, unsigned int ambient_dim = d_,
             class Traits_ = Linear_cell_complex_traits<ambient_dim>,
#if defined(CGAL_CMAP_DART_DEPRECATED)
             class Items_ = Linear_cell_complex_min_items<d_>,
#else
             class Items_ = Linear_cell_complex_min_items,
#endif
             class Alloc_ = CGAL_ALLOCATOR(int),
             template<unsigned int,class,class,class,class>
             class CMap = Combinatorial_map_base,
             class Storage_ = CMap_linear_cell_complex_storage_1<d_, ambient_dim,
                                                                 Traits_, Items_,
                                                                 Alloc_> >
  class Linear_cell_complex;
#endif

} // CGAL

#endif // LINEAR_CELL_COMPLEX_FWD_H
