// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#include <CGAL/Skin_surface_traits_3.h>

#ifndef CGAL_SKIN_SURFACE_FILTERED_TRAITS_3_H
#define CGAL_SKIN_SURFACE_FILTERED_TRAITS_3_H

#include <CGAL/license/Skin_surface_3.h>

#include <CGAL/basic.h>
#include <CGAL/Skin_surface_traits_3.h>
#include <CGAL/Filtered_predicate.h>

namespace CGAL {

template <class K>
class Skin_surface_filtered_traits_3
  : public Skin_surface_traits_base_3<K>
{
  // Exact traits is based on the exact kernel.
  typedef Skin_surface_traits_3<typename K::Exact_kernel>         Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Skin_surface_traits_3<typename K::Approximate_kernel>   Filtering_traits;

  typedef typename K::C2E C2E;
  typedef typename K::C2F C2F;

  typedef Skin_surface_traits_base_3<K>           Base;
public:
  typedef Filtered_predicate<typename Exact_traits::Side_of_mixed_cell_3,
                             typename Filtering_traits::Side_of_mixed_cell_3,
                             C2E,
                             C2F >  Side_of_mixed_cell_3;

  enum { Has_filtered_predicates=true };
  enum { Has_static_filters=false };

  Skin_surface_filtered_traits_3() {}
  Skin_surface_filtered_traits_3(typename Base::FT s) : Base(s) {}

  // Only make the predicates filtered, not the constructions:
  Side_of_mixed_cell_3
  side_of_mixed_cell_3_object() const
  {
    return Side_of_mixed_cell_3(Base::get_shrink());
  }
};

} // namespace CGAL
#endif // CGAL_SKIN_SURFACE_FILTERED_TRAITS_3_H
