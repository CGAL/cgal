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
#ifndef CGAL_LINEAR_CELL_COMPLEX_MIN_ITEMS_H
#define CGAL_LINEAR_CELL_COMPLEX_MIN_ITEMS_H 1

#if defined(CGAL_CMAP_DART_DEPRECATED) && !defined(CGAL_NO_DEPRECATED_CODE)
#include <CGAL/Dart.h>
#endif

#include <CGAL/Linear_cell_complex_fwd.h>
#include <CGAL/Cell_attribute_with_point.h>
#include <CGAL/tuple.h>

namespace CGAL {

  /** @file Linear_cell_complex_min_items.h
   * Definition of min item class for map with points.
   */

  /** Minimal items for linear cell complexes.
   * Linear_cell_complex_min_items defines what is the item class
   * for a linear cell complex. It provides definitions for attributes
   * associated to vertices (containing points), and information associated with darts.
   */
#if defined(CGAL_CMAP_DART_DEPRECATED) && !defined(CGAL_NO_DEPRECATED_CODE)
  template <unsigned int d>
  struct CGAL_DEPRECATED Linear_cell_complex_min_items
  {
    /// Dart_wrapper defines the type of darts used.
    template <class LCC>
    struct Dart_wrapper
    {
      typedef CGAL::Dart<d, LCC> Dart;

      typedef CGAL::Cell_attribute_with_point<LCC> Vertex_attrib;
      typedef std::tuple<Vertex_attrib>    Attributes;
    };
  };
#else
  struct Linear_cell_complex_min_items
  {
    /// Dart_wrapper defines the type of darts used.
    template <class LCC>
    struct Dart_wrapper
    {
      typedef CGAL::Cell_attribute_with_point<LCC> Vertex_attrib;
      typedef std::tuple<Vertex_attrib>    Attributes;
    };
  };
#endif

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_MIN_ITEMS_H //
// EOF //
