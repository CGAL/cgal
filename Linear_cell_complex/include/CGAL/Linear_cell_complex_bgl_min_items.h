// Copyright (c) 2017 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_LINEAR_CELL_COMPLEX_BGL_MIN_ITEMS_H
#define CGAL_LINEAR_CELL_COMPLEX_BGL_MIN_ITEMS_H 1

#include <CGAL/Cell_attribute_with_point_and_id.h>
#include <CGAL/Cell_attribute_with_id.h>
#include <CGAL/tuple.h>

namespace CGAL {

  struct Linear_cell_complex_bgl_min_items
  {
    /// Dart_wrapper defines the type of darts used.
    template <class LCC>
    struct Dart_wrapper
    {
      typedef CGAL::Tag_true Darts_with_id;
      typedef CGAL::Cell_attribute_with_point_and_id<LCC> Vertex_attribute;
      typedef CGAL::Cell_attribute_with_id<LCC> Face_attribute;
      typedef std::tuple<Vertex_attribute, void, Face_attribute> Attributes;
    };
  };

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_BGL_MIN_ITEMS_H //
// EOF //
