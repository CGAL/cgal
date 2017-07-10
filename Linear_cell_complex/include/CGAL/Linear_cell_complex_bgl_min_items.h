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
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_LINEAR_CELL_COMPLEX_BGL_MIN_ITEMS_H
#define CGAL_LINEAR_CELL_COMPLEX_BGL_MIN_ITEMS_H 1

#include <CGAL/Cell_attribute_with_point_and_id.h>
#include <CGAL/Cell_attribute_with_id.h>

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
      typedef CGAL::cpp11::tuple<Vertex_attribute, void, Face_attribute> Attributes;
    };
  };
  
} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_BGL_MIN_ITEMS_H //
// EOF //
