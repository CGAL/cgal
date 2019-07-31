// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Laurent Saboret

#ifndef CGAL_RANDOM_SIMPLIFY_POINT_SET_H
#define CGAL_RANDOM_SIMPLIFY_POINT_SET_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/algorithm.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Iterator_range.h>

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iterator>
#include <set>
#include <algorithm>
#include <cmath>

namespace CGAL {

/**
   \ingroup PkgPointSetProcessing3Algorithms
   Randomly deletes a user-specified fraction of the input points.

   This method modifies the order of input points so as to pack all remaining points first,
   and returns an iterator over the first point to remove (see erase-remove idiom).
   For this reason it should not be called on sorted containers.

   \tparam PointRange is a model of `Range`.

   \param points input point range.
   \param removed_percentage percentage of points to remove.

   \return iterator over the first point to remove.
*/
template <typename PointRange>
typename PointRange::iterator
random_simplify_point_set(
  PointRange& points,
  double removed_percentage)
{
  CGAL_point_set_processing_precondition(removed_percentage >= 0 && removed_percentage <= 100);

  // Random shuffle
  CGAL::cpp98::random_shuffle (points.begin(), points.end());

  // Computes first iterator to remove
  std::size_t nb_points = std::distance(points.begin(), points.end());
  std::size_t first_index_to_remove = (std::size_t)(double(nb_points) * ((100.0-removed_percentage)/100.0));
  typename PointRange::iterator first_point_to_remove = points.begin();
  std::advance(first_point_to_remove, first_index_to_remove);

  return first_point_to_remove;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_RANDOM_SIMPLIFY_POINT_SET_H
