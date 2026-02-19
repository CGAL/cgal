// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Efi Fogel       <efif@post.tau.ac.il>

#ifndef CGAL_SURFACE_SWEEP_2_ALGORITHMS_H
#define CGAL_SURFACE_SWEEP_2_ALGORITHMS_H

#include <CGAL/license/Surface_sweep_2.h>

/*! File
 *
 * \file Definition of the surface-sweep intersect function.
 */
#include <CGAL/Surface_sweep_2.h>
#include <CGAL/Surface_sweep_2/Intersection_polylines_visitor.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! Subdivide a sequence of input segments according to their pairwise intersections.
 * Each segment is subdivided into sub-segments, referred to as polylines.
 */
template <typename InputIterator, typename Points_, typename Polylines_, typename Base_traits_>
void compute_intersection_polylines(InputIterator begin, InputIterator end, Points_& points, Polylines_& polylines,
                                    Base_traits_& base_traits) {
  using Traits = Arr_curve_data_traits_2<Base_traits_, std::size_t>;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using Points = Points_;
  using Polylines = Polylines_;
  using Visitor = Intersection_polylines_visitor<Traits, Points, Polylines>;
  using Surface_sweep = Surface_sweep_2<Visitor>;

  // TODO currently polylines is forced to have [] access

  std::vector<X_monotone_curve_2> curves;
  curves.reserve(std::distance(begin, end));
  std::size_t i=0;
  for(auto it=begin; it!=end; ++it)
    curves.emplace_back(*it, i++);

  Traits traits(base_traits);
  Visitor visitor(points, polylines);
  Surface_sweep surface_sweep(&traits, &visitor);
  visitor.sweep(curves.begin(), curves.end());
}

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
