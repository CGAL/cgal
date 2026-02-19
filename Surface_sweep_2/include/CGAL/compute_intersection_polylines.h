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
template <typename InputIterator, typename Points_, typename Polylines_, typename Traits_>
void compute_intersection_polylines(InputIterator begin, InputIterator end, Points_& points, Polylines_& polylines,
                                    Traits_& traits) {
  using Traits = Traits_;
  using Points = Points_;
  using Polylines = Polylines_;
  using Visitor = Intersection_polylines_visitor<Traits, Points, Polylines>;
  using Surface_sweep = Surface_sweep_2<Visitor>;

  Visitor visitor(points, polylines);
  Surface_sweep surface_sweep(&traits, &visitor);
  visitor.sweep(begin, end);
}

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
