// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * file   algo/3d/KernelWrapper.h
 * author Gernot Walzl
 * date   2012-03-08
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_KERNEL_KERNEL_WRAPPER_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_KERNEL_KERNEL_WRAPPER_H

#include <CGAL/license/Straight_skeleton_3.h>

#include <CGAL/Straight_skeleton_3/internal/debug.h>

#include <CGAL/enum.h>
#include <CGAL/determinant.h>
#include <CGAL/intersections.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_3.h>

#include <cmath>
#include <optional>
#include <utility>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace kernel {

template <typename K>
class Kernel_wrapper
{
  using FT = typename K::FT;
  using Point_3 = typename K::Point_3;
  using Vector_3 = typename K::Vector_3;
  using Line_3 = typename K::Line_3;
  using Plane_3 = typename K::Plane_3;

public:
  static std::optional<Point_3> intersection(const Plane_3& plane1,
                                             const Plane_3& plane2,
                                             const Plane_3& plane3)
  {
    auto res = CGAL::intersection(plane1, plane2, plane3);
    if (!res) {
      CGAL_SS3_TRAITS_TRACE("Intersection of 3 planes is... not?");
    } else if (const Point_3* ipoint = std::get_if<Point_3>(&*res)) {
      return *ipoint;
    } else {
      CGAL_SS3_TRAITS_TRACE_CODE(if (const Line_3* iline = std::get_if<Line_3>(&*res)) {)
      CGAL_SS3_TRAITS_TRACE_CODE(CGAL_USE(iline);)
      CGAL_SS3_TRAITS_TRACE("Intersection of 3 planes is a line");
      CGAL_SS3_TRAITS_TRACE_CODE(} else if(const Plane_3* iplane = std::get_if<Plane_3>(&*res)) {)
      CGAL_SS3_TRAITS_TRACE_CODE(CGAL_USE(iplane);)
      CGAL_SS3_TRAITS_TRACE("Intersection of 3 planes is a plane");
      CGAL_SS3_TRAITS_TRACE_CODE(} else {)
      CGAL_SS3_TRAITS_TRACE("Intersection of 3 planes is... something else?");
      CGAL_SS3_TRAITS_TRACE_CODE(})
    }
    return std::nullopt;
  }

  static std::optional<Line_3> intersection(const Plane_3& plane1,
                                            const Plane_3& plane2)
  {
    auto res = CGAL::intersection(plane1, plane2);
    if (!res) {
      CGAL_SS3_TRAITS_TRACE("Intersection of 2 planes is... not?");
    } else if (const Line_3 *iline = std::get_if<Line_3>(&*res)) {
      return *iline;
    } else {
      CGAL_SS3_TRAITS_TRACE_CODE(if (const Plane_3* iplane = std::get_if<Plane_3>(&*res)) {)
      CGAL_SS3_TRAITS_TRACE_CODE(CGAL_USE(iplane);)
      CGAL_SS3_TRAITS_TRACE("Intersection of 2 planes is a plane");
      CGAL_SS3_TRAITS_TRACE_CODE(} else {)
      CGAL_SS3_TRAITS_TRACE("Intersection of plane and line is... something else?");
      CGAL_SS3_TRAITS_TRACE_CODE(})
    }
    return std::nullopt;
  }

  static std::optional<Point_3> intersection(const Plane_3& plane,
                                             const Line_3& line)
  {
    auto res = CGAL::intersection(plane, line);
    if (!res) {
      CGAL_SS3_TRAITS_TRACE("Intersection of plane and line is... not?");
    } else if (const Point_3 *ipoint = std::get_if<Point_3>(&*res)) {
      return *ipoint;
    } else {
      CGAL_SS3_TRAITS_TRACE_CODE(if (const Line_3 *iline = std::get_if<Line_3>(&*res)) {)
      CGAL_SS3_TRAITS_TRACE_CODE(CGAL_USE(iline);)
      CGAL_SS3_TRAITS_TRACE("Intersection of plane and line is the line itself");
      CGAL_SS3_TRAITS_TRACE_CODE(} else {)
      CGAL_SS3_TRAITS_TRACE("Intersection of plane and line is... something else?");
      CGAL_SS3_TRAITS_TRACE_CODE(})

      CGAL_warning_msg(false, "intersection of plane and line failed to produce a point");
    }
    return std::nullopt;
  }

  static Plane_3 bisector(const Plane_3& plane1,
                          const Plane_3& plane2)
  {
    CGAL_SS3_TRAITS_TRACE("Warning: bisector() calls a square root!");
    return CGAL::bisector(plane1, plane2);
  }

  static bool has_normalized_plane(const Plane_3& plane)
  {
    CGAL_SS3_DEBUG_SPTR(plane);

    const FT& a = plane.a();
    const FT& b = plane.b();
    const FT& c = plane.c();

    // inaccuracies during normalization since the sqrt is (usually) not exact
    return (square(a) + square(b) + square(c) - FT(1)) <= 1e-5;
  }

  static int side(const Plane_3& plane, const Point_3& point)
  {
    int result = 0;
    CGAL::Oriented_side side = plane.oriented_side(point);
    if (side == CGAL::ON_POSITIVE_SIDE) result = 1;
    if (side == CGAL::ON_NEGATIVE_SIDE) result = -1;
    return result;
  }

  static int orientation(const Line_3& line1,
                         const Line_3& line2)
  {
    int result = 0;
    Vector_3 dir1 = line1.to_vector();
    Vector_3 dir2 = line2.to_vector();
    Point_3 p0 = line1.point();
    Point_3 p1 = p0 + dir1;
    Point_3 p2 = line2.point();
    Plane_3 plane(p0, p1, p2);
    Point_3 point = p2 + dir2;
    CGAL::Oriented_side side = plane.oriented_side(point);
    if (side == CGAL::ON_POSITIVE_SIDE) {
      result = 1;
    } else if (side == CGAL::ON_NEGATIVE_SIDE) {
      result = -1;
    }
    return result;
  }
};

} // namespace kernel
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_KERNEL_KERNEL_WRAPPER_H */
