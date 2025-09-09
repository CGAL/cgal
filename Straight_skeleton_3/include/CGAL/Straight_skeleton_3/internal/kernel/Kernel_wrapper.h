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
 * @file   algo/3d/KernelWrapper.h
 * @author Gernot Walzl
 * @date   2012-03-08
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_KERNEL_KERNEL_WRAPPER_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_KERNEL_KERNEL_WRAPPER_H

#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_factory.h>

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
class KernelWrapper
{
  using FT = typename K::FT;
  using Point_3 = typename K::Point_3;
  using Segment_3 = typename K::Segment_3;
  using Vector_3 = typename K::Vector_3;
  using Line_3 = typename K::Line_3;
  using Plane_3 = typename K::Plane_3;

  using Point3SPtr = std::shared_ptr<Point_3>;
  using Segment3SPtr = std::shared_ptr<Segment_3>;
  using Vector3SPtr = std::shared_ptr<Vector_3>;
  using Line3SPtr = std::shared_ptr<Line_3>;
  using Plane3SPtr = std::shared_ptr<Plane_3>;

  using KernelFactory = kernel::KernelFactory<K>;

public:
  static Point3SPtr intersection(Plane3SPtr plane1, Plane3SPtr plane2, Plane3SPtr plane3)
  {
    Point3SPtr result = Point3SPtr();

    auto res = CGAL::intersection(*plane1, *plane2, *plane3);
    if (!res) {
      CGAL_SS3_TRAITS_TRACE("Intersection of 3 planes is... not?");
    } else if (const Point_3* ipoint = std::get_if<Point_3>(&*res)) {
      result = KernelFactory::createPoint3(*ipoint);
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

    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Line3SPtr intersection(Plane3SPtr plane1, Plane3SPtr plane2)
  {
    Line3SPtr result = Line3SPtr();

    auto res = CGAL::intersection(*plane1, *plane2);
    if (!res) {
      CGAL_SS3_TRAITS_TRACE("Intersection of 2 planes is... not?");
    } else if (const Line_3 *iline = std::get_if<Line_3>(&*res)) {
      result = KernelFactory::createLine3(*iline);
    } else {
      CGAL_SS3_TRAITS_TRACE_CODE(if (const Plane_3* iplane = std::get_if<Plane_3>(&*res)) {)
      CGAL_SS3_TRAITS_TRACE_CODE(CGAL_USE(iplane);)
      CGAL_SS3_TRAITS_TRACE("Intersection of 2 planes is a plane");
      CGAL_SS3_TRAITS_TRACE_CODE(} else {)
      CGAL_SS3_TRAITS_TRACE("Intersection of plane and line is... something else?");
      CGAL_SS3_TRAITS_TRACE_CODE(})
    }

    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Point3SPtr intersection(Plane3SPtr plane, Line3SPtr line)
  {
    Point3SPtr result = Point3SPtr();

    auto res = CGAL::intersection(*plane, *line);
    if (!res) {
      CGAL_SS3_TRAITS_TRACE("Intersection of plane and line is... not?");
    } else if (const Point_3 *ipoint = std::get_if<Point_3>(&*res)) {
      result = KernelFactory::createPoint3(*ipoint);
    } else {
      CGAL_SS3_TRAITS_TRACE_CODE(if (const Line_3 *iline = std::get_if<Line_3>(&*res)) {)
      CGAL_SS3_TRAITS_TRACE_CODE(CGAL_USE(iline);)
      CGAL_SS3_TRAITS_TRACE("Intersection of plane and line is the line itself");
      CGAL_SS3_TRAITS_TRACE_CODE(} else {)
      CGAL_SS3_TRAITS_TRACE("Intersection of plane and line is... something else?");
      CGAL_SS3_TRAITS_TRACE_CODE(})

      CGAL_warning_msg(false, "intersection of plane and line failed to produce a point");
    }

    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Plane3SPtr bisector(Plane3SPtr plane1, Plane3SPtr plane2)
  {
    Plane3SPtr result = Plane3SPtr();

    CGAL_SS3_TRAITS_TRACE("Warning: bisector call brings a SQRT");

    // @tmp Hardcore disable the SQRT
    std::exit(1);

    result = KernelFactory::createPlane3(CGAL::bisector(*plane1, *plane2));
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static FT squared_distance(Point3SPtr p1, Point3SPtr p2)
  {
    return CGAL::squared_distance(*p1, *p2);
  }

  static FT squared_distance(Segment3SPtr segment, Point3SPtr point)
  {
    return CGAL::squared_distance(*segment, *point);
  }

  static FT squared_distance(Line3SPtr line, Point3SPtr point)
  {
    return CGAL::squared_distance(*line, *point);
  }

  static FT squared_distance(Plane3SPtr plane, Point3SPtr point)
  {
    return CGAL::squared_distance(*plane, *point);
  }

  static Plane3SPtr opposite(Plane3SPtr plane) {
    Plane3SPtr result = KernelFactory::createPlane3(plane->opposite());
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Line3SPtr opposite(Line3SPtr line)
  {
    Line3SPtr result = KernelFactory::createLine3(line->opposite());
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static bool isNormalizedPlane(Plane3SPtr plane)
  {
    const FT& a = plane->a();
    const FT& b = plane->b();
    const FT& c = plane->c();

    // inaccuracies during normalization since the sqrt is (usually) not exact
    return (a*a + b*b + c*c - 1) <= 1e-5;
  }

  static int side(Plane3SPtr plane, Point3SPtr point)
  {
    int result = 0;
    CGAL::Oriented_side side = plane->oriented_side(*point);
    if (side == CGAL::ON_POSITIVE_SIDE) result = 1;
    if (side == CGAL::ON_NEGATIVE_SIDE) result = -1;
    return result;
  }

  static int orientation(Line3SPtr line1, Line3SPtr line2)
  {
    int result = 0;
    Vector_3 dir1 = line1->to_vector();
    Vector_3 dir2 = line2->to_vector();
    Point_3 p0 = line1->point();
    Point_3 p1 = p0 + dir1;
    Point_3 p2 = line2->point();
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

  double angle(Vector3SPtr v1, Vector3SPtr v2)
  {
    double result = 0.0;
    FT arg = 0.0;
    arg = ((*v1)*(*v2)) / CGAL::disallowed_sqrt(v1->squared_length() * v2->squared_length());
    // fixes issues with floating point precision
    if (arg <= -1.0) {
      result = CGAL_PI;
    } else if (arg >= 1.0) {
      result = 0.0;
    } else {
      result = acos(CGAL::to_double(arg));
    }
    return result;
  }

  double angle(Line3SPtr line1, Line3SPtr line2)
  {
    double result = 0.0;
    Vector3SPtr v1 = KernelFactory::createVector3(line1);
    Vector3SPtr v2 = KernelFactory::createVector3(line2);
    result = angle(v1, v2);
    return result;
  }

  double angle(Plane3SPtr plane, Line3SPtr line)
  {
    double result = 0.0;
    Vector3SPtr v_plane = KernelFactory::createVector3(plane);
    Vector3SPtr v_line = KernelFactory::createVector3(line);
    result = angle(v_plane, v_line);
    if (result > CGAL_PI/2.0) {
      result = result - CGAL_PI/2.0;
    } else {
      result = CGAL_PI/2.0 - result;
    }
    return result;
  }

  /**
    * Computes the angle between the normal vectors of given planes.
    */
  static double angle(Plane3SPtr plane1, Plane3SPtr plane2)
  {
    double result = 0.0;
    Vector3SPtr v1 = KernelFactory::createVector3(plane1);
    Vector3SPtr v2 = KernelFactory::createVector3(plane2);
    result = angle(v1, v2);
    return result;
  }

  static bool isInside(Point3SPtr p, Point3SPtr p_box_1, Point3SPtr p_box_2)
  {
    bool result = true;
    for (unsigned int i = 0; i < 3; ++i) {
      if ((*p_box_1)[i] < (*p_box_2)[i]) {
        if (!( (*p_box_1)[i] <= (*p)[i] && (*p)[i] <= (*p_box_2)[i] )) {
          result = false;
        }
      } else {
        if (!( (*p_box_2)[i] <= (*p)[i] && (*p)[i] <= (*p_box_1)[i] )) {
          result = false;
        }
      }
    }
    return result;
  }

  static Vector3SPtr cross(Vector3SPtr v1, Vector3SPtr v2)
  {
    Vector3SPtr result = Vector3SPtr();
    result = KernelFactory::createVector3(CGAL::cross_product(*v1, *v2));
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Point3SPtr projection(Line3SPtr line, Point3SPtr point)
  {
    Point3SPtr result = Point3SPtr();
    result = KernelFactory::createPoint3(line->projection(*point));
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Point3SPtr projection(Plane3SPtr plane, Point3SPtr point)
  {
    Point3SPtr result = Point3SPtr();
    result = KernelFactory::createPoint3(plane->projection(*point));
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }
};

} // namespace kernel
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_KERNEL_KERNEL_WRAPPER_H */
