// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * file   data/3d/KernelFactory.h
 * author Gernot Walzl
 * date   2011-11-26
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_KERNEL_FACTORY_H
#define CGAL_STRAIGHT_SKELETON_3_KERNEL_FACTORY_H

#include <CGAL/Straight_skeleton_3/internal/debug.h>

#include <memory>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace kernel {

template <typename K>
class KernelFactory
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

public:
  static Point3SPtr createPoint3(const FT& x, const FT& y, const FT& z)
  {
    Point3SPtr result = std::make_shared<Point_3>(x, y, z);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Point3SPtr createPoint3(const Point_3& point)
  {
    Point3SPtr result = std::make_shared<Point_3>(point);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Point3SPtr createPoint3(Vector3SPtr vector)
  {
    CGAL_SS3_DEBUG_SPTR(vector);
    Point3SPtr result = std::make_shared<Point_3>((*vector)[0], (*vector)[1], (*vector)[2]);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Segment3SPtr createSegment3(Point3SPtr src, Point3SPtr dst)
  {
    CGAL_SS3_DEBUG_SPTR(src);
    CGAL_SS3_DEBUG_SPTR(dst);
    Segment3SPtr result = std::make_shared<Segment_3>(*src, *dst);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Segment3SPtr createSegment3(const Segment_3& seg)
  {
    Segment3SPtr result = std::make_shared<Segment_3>(seg);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Line3SPtr createLine3(Point3SPtr p, Point3SPtr q)
  {
    CGAL_SS3_DEBUG_SPTR(p);
    CGAL_SS3_DEBUG_SPTR(q);
    Line3SPtr result = std::make_shared<Line_3>(*p, *q);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Line3SPtr createLine3(const Line_3& line)
  {
    Line3SPtr result = std::make_shared<Line_3>(line);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Line3SPtr createLine3(Point3SPtr p, Vector3SPtr direction)
  {
    CGAL_SS3_DEBUG_SPTR(p);
    CGAL_SS3_DEBUG_SPTR(direction);
    Line3SPtr result = std::make_shared<Line_3>(*p, *direction);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Plane3SPtr createPlane3(const FT& a, const FT& b, const FT& c, const FT& d)
  {
    Plane3SPtr result = std::make_shared<Plane_3>(a, b, c, d);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Plane3SPtr createPlane3(Point3SPtr p, Point3SPtr q, Point3SPtr r)
  {
    CGAL_SS3_DEBUG_SPTR(p);
    CGAL_SS3_DEBUG_SPTR(q);
    CGAL_SS3_DEBUG_SPTR(r);
    CGAL_assertion(p != q && q != r && r != p);
    Plane3SPtr result = std::make_shared<Plane_3>(*p, *q, *r);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Plane3SPtr createPlane3(Point3SPtr p, Vector3SPtr normal)
  {
    Plane3SPtr result = std::make_shared<Plane_3>(*p, *normal);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Plane3SPtr createPlane3(const Plane_3& plane)
  {
    Plane3SPtr result = std::make_shared<Plane_3>(plane);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Vector3SPtr createVector3(const FT& x, const FT& y, const FT& z)
  {
    Vector3SPtr result = std::make_shared<Vector_3>(x, y, z);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Vector3SPtr createVector3(const Vector_3& vector)
  {
    Vector3SPtr result = std::make_shared<Vector_3>(vector);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Vector3SPtr createVector3(Point3SPtr point)
  {
    Vector3SPtr result = std::make_shared<Vector_3>(point->x(), point->y(), point->z());
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  /**
   * creates a vector with the same direction as the given line.
   */
  static Vector3SPtr createVector3(Line3SPtr line)
  {
    Vector3SPtr result = createVector3(line->to_vector());
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  /**
   * returns the normal vector of the plane.
   */
  static Vector3SPtr createVector3(Plane3SPtr plane)
  {
    Vector3SPtr result = createVector3(plane->orthogonal_vector());
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }
};

} // namespace kernel
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_KERNEL_FACTORY_H */
