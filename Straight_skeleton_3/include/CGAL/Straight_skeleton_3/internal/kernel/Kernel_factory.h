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
 * @file   data/3d/KernelFactory.h
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_KERNEL_FACTORY_H
#define CGAL_STRAIGHT_SKELETON_3_KERNEL_FACTORY_H

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
  using Sphere_3 = typename K::Sphere_3;

  using Point3SPtr = std::shared_ptr<Point_3>;
  using Segment3SPtr = std::shared_ptr<Segment_3>;
  using Vector3SPtr = std::shared_ptr<Vector_3>;
  using Line3SPtr = std::shared_ptr<Line_3>;
  using Plane3SPtr = std::shared_ptr<Plane_3>;
  using Sphere3SPtr = std::shared_ptr<Sphere_3>;

public:
  static Point3SPtr createPoint3(const FT& x, const FT& y, const FT& z)
  {
    return std::make_shared<Point_3>(x, y, z);
  }

  static Point3SPtr createPoint3(const Point_3& point)
  {
    return std::make_shared<Point_3>(point);
  }

  static Point3SPtr createPoint3(Vector3SPtr vector)
  {
    Point3SPtr result;
    result = Point3SPtr(new Point_3((*vector)[0], (*vector)[1], (*vector)[2]));
    return result;
  }

  /**
  * returns the center of the sphere.
  */
  static Point3SPtr createPoint3(Sphere3SPtr sphere)
  {
    return createPoint3(sphere->center());
  }

  static Segment3SPtr createSegment3(Point3SPtr src, Point3SPtr dst)
  {
    Segment3SPtr result = Segment3SPtr();
    if (src != dst) {
      result = Segment3SPtr(new Segment_3(*src, *dst));
    }
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Segment3SPtr createSegment3(const Segment_3& seg)
  {
    return std::make_shared<Segment_3>(seg);
  }

  static Line3SPtr createLine3(Point3SPtr p, Point3SPtr q)
  {
    Line3SPtr result = Line3SPtr();
    if (p != q) {
      result = Line3SPtr(new Line_3(*p, *q));
    }
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Line3SPtr createLine3(const Line_3& line)
  {
    return std::make_shared<Line_3>(line);
  }

  static Line3SPtr createLine3(Point3SPtr p, Vector3SPtr direction)
  {
    return std::make_shared<Line_3>(*p, *direction);
  }

  static Plane3SPtr createPlane3(const FT& a, const FT& b, const FT& c, const FT& d)
  {
    return std::make_shared<Plane_3>(a, b, c, d);
  }

  static Plane3SPtr createPlane3(Point3SPtr p, Point3SPtr q, Point3SPtr r)
  {
    Plane3SPtr result = Plane3SPtr();
    if (p != q && q != r && r != p) {
      result = Plane3SPtr(new Plane_3(*p, *q, *r));
    }
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Plane3SPtr createPlane3(Point3SPtr p, Vector3SPtr normal)
  {
    return std::make_shared<Plane_3>(*p, *normal);
  }

  static Plane3SPtr createPlane3(const Plane_3& plane)
  {
    return std::make_shared<Plane_3>(plane);
  }

  static Sphere3SPtr createSphere3(Point3SPtr center, const FT& radius)
  {
    return std::make_shared<Sphere_3>(*center, radius);
  }

  static Sphere3SPtr createSphere3(const Sphere_3& sphere)
  {
    return std::make_shared<Sphere_3>(sphere);
  }

  static Vector3SPtr createVector3(const FT& x, const FT& y, const FT& z)
  {
    return std::make_shared<Vector_3>(x, y, z);
  }

  static Vector3SPtr createVector3(const Vector_3& vector)
  {
    return std::make_shared<Vector_3>(vector);
  }

  static Vector3SPtr createVector3(Point3SPtr point)
  {
    return std::make_shared<Vector_3>(point->x(), point->y(), point->z());
  }

  /**
   * creates a vector with the same direction as the given line.
   */
  static Vector3SPtr createVector3(Line3SPtr line)
  {
    return createVector3(line->to_vector());
  }

  /**
   * returns the normal vector of the plane.
   */
  static Vector3SPtr createVector3(Plane3SPtr plane)
  {
    return createVector3(plane->orthogonal_vector());
  }
};

} // namespace kernel
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_KERNEL_FACTORY_H */
