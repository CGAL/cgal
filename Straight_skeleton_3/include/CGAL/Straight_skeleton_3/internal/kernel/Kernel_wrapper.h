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
  using Sphere_3 = typename K::Sphere_3;

  using Point3SPtr = std::shared_ptr<Point_3>;
  using Segment3SPtr = std::shared_ptr<Segment_3>;
  using Vector3SPtr = std::shared_ptr<Vector_3>;
  using Line3SPtr = std::shared_ptr<Line_3>;
  using Plane3SPtr = std::shared_ptr<Plane_3>;
  using Sphere3SPtr = std::shared_ptr<Sphere_3>;

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

  /**
    * If a line intersects a sphere, there are 2 intersection points.
    * The first one is returned here.
    */
  static Point3SPtr intersection(Sphere3SPtr sphere, Line3SPtr line)
  {
    Point3SPtr result = Point3SPtr();
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    FT radius;

    radius = CGAL::disallowed_sqrt(sphere->squared_radius());
    Vector3SPtr dir = normalize(KernelFactory::createVector3(line));
    Plane3SPtr plane = KernelFactory::createPlane3(p_center, dir);
    Point3SPtr p_intersect = intersection(plane, line);
    FT dist = distance(p_center, p_intersect);
    if (dist == radius) {
      result = p_intersect;
    } else if (dist < radius) {
      FT amount = - CGAL::disallowed_sqrt(radius*radius - dist*dist);
      result = KernelFactory::createPoint3((*p_intersect) + ((*dir)*amount));
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

  static FT distance(Point3SPtr p1, Point3SPtr p2)
  {
    return CGAL::disallowed_sqrt(CGAL::squared_distance(*p1, *p2));
  }

  static FT distance(Plane3SPtr plane, Point3SPtr point)
  {
    return CGAL::disallowed_sqrt(CGAL::squared_distance(*plane, *point));
  }

  static FT distance(Line3SPtr line, Point3SPtr point)
  {
    return CGAL::disallowed_sqrt(CGAL::squared_distance(*line, *point));
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

  static Vector3SPtr normalize(Vector3SPtr v)
  {
    Vector3SPtr result;
    result = KernelFactory::createVector3(*v / CGAL::disallowed_sqrt(v->squared_length()));
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

  static Plane3SPtr offsetPlane(Plane3SPtr plane, const FT& offset)
  {
    Plane3SPtr result = Plane3SPtr();
    CGAL_precondition(isNormalizedPlane(plane));

    const FT& a = plane->a();
    const FT& b = plane->b();
    const FT& c = plane->c();
    const FT& d = plane->d();

    CGAL_SS3_TRAITS_TRACE("Plane offset from: " << *plane);
    result = KernelFactory::createPlane3(a, b, c, d - offset);
    CGAL_SS3_TRAITS_TRACE("Plane offset to: " << *result);

    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Point3SPtr offsetPoint(Point3SPtr point, Vector3SPtr dir, const FT& offset)
  {
    Point3SPtr result;
    Vector_3 dir_normalized = *dir / CGAL::disallowed_sqrt(dir->squared_length());
    Point_3 p_moved = *point + (dir_normalized * offset);
    result = KernelFactory::createPoint3(p_moved);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Point3SPtr intersectionPointOffsetPlanes(Plane3SPtr plane_0, const FT& w0,
                                                  Plane3SPtr plane_1, const FT& w1,
                                                  Plane3SPtr plane_2, const FT& w2,
                                                  Plane3SPtr plane_3, const FT& w3)
  {
    CGAL_precondition(!(is_zero(w0) && is_zero(w1) && is_zero(w2) && is_zero(w3)));

    const FT& a0 = plane_0->a();
    const FT& b0 = plane_0->b();
    const FT& c0 = plane_0->c();
    const FT& d0 = plane_0->d();
    const FT& a1 = plane_1->a();
    const FT& b1 = plane_1->b();
    const FT& c1 = plane_1->c();
    const FT& d1 = plane_1->d();
    const FT& a2 = plane_2->a();
    const FT& b2 = plane_2->b();
    const FT& c2 = plane_2->c();
    const FT& d2 = plane_2->d();
    const FT& a3 = plane_3->a();
    const FT& b3 = plane_3->b();
    const FT& c3 = plane_3->c();
    const FT& d3 = plane_3->d();

    CGAL_SS3_TRAITS_TRACE("Coefficients\n" << a0 << " " << b0 << " " << c0 << " " << d0 << "\n"
                                           << a1 << " " << b1 << " " << c1 << " " << d1 << "\n"
                                           << a2 << " " << b2 << " " << c2 << " " << d2 << "\n"
                                           << a3 << " " << b3 << " " << c3 << " " << d3);
    CGAL_SS3_TRAITS_TRACE("Weights\n" << w0 << " " << w1 << " " << w2 << " " << w3);
    CGAL_SS3_TRAITS_TRACE("CHECK det " << CGAL::determinant(a0, b0, c0, d0,
                                                            a1, b1, c1, d1,
                                                            a2, b2, c2, d2,
                                                            a3, b3, c3, d3));

    CGAL_assertion(isNormalizedPlane(plane_0));
    CGAL_assertion(isNormalizedPlane(plane_1));
    CGAL_assertion(isNormalizedPlane(plane_2));
    CGAL_assertion(isNormalizedPlane(plane_3));

    FT den = (-a0*b1*c2*w3 + a0*b1*c3*w2 + a0*b2*c1*w3 - a0*b2*c3*w1 - a0*b3*c1*w2 + a0*b3*c2*w1 + a1*b0*c2*w3 - a1*b0*c3*w2 - a1*b2*c0*w3 + a1*b2*c3*w0 + a1*b3*c0*w2 - a1*b3*c2*w0 - a2*b0*c1*w3 + a2*b0*c3*w1 + a2*b1*c0*w3 - a2*b1*c3*w0 - a2*b3*c0*w1 + a2*b3*c1*w0 + a3*b0*c1*w2 - a3*b0*c2*w1 - a3*b1*c0*w2 + a3*b1*c2*w0 + a3*b2*c0*w1 - a3*b2*c1*w0);

    // note that below is only valid for normalized coefficients
    FT x = (b0*c1*d2*w3 - b0*c1*d3*w2 - b0*c2*d1*w3 + b0*c2*d3*w1 + b0*c3*d1*w2 - b0*c3*d2*w1 - b1*c0*d2*w3 + b1*c0*d3*w2 + b1*c2*d0*w3 - b1*c2*d3*w0 - b1*c3*d0*w2 + b1*c3*d2*w0 + b2*c0*d1*w3 - b2*c0*d3*w1 - b2*c1*d0*w3 + b2*c1*d3*w0 + b2*c3*d0*w1 - b2*c3*d1*w0 - b3*c0*d1*w2 + b3*c0*d2*w1 + b3*c1*d0*w2 - b3*c1*d2*w0 - b3*c2*d0*w1 + b3*c2*d1*w0) / den;

    FT y = (-a0*c1*d2*w3 + a0*c1*d3*w2 + a0*c2*d1*w3 - a0*c2*d3*w1 - a0*c3*d1*w2 + a0*c3*d2*w1 + a1*c0*d2*w3 - a1*c0*d3*w2 - a1*c2*d0*w3 + a1*c2*d3*w0 + a1*c3*d0*w2 - a1*c3*d2*w0 - a2*c0*d1*w3 + a2*c0*d3*w1 + a2*c1*d0*w3 - a2*c1*d3*w0 - a2*c3*d0*w1 + a2*c3*d1*w0 + a3*c0*d1*w2 - a3*c0*d2*w1 - a3*c1*d0*w2 + a3*c1*d2*w0 + a3*c2*d0*w1 - a3*c2*d1*w0) / den;

    FT z = (a0*b1*d2*w3 - a0*b1*d3*w2 - a0*b2*d1*w3 + a0*b2*d3*w1 + a0*b3*d1*w2 - a0*b3*d2*w1 - a1*b0*d2*w3 + a1*b0*d3*w2 + a1*b2*d0*w3 - a1*b2*d3*w0 - a1*b3*d0*w2 + a1*b3*d2*w0 + a2*b0*d1*w3 - a2*b0*d3*w1 - a2*b1*d0*w3 + a2*b1*d3*w0 + a2*b3*d0*w1 - a2*b3*d1*w0 - a3*b0*d1*w2 + a3*b0*d2*w1 + a3*b1*d0*w2 - a3*b1*d2*w0 - a3*b2*d0*w1 + a3*b2*d1*w0) / den;

    Point3SPtr point = KernelFactory::createPoint3(x, y, z);

    // @todo post condition that the points are at the same (weighted) time from the faces

    return point;
  }

  static FT intersectionTimeOffsetPlanes(Plane3SPtr plane_0, const FT& w0,
                                         Plane3SPtr plane_1, const FT& w1,
                                         Plane3SPtr plane_2, const FT& w2,
                                         Plane3SPtr plane_3, const FT& w3,
                                         const std::optional<FT>& past_bound = std::nullopt,
                                         const std::optional<FT>& future_bound = std::nullopt)
  {
    CGAL_precondition(!(is_zero(w0) && is_zero(w1) && is_zero(w2) && is_zero(w3)));

    const FT& a0 = plane_0->a();
    const FT& b0 = plane_0->b();
    const FT& c0 = plane_0->c();
    const FT& d0 = plane_0->d();
    const FT& a1 = plane_1->a();
    const FT& b1 = plane_1->b();
    const FT& c1 = plane_1->c();
    const FT& d1 = plane_1->d();
    const FT& a2 = plane_2->a();
    const FT& b2 = plane_2->b();
    const FT& c2 = plane_2->c();
    const FT& d2 = plane_2->d();
    const FT& a3 = plane_3->a();
    const FT& b3 = plane_3->b();
    const FT& c3 = plane_3->c();
    const FT& d3 = plane_3->d();

    CGAL_SS3_TRAITS_TRACE("Coefficients\n" << a0 << " " << b0 << " " << c0 << " " << d0 << "\n"
                                           << a1 << " " << b1 << " " << c1 << " " << d1 << "\n"
                                           << a2 << " " << b2 << " " << c2 << " " << d2 << "\n"
                                           << a3 << " " << b3 << " " << c3 << " " << d3);
    CGAL_SS3_TRAITS_TRACE("Weights\n" << w0 << " " << w1 << " " << w2 << " " << w3);
    CGAL_SS3_TRAITS_TRACE("CHECK det " << CGAL::determinant(a0, b0, c0, d0,
                                                            a1, b1, c1, d1,
                                                            a2, b2, c2, d2,
                                                            a3, b3, c3, d3));

    CGAL_assertion(isNormalizedPlane(plane_0));
    CGAL_assertion(isNormalizedPlane(plane_1));
    CGAL_assertion(isNormalizedPlane(plane_2));
    CGAL_assertion(isNormalizedPlane(plane_3));

    FT den = (-a0*b1*c2*w3 + a0*b1*c3*w2 + a0*b2*c1*w3 - a0*b2*c3*w1 - a0*b3*c1*w2 + a0*b3*c2*w1 + a1*b0*c2*w3 - a1*b0*c3*w2 - a1*b2*c0*w3 + a1*b2*c3*w0 + a1*b3*c0*w2 - a1*b3*c2*w0 - a2*b0*c1*w3 + a2*b0*c3*w1 + a2*b1*c0*w3 - a2*b1*c3*w0 - a2*b3*c0*w1 + a2*b3*c1*w0 + a3*b0*c1*w2 - a3*b0*c2*w1 - a3*b1*c0*w2 + a3*b1*c2*w0 + a3*b2*c0*w1 - a3*b2*c1*w0);

    FT t = (-a0*b1*c2*d3 + a0*b1*c3*d2 + a0*b2*c1*d3 - a0*b2*c3*d1 - a0*b3*c1*d2 + a0*b3*c2*d1 + a1*b0*c2*d3 - a1*b0*c3*d2 - a1*b2*c0*d3 + a1*b2*c3*d0 + a1*b3*c0*d2 - a1*b3*c2*d0 - a2*b0*c1*d3 + a2*b0*c3*d1 + a2*b1*c0*d3 - a2*b1*c3*d0 - a2*b3*c0*d1 + a2*b3*c1*d0 + a3*b0*c1*d2 - a3*b0*c2*d1 - a3*b1*c0*d2 + a3*b1*c2*d0 + a3*b2*c0*d1 - a3*b2*c1*d0) / den;

    if (past_bound && t >= *past_bound)  {
      CGAL_SS3_TRAITS_TRACE("Event is strictly in the past");
      return { };
    }

    if (future_bound && t <= *future_bound) {
      CGAL_SS3_TRAITS_TRACE("Event is too far in the future");
      return { };
    }

    return t;
  }

  static std::pair<Point3SPtr, FT> intersectionPointAndTimeOffsetPlanes(Plane3SPtr plane_0, const FT& w0,
                                                                        Plane3SPtr plane_1, const FT& w1,
                                                                        Plane3SPtr plane_2, const FT& w2,
                                                                        Plane3SPtr plane_3, const FT& w3,
                                                                        const std::optional<FT>& past_bound = std::nullopt,
                                                                        const std::optional<FT>& future_bound = std::nullopt)
  {
    CGAL_precondition(!(is_zero(w0) && is_zero(w1) && is_zero(w2) && is_zero(w3)));

    const FT& a0 = plane_0->a();
    const FT& b0 = plane_0->b();
    const FT& c0 = plane_0->c();
    const FT& d0 = plane_0->d();
    const FT& a1 = plane_1->a();
    const FT& b1 = plane_1->b();
    const FT& c1 = plane_1->c();
    const FT& d1 = plane_1->d();
    const FT& a2 = plane_2->a();
    const FT& b2 = plane_2->b();
    const FT& c2 = plane_2->c();
    const FT& d2 = plane_2->d();
    const FT& a3 = plane_3->a();
    const FT& b3 = plane_3->b();
    const FT& c3 = plane_3->c();
    const FT& d3 = plane_3->d();

    CGAL_SS3_TRAITS_TRACE("Coefficients\n" << a0 << " " << b0 << " " << c0 << " " << d0 << "\n"
                                           << a1 << " " << b1 << " " << c1 << " " << d1 << "\n"
                                           << a2 << " " << b2 << " " << c2 << " " << d2 << "\n"
                                           << a3 << " " << b3 << " " << c3 << " " << d3);
    CGAL_SS3_TRAITS_TRACE("Weights\n" << w0 << " " << w1 << " " << w2 << " " << w3);
    CGAL_SS3_TRAITS_TRACE("CHECK det " << CGAL::determinant(a0, b0, c0, d0,
                                                            a1, b1, c1, d1,
                                                            a2, b2, c2, d2,
                                                            a3, b3, c3, d3));

    CGAL_assertion(isNormalizedPlane(plane_0));
    CGAL_assertion(isNormalizedPlane(plane_1));
    CGAL_assertion(isNormalizedPlane(plane_2));
    CGAL_assertion(isNormalizedPlane(plane_3));

    FT tn = (-a0*b1*c2*d3 + a0*b1*c3*d2 + a0*b2*c1*d3 - a0*b2*c3*d1 - a0*b3*c1*d2 + a0*b3*c2*d1 + a1*b0*c2*d3 - a1*b0*c3*d2 - a1*b2*c0*d3 + a1*b2*c3*d0 + a1*b3*c0*d2 - a1*b3*c2*d0 - a2*b0*c1*d3 + a2*b0*c3*d1 + a2*b1*c0*d3 - a2*b1*c3*d0 - a2*b3*c0*d1 + a2*b3*c1*d0 + a3*b0*c1*d2 - a3*b0*c2*d1 - a3*b1*c0*d2 + a3*b1*c2*d0 + a3*b2*c0*d1 - a3*b2*c1*d0);

    FT den = (-a0*b1*c2*w3 + a0*b1*c3*w2 + a0*b2*c1*w3 - a0*b2*c3*w1 - a0*b3*c1*w2 + a0*b3*c2*w1 + a1*b0*c2*w3 - a1*b0*c3*w2 - a1*b2*c0*w3 + a1*b2*c3*w0 + a1*b3*c0*w2 - a1*b3*c2*w0 - a2*b0*c1*w3 + a2*b0*c3*w1 + a2*b1*c0*w3 - a2*b1*c3*w0 - a2*b3*c0*w1 + a2*b3*c1*w0 + a3*b0*c1*w2 - a3*b0*c2*w1 - a3*b1*c0*w2 + a3*b1*c2*w0 + a3*b2*c0*w1 - a3*b2*c1*w0);

    CGAL::Sign s = CGAL::sign(den);

    // Bound checks. The algorithm works is a decreasing t, so past is greater than future.
    //
    // Empirically:
    // - It's about as likely to be greater than 'past' than it is to be lower than 'future'
    // - Avoiding the division before the checks does not yield observable gains
    // - Adding a check past_bound == 0 to simply check tn and den signs also gains nothing...
    //
    // t >= past_bound
    // tn/den - past_bound >= 0
    //   | tn - den*past_bound >= 0 if den > 0 => + (tn - den*past_bound) >= 0
    //   | tn - den*past_bound <= 0 if den < 0 => - (tn - den*past_bound) >= 0
    // sign(den) * (tn - den*past_bound) >= 0
    if (past_bound && !CGAL::is_negative(s * (tn - *past_bound * den))) {
      CGAL_SS3_TRAITS_TRACE("Event is strictly in the past");
      return { };
    }

    if (future_bound && !CGAL::is_positive(s * (tn - *future_bound * den))) {
      CGAL_SS3_TRAITS_TRACE("Event is too far in the future");
      return { };
    }

    FT t = tn / den;

    FT x = (b0*c1*d2*w3 - b0*c1*d3*w2 - b0*c2*d1*w3 + b0*c2*d3*w1 + b0*c3*d1*w2 - b0*c3*d2*w1 - b1*c0*d2*w3 + b1*c0*d3*w2 + b1*c2*d0*w3 - b1*c2*d3*w0 - b1*c3*d0*w2 + b1*c3*d2*w0 + b2*c0*d1*w3 - b2*c0*d3*w1 - b2*c1*d0*w3 + b2*c1*d3*w0 + b2*c3*d0*w1 - b2*c3*d1*w0 - b3*c0*d1*w2 + b3*c0*d2*w1 + b3*c1*d0*w2 - b3*c1*d2*w0 - b3*c2*d0*w1 + b3*c2*d1*w0) / den;

    FT y = (-a0*c1*d2*w3 + a0*c1*d3*w2 + a0*c2*d1*w3 - a0*c2*d3*w1 - a0*c3*d1*w2 + a0*c3*d2*w1 + a1*c0*d2*w3 - a1*c0*d3*w2 - a1*c2*d0*w3 + a1*c2*d3*w0 + a1*c3*d0*w2 - a1*c3*d2*w0 - a2*c0*d1*w3 + a2*c0*d3*w1 + a2*c1*d0*w3 - a2*c1*d3*w0 - a2*c3*d0*w1 + a2*c3*d1*w0 + a3*c0*d1*w2 - a3*c0*d2*w1 - a3*c1*d0*w2 + a3*c1*d2*w0 + a3*c2*d0*w1 - a3*c2*d1*w0) / den;

    FT z = (a0*b1*d2*w3 - a0*b1*d3*w2 - a0*b2*d1*w3 + a0*b2*d3*w1 + a0*b3*d1*w2 - a0*b3*d2*w1 - a1*b0*d2*w3 + a1*b0*d3*w2 + a1*b2*d0*w3 - a1*b2*d3*w0 - a1*b3*d0*w2 + a1*b3*d2*w0 + a2*b0*d1*w3 - a2*b0*d3*w1 - a2*b1*d0*w3 + a2*b1*d3*w0 + a2*b3*d0*w1 - a2*b3*d1*w0 - a3*b0*d1*w2 + a3*b0*d2*w1 + a3*b1*d0*w2 - a3*b1*d2*w0 - a3*b2*d0*w1 + a3*b2*d1*w0) / den;

    Point3SPtr result = KernelFactory::createPoint3(x, y, z);

    CGAL_SS3_TRAITS_TRACE("CHECK x|y|z|t " << x << " " << y << " " << z << " " << t);
    CGAL_SS3_TRAITS_TRACE("CHECK 0: " << a0*x + b0*y + c0*z + d0 - w0*t);
    CGAL_SS3_TRAITS_TRACE("CHECK 1: " << a1*x + b1*y + c1*z + d1 - w1*t);
    CGAL_SS3_TRAITS_TRACE("CHECK 2: " << a2*x + b2*y + c2*z + d2 - w2*t);
    CGAL_SS3_TRAITS_TRACE("CHECK 3: " << a3*x + b3*y + c3*z + d3 - w3*t);

    CGAL_postcondition(a0*x + b0*y + c0*z + d0 - w0*t == 0);
    CGAL_postcondition(a1*x + b1*y + c1*z + d1 - w1*t == 0);
    CGAL_postcondition(a2*x + b2*y + c2*z + d2 - w2*t == 0);
    CGAL_postcondition(a3*x + b3*y + c3*z + d3 - w3*t == 0);

    return { result, t };
  }

  /**
    * http://de.wikipedia.org/wiki/Drehmatrix
    *
    *              [ n_x^2 (1 - cos(alpha)) + cos(alpha)         n_x n_y (1 - cos(alpha)) - n_z sin(alpha)   n_x n_z (1 - cos(alpha)) + n_y sin(alpha) ]
    * R_n(alpha) = [ n_y n_x (1 - cos(alpha)) + n_z sin(alpha)   n_y^2 (1 - cos(alpha)) + cos(alpha)         n_y n_z (1 - cos(alpha)) - n_x sin(alpha) ]
    *              [ n_z n_x (1 - cos(alpha)) - n_y sin(alpha)   n_z n_y (1 - cos(alpha)) + n_x sin(alpha)   n_z^2 (1 - cos(alpha)) + cos(alpha)       ]
    */
  static Vector3SPtr rotateVector(Vector3SPtr vector, Vector3SPtr axis, const FT& angle)
  {
    Vector3SPtr result;
    Vector3SPtr v_n = KernelWrapper::normalize(axis);
    FT n[3];
    for (unsigned int i = 0; i < 3; ++i) {
      n[i] = (*v_n)[i];
    }

    FT cos_angle = std::cos(CGAL::to_double(angle));
    FT sin_angle = std::sin(CGAL::to_double(angle));

    FT rotation[3][3];   // http://de.wikipedia.org/wiki/Drehmatrix
    rotation[0][0] = n[0]*n[0] * (1.0-cos_angle) + cos_angle;
    rotation[0][1] = n[0]*n[1] * (1.0-cos_angle) - n[2] * sin_angle;
    rotation[0][2] = n[0]*n[2] * (1.0-cos_angle) + n[1] * sin_angle;
    rotation[1][0] = n[1]*n[0] * (1.0-cos_angle) + n[2] * sin_angle;
    rotation[1][1] = n[1]*n[1] * (1.0-cos_angle) + cos_angle;
    rotation[1][2] = n[1]*n[2] * (1.0-cos_angle) - n[0] * sin_angle;
    rotation[2][0] = n[2]*n[0] * (1.0-cos_angle) - n[1] * sin_angle;
    rotation[2][1] = n[2]*n[1] * (1.0-cos_angle) + n[0] * sin_angle;
    rotation[2][2] = n[2]*n[2] * (1.0-cos_angle) + cos_angle;

    FT rotated[3];
    for (unsigned int r = 0; r < 3; r++) {
      rotated[r] = 0.0;
      for (unsigned int c = 0; c < 3; c++) {
        rotated[r] += rotation[r][c] * (*vector)[c];
      }
    }
    result = KernelFactory::createVector3(rotated[0], rotated[1], rotated[2]);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Plane3SPtr rotatePlane(Plane3SPtr plane, Line3SPtr line, const FT& angle)
  {
    Plane3SPtr result;
    Point3SPtr point;
    point = KernelFactory::createPoint3(line->point(0));
    Vector3SPtr dir = KernelFactory::createVector3(line);
    Vector3SPtr normal = KernelFactory::createVector3(plane);
    Vector3SPtr normal_rotated = rotateVector(normal, dir, angle);
    result = KernelFactory::createPlane3(point, normal_rotated);
    CGAL_SS3_DEBUG_SPTR(result);
    return result;
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

  static int comparePoints(Vector3SPtr v_dir, Point3SPtr p_1, Point3SPtr p_2)
  {
    int result = 0;
    FT value = *v_dir * (*p_2 - *p_1);
    if (value > 0.0) {         // angle < CGAL_PI/2.0
      result = -1;
    } else if (value < 0.0) {  // angle > CGAL_PI/2.0
      result = 1;
    }
    return result;
  }

  static Point3SPtr replaceCoord(Point3SPtr point,
                                 Point3SPtr replacement,
                                 unsigned int coord)
  {
    Point3SPtr result = Point3SPtr();

    if (coord == 0) {
      result = Point3SPtr(new Point_3(replacement->x(), point->y(), point->z()));
    } else if (coord == 1) {
      result = Point3SPtr(new Point_3(point->x(), replacement->y(), point->z()));
    } else if (coord == 2) {
      result = Point3SPtr(new Point_3(point->x(), point->y(), replacement->z()));
    }

    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }
};

} // namespace kernel
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_KERNEL_KERNEL_WRAPPER_H */
