// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_GEOM_UTILS_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_GEOM_UTILS_H

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_factory.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_wrapper.h>

#include <CGAL/assertions.h>
#include <CGAL/enum.h>
#include <CGAL/determinant.h>
#include <CGAL/number_utils.h>

#include <cmath>
#include <optional>
#include <utility>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename K>
class GeomUtils
{
  using FT = typename K::FT;
  using Point_3 = typename K::Point_3;
  using Plane_3 = typename K::Plane_3;

  using Point3SPtr = std::shared_ptr<Point_3>;
  using Plane3SPtr = std::shared_ptr<Plane_3>;

  using KernelWrapper = kernel::KernelWrapper<K>;
  using KernelFactory = kernel::KernelFactory<K>;

public:
  static Plane3SPtr offsetPlane(const Plane3SPtr& plane, const FT& shift)
  {
    CGAL_SS3_DEBUG_SPTR(plane);
    CGAL_precondition(KernelWrapper::hasNormalizedPlane(plane));

    CGAL_SS3_TRAITS_TRACE("Plane offset from: " << *plane);

    const FT& a = plane->a();
    const FT& b = plane->b();
    const FT& c = plane->c();
    const FT& d = plane->d();
    Plane3SPtr result = KernelFactory::createPlane3(a, b, c, d - shift);
    CGAL_SS3_TRAITS_TRACE("Plane offset to: " << *result);

    CGAL_SS3_DEBUG_SPTR(result);
    return result;
  }

  static Point3SPtr intersectionPointOffsetPlanes(const Plane3SPtr& plane_0, const FT& w0,
                                                  const Plane3SPtr& plane_1, const FT& w1,
                                                  const Plane3SPtr& plane_2, const FT& w2,
                                                  const Plane3SPtr& plane_3, const FT& w3)
  {
    CGAL_SS3_DEBUG_SPTR(plane_0);
    CGAL_SS3_DEBUG_SPTR(plane_1);
    CGAL_SS3_DEBUG_SPTR(plane_2);
    CGAL_SS3_DEBUG_SPTR(plane_3);
    CGAL_precondition(!is_zero(w0) && !is_zero(w1) && !is_zero(w2) && !is_zero(w3));

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

    CGAL_assertion(KernelWrapper::hasNormalizedPlane(plane_0));
    CGAL_assertion(KernelWrapper::hasNormalizedPlane(plane_1));
    CGAL_assertion(KernelWrapper::hasNormalizedPlane(plane_2));
    CGAL_assertion(KernelWrapper::hasNormalizedPlane(plane_3));

    FT den = (-a0*b1*c2*w3 + a0*b1*c3*w2 + a0*b2*c1*w3 - a0*b2*c3*w1 - a0*b3*c1*w2 + a0*b3*c2*w1 + a1*b0*c2*w3 - a1*b0*c3*w2 - a1*b2*c0*w3 + a1*b2*c3*w0 + a1*b3*c0*w2 - a1*b3*c2*w0 - a2*b0*c1*w3 + a2*b0*c3*w1 + a2*b1*c0*w3 - a2*b1*c3*w0 - a2*b3*c0*w1 + a2*b3*c1*w0 + a3*b0*c1*w2 - a3*b0*c2*w1 - a3*b1*c0*w2 + a3*b1*c2*w0 + a3*b2*c0*w1 - a3*b2*c1*w0);

    // note that below is only valid for normalized coefficients
    FT x = (b0*c1*d2*w3 - b0*c1*d3*w2 - b0*c2*d1*w3 + b0*c2*d3*w1 + b0*c3*d1*w2 - b0*c3*d2*w1 - b1*c0*d2*w3 + b1*c0*d3*w2 + b1*c2*d0*w3 - b1*c2*d3*w0 - b1*c3*d0*w2 + b1*c3*d2*w0 + b2*c0*d1*w3 - b2*c0*d3*w1 - b2*c1*d0*w3 + b2*c1*d3*w0 + b2*c3*d0*w1 - b2*c3*d1*w0 - b3*c0*d1*w2 + b3*c0*d2*w1 + b3*c1*d0*w2 - b3*c1*d2*w0 - b3*c2*d0*w1 + b3*c2*d1*w0) / den;

    FT y = (-a0*c1*d2*w3 + a0*c1*d3*w2 + a0*c2*d1*w3 - a0*c2*d3*w1 - a0*c3*d1*w2 + a0*c3*d2*w1 + a1*c0*d2*w3 - a1*c0*d3*w2 - a1*c2*d0*w3 + a1*c2*d3*w0 + a1*c3*d0*w2 - a1*c3*d2*w0 - a2*c0*d1*w3 + a2*c0*d3*w1 + a2*c1*d0*w3 - a2*c1*d3*w0 - a2*c3*d0*w1 + a2*c3*d1*w0 + a3*c0*d1*w2 - a3*c0*d2*w1 - a3*c1*d0*w2 + a3*c1*d2*w0 + a3*c2*d0*w1 - a3*c2*d1*w0) / den;

    FT z = (a0*b1*d2*w3 - a0*b1*d3*w2 - a0*b2*d1*w3 + a0*b2*d3*w1 + a0*b3*d1*w2 - a0*b3*d2*w1 - a1*b0*d2*w3 + a1*b0*d3*w2 + a1*b2*d0*w3 - a1*b2*d3*w0 - a1*b3*d0*w2 + a1*b3*d2*w0 + a2*b0*d1*w3 - a2*b0*d3*w1 - a2*b1*d0*w3 + a2*b1*d3*w0 + a2*b3*d0*w1 - a2*b3*d1*w0 - a3*b0*d1*w2 + a3*b0*d2*w1 + a3*b1*d0*w2 - a3*b1*d2*w0 - a3*b2*d0*w1 + a3*b2*d1*w0) / den;

    Point3SPtr point = KernelFactory::createPoint3(x, y, z);

    return point;
  }

  static FT intersectionTimeOffsetPlanes(const Plane3SPtr& plane_0, const FT& w0,
                                         const Plane3SPtr& plane_1, const FT& w1,
                                         const Plane3SPtr& plane_2, const FT& w2,
                                         const Plane3SPtr& plane_3, const FT& w3)
  {
    CGAL_SS3_DEBUG_SPTR(plane_0);
    CGAL_SS3_DEBUG_SPTR(plane_1);
    CGAL_SS3_DEBUG_SPTR(plane_2);
    CGAL_SS3_DEBUG_SPTR(plane_3);
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

    CGAL_assertion(KernelWrapper::hasNormalizedPlane(plane_0));
    CGAL_assertion(KernelWrapper::hasNormalizedPlane(plane_1));
    CGAL_assertion(KernelWrapper::hasNormalizedPlane(plane_2));
    CGAL_assertion(KernelWrapper::hasNormalizedPlane(plane_3));

    FT den = (-a0*b1*c2*w3 + a0*b1*c3*w2 + a0*b2*c1*w3 - a0*b2*c3*w1 - a0*b3*c1*w2 + a0*b3*c2*w1 + a1*b0*c2*w3 - a1*b0*c3*w2 - a1*b2*c0*w3 + a1*b2*c3*w0 + a1*b3*c0*w2 - a1*b3*c2*w0 - a2*b0*c1*w3 + a2*b0*c3*w1 + a2*b1*c0*w3 - a2*b1*c3*w0 - a2*b3*c0*w1 + a2*b3*c1*w0 + a3*b0*c1*w2 - a3*b0*c2*w1 - a3*b1*c0*w2 + a3*b1*c2*w0 + a3*b2*c0*w1 - a3*b2*c1*w0);

    FT tn = (-a0*b1*c2*d3 + a0*b1*c3*d2 + a0*b2*c1*d3 - a0*b2*c3*d1 - a0*b3*c1*d2 + a0*b3*c2*d1 + a1*b0*c2*d3 - a1*b0*c3*d2 - a1*b2*c0*d3 + a1*b2*c3*d0 + a1*b3*c0*d2 - a1*b3*c2*d0 - a2*b0*c1*d3 + a2*b0*c3*d1 + a2*b1*c0*d3 - a2*b1*c3*d0 - a2*b3*c0*d1 + a2*b3*c1*d0 + a3*b0*c1*d2 - a3*b0*c2*d1 - a3*b1*c0*d2 + a3*b1*c2*d0 + a3*b2*c0*d1 - a3*b2*c1*d0);

    return tn / den;
  }
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_GEOM_UTILS_H */
