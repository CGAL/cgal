// Copyright (c) 2024 Max-Planck-Institute Saarbruecken (Germany), GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//                 Karl Bringmann <kbringma@mpi-inf.mpg.de>
//                 Andreas Fabri
// =============================================================================

#ifndef CGAL_FRECHET_DISTANCE_INTERNAL_CERTIFICATE_H
#define CGAL_FRECHET_DISTANCE_INTERNAL_CERTIFICATE_H
#include <CGAL/license/Frechet_distance.h>
#include <CGAL/Frechet_distance/internal/Curve.h>
#include <CGAL/Frechet_distance/internal/geometry_basics.h>

namespace CGAL {
namespace Frechet_distance {
namespace internal {

template <typename C>
class Certificate
{
    using Curve = C;
    using CPositions = CGAL::Frechet_distance::internal::CPositions<C>;
    using CPosition = CGAL::Frechet_distance::internal::CPosition<C>;
    using CPoint = CGAL::Frechet_distance::internal::CPoint<C>;
public:
    Certificate() = default;

    bool isYes() const
    {
        assert(isValid());
        return lessThan;
    }

    bool isValid() const { return valid; }

    bool check() const;
    const CPositions& getTraversal() const { return traversal; }

    void addPoint(const CPosition& pos) { traversal.push_back(pos); }
    void setAnswer(bool answer) { lessThan = answer; }
    void setCurves(const Curve* curve1, const Curve* curve2)
    {
        curve_pair[0] = curve1;
        curve_pair[1] = curve2;
    }
    void setDistance(typename Curve::IFT distance)
    {
        dist = distance;
    }
    void validate() { valid = true; };
    void reset()
    {
        valid = false;
        traversal.clear();
    }
    void clear() { traversal.clear(); }

    void dump_certificate() const;

private:
    // if YES certificate: (traversal1[0], traversal2[0]), ..., (traversal1[T],
    // traversal2[T]) should be feasible traversal of curve1 and curve2
    // (interpolate between discrete points) if NO certificate: certificate
    // should be as described in certificate-outline.pdf

    CPositions  traversal;
    std::array<const Curve*, 2> curve_pair;
    typename Curve::IFT dist;

    bool lessThan;
    bool valid = false;

    bool feasible(const CPosition& pt) const;
    bool feasible(const CPoint& pt1, const CPoint& pt2) const;
    bool nonEmpty(CurveID fixed_curve, const CPoint& fixed_point,
                  const CPoint& start_pt, const CPoint& end_point) const;
};

} // namespace internal
} // namespace Frechet_distance_
} // namespace CGAL

#endif // CGAL_FRECHET_DISTANCE_INTERNAL_CERTIFICATE_H
