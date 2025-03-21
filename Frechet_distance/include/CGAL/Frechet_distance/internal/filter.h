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

#pragma once
#include <CGAL/license/Frechet_distance.h>
#include <CGAL/Frechet_distance/internal/certificate.h>
#include <CGAL/Frechet_distance/internal/curve.h>
#include <CGAL/Frechet_distance/internal/geometry_basics.h>

namespace CGAL {
namespace Frechet_distance {
namespace internal {

// TODO: we can use Cartesian_converter here when we have one-sided approximate
// decisions

/*!
 * \ingroup PkgFrechetDistanceFunctions
 * A class representing a
*/
template <typename K>
class Filter
{
private:
    using distance_t = typename K::distance_t;
    using PointID = typename K::PointID;
    using Point = typename K::Point;
    using Certificate = CGAL::Frechet_distance::internal::Certificate<K>;
    using Curve = K;
    using CPoint = CGAL::Frechet_distance::internal::CPoint<K>;

    Certificate cert;
    const Curve *curve1_pt, *curve2_pt;
    typename Curve::IFT distance;

public:
    Filter(const Curve& curve1, const Curve& curve2, typename Curve::IFT const& distance)
      : curve1_pt(&curve1), curve2_pt(&curve2), distance(distance)
    {
      cert.setCurves(&curve1, &curve2);
      cert.setDistance(distance);
    }

    Certificate const& getCertificate() { return cert; };

    bool bichromaticFarthestDistance();
    bool greedy();
    bool adaptiveGreedy(PointID& pos1, PointID& pos2);
    bool adaptiveSimultaneousGreedy();
    bool negative(PointID pos1, PointID pos2);

private:
    bool isPointTooFarFromCurve(Point const& fixed, const Curve& curve);
    bool isFree(Point const& fixed, Curve const& var_curve,
                       PointID start, PointID end);
    bool isFree(Curve const& curve1, PointID start1, PointID end1,
                       Curve const& curve2, PointID start2, PointID end2);
    void increase(size_t& step);
    void decrease(size_t& step);
};

template <typename K>
bool Filter<K>::isPointTooFarFromCurve(Point const& fixed, const Curve& curve)
{
    if (possibly(Curve::distance(fixed, curve.front(), curve.traits()) <= distance) || // Uncertain (A)
        possibly(Curve::distance(fixed, curve.back(), curve.traits()) <= distance)) {
        return false;
    }
    std::size_t stepsize = 1;
    for (PointID pt = 0; pt < curve.size() - 1;) {
        stepsize = std::min<std::size_t>(stepsize, curve.size() - 1 - pt);
        auto mid = pt + (stepsize + 1) / 2;
        auto mid_dist = Curve::distance(fixed, curve[mid], curve.traits());
        auto maxdist = (CGAL::max)(curve.curve_length(pt, mid),
                                  curve.curve_length(mid, pt + stepsize));
        auto comp_dist = distance + maxdist;

        if (certainly(mid_dist > comp_dist)) { // Uncertain (A)
            pt += stepsize;
            stepsize *= 2;
        } else if (stepsize > 1) {
            stepsize /= 2;
        } else {
            return false;
        }
    }
    return true;
}

template <typename K>
bool Filter<K>::isFree(Point const& fixed, Curve const& var_curve, PointID start,
                    PointID end)
{
    auto mid = (start + end + 1) / 2;
    auto max = (CGAL::max)(var_curve.curve_length(start + 1, mid),
                           var_curve.curve_length(mid, end));

    if (certainly(distance > max)){
      auto mid_dist = Curve::distance(fixed, var_curve[mid], var_curve.traits());
      if(certainly(mid_dist <= distance - max)) { // Uncertain (A)
        return true;
      }
    }
    return false;
}

template <typename K>
bool Filter<K>::isFree(Curve const& curve1, PointID start1, PointID end1,
                       Curve const& curve2, PointID start2, PointID end2)
{
    auto mid1 = (start1 + end1 + 1) / 2;
    auto mid2 = (start2 + end2 + 1) / 2;
    auto max1 = (CGAL::max)(curve1.curve_length(start1 + 1, mid1),
                           curve1.curve_length(mid1, end1));
    auto max2 = (CGAL::max)(curve2.curve_length(start2 + 1, mid2),
                           curve2.curve_length(mid2, end2));
    auto mid_dist = Curve::distance(curve1[mid1], curve2[mid2], curve1.traits());

    // auto comp_dist = distance - max1 - max2;  // TODO: understand why if we use auto with gmpxx we get a use after free with the sanitizer (expression template I guess but why)
    typename Curve::IFT comp_dist = distance - max1 - max2;
    return certainly(! is_negative(comp_dist)) && certainly(mid_dist <= comp_dist); // Uncertain (A)
}

template <typename K>
void Filter<K>::increase(size_t& step) { step = static_cast<size_t>(std::ceil(1.5 * step)); }

template <typename K>
void Filter<K>::decrease(size_t& step) { step /= 2; }

// TODO: is it better to use #ifdef CERTIFY blocks here to avoid constructing
// CPosition objects?
// Also, calls to cert.XXX() seem to do unnecessary work right now
// uncritical, but could be optimized
template <typename K>
bool Filter<K>::bichromaticFarthestDistance()
{
    cert.reset();

    auto& curve1 = *curve1_pt;
    auto& curve2 = *curve2_pt;
    auto const& bbox1 = curve1.bbox();
    auto const& bbox2 = curve2.bbox();

    // Computes the furthest pair of points in pair of bounding boxes:
    // This can be computed coordinate wise due to the symmetry of the boxes.
    typename Curve::IFT squared_max_dist = 0;
    for (int i = 0; i < K::dimension; ++i) {
        auto d1 = CGAL::square(typename Curve::IFT((bbox1.max)(i)) - typename Curve::IFT((bbox2.min)(i)));
        auto d2 = CGAL::square(typename Curve::IFT((bbox2.max)(i)) - typename Curve::IFT((bbox1.min)(i)));
        squared_max_dist += (CGAL::max)(d1, d2);
    }

    if (certainly(sqrt(squared_max_dist) <= distance)) {
        cert.setAnswer(true);
        cert.addPoint({CPoint(0, 0), CPoint(0, 0)});
        if (curve2.size() > 1) {
            cert.addPoint({CPoint(curve1.size() - 1, 0), CPoint(0, 0)});
        }
        cert.addPoint(
            {CPoint(curve1.size() - 1, 0), CPoint(curve2.size() - 1, 0)});
        cert.validate();

        return true;
    }

    return false;
}

template <typename K>
bool Filter<K>::greedy()
{
    cert.reset();
    auto& curve1 = *curve1_pt;
    auto& curve2 = *curve2_pt;

    PointID pos1 = 0;
    PointID pos2 = 0;

    if (possibly(Curve::distance(curve1[0], curve2[0], curve1.traits()) > distance) ||  // Uncertain (A)
        possibly(Curve::distance(curve1.back(), curve2.back(), curve1.traits()) > distance)) {
        return false;
    }

    auto d_ = Curve::distance(curve1.back(), curve2.back(), curve1.traits());
    // Note that we only exit this loop if we reached the endpoints, which were
    // already checked to be close.
    while (pos1 + pos2 < curve1.size() + curve2.size() - 2) {
        d_ =
          (CGAL::max)(d_, Curve::distance(curve1[pos1], curve2[pos2], curve1.traits()));

        if (possibly(d_ > distance)) { // Uncertain (A)
            return false;
        }

        if (curve1.size() - 1 == pos1) {
            ++pos2;
        } else if (curve2.size() - 1 == pos2) {
            ++pos1;
        } else {
            auto dist1 = Curve::squared_distance(curve1[pos1 + 1], curve2[pos2], curve1.traits());
            auto dist2 = Curve::squared_distance(curve1[pos1], curve2[pos2 + 1], curve1.traits());
            auto dist12 = Curve::squared_distance(curve1[pos1 + 1], curve2[pos2 + 1], curve1.traits());

            if (possibly(dist1 < dist2) && possibly(dist1 < dist12)) { // Uncertain (A)
                ++pos1;
            } else if (possibly(dist2 < dist12)) { // Uncertain (A)
                ++pos2;
            } else {
                ++pos1;
                ++pos2;
            }
        }
    }

    return true;
}

template <typename K>
bool Filter<K>::adaptiveGreedy(PointID& pos1, PointID& pos2)
{
    cert.reset();
    auto& curve1 = *curve1_pt;
    auto& curve2 = *curve2_pt;
    // PointID pos1 = 0;
    // PointID pos2 = 0;
    pos1 = 0;
    pos2 = 0;
    cert.addPoint({CPoint(pos1, 0), CPoint(pos2, 0)});

    if (possibly(Curve::distance(curve1[0], curve2[0], curve1.traits()) > distance) ||  // Uncertain (A)
        possibly(Curve::distance(curve1.back(), curve2.back(), curve1.traits()) > distance)) {
        return false;
    }

    std::size_t step = (std::max)(curve1.size(), curve2.size());
    while (pos1 + pos2 < curve1.size() + curve2.size() - 2) {
        //++numSteps;
        // if we have to do the step on curve 2
        if (curve1.size() - 1 == pos1) {
            auto new_pos2 =
                std::min<typename PointID::IDType>(pos2 + step, curve2.size() - 1);
            if (isFree(curve1[pos1], curve2, pos2, new_pos2)) {
                pos2 = new_pos2;
                cert.addPoint({CPoint(pos1, 0), CPoint(pos2, 0)});
                // step *= 2;
                increase(step);
            } else if (step == 1) {
                return false;
            } else {
                decrease(step);
            }
        }
        // if we have to do the step on curve 1
        else if (curve2.size() - 1 == pos2) {
            auto new_pos1 =
                std::min<typename PointID::IDType>(pos1 + step, curve1.size() - 1);
            if (isFree(curve2[pos2], curve1, pos1, new_pos1)) {
                pos1 = new_pos1;
                cert.addPoint({CPoint(pos1, 0), CPoint(pos2, 0)});
                // step *= 2;
                increase(step);
            } else if (step == 1) {
                return false;
            } else {
                decrease(step);
            }
        }
        // if we cannot reduce the step size
        else if (step == 1) {
            auto dist1 = Curve::distance(curve1[pos1 + 1], curve2[pos2], curve1.traits());
            auto dist2 = Curve::distance(curve1[pos1], curve2[pos2 + 1], curve1.traits());
            auto dist12 = Curve::distance(curve1[pos1 + 1], curve2[pos2 + 1], curve1.traits());

            if (certainly(dist1 <= distance) && possibly(dist1 < dist2) && possibly(dist1 < dist12)) { // Uncertain (A)
                ++pos1;
            } else if (certainly(dist2 <= distance) && possibly(dist2 < dist12)) {
                ++pos2;
            } else if (certainly(dist12 <= distance)) {
                ++pos1;
                ++pos2;
            } else {
                return false;
            }

            cert.addPoint({CPoint(pos1, 0), CPoint(pos2, 0)});

            step = 2;
        } else {
            auto new_pos1 =
                std::min<typename PointID::IDType>(pos1 + step, curve1.size() - 1);
            auto new_pos2 =
                std::min<typename PointID::IDType>(pos2 + step, curve2.size() - 1);
            bool step1_possible =
                isFree(curve2[pos2], curve1, pos1, new_pos1);
            bool step2_possible =
                isFree(curve1[pos1], curve2, pos2, new_pos2);

            if (step1_possible && step2_possible) {
                auto dist_after_step1 =
                    Curve::squared_distance(curve1[new_pos1], curve2[pos2], curve1.traits());
                auto dist_after_step2 =
                    Curve::squared_distance(curve1[pos1], curve2[new_pos2], curve1.traits());
                if (possibly(dist_after_step1 <= dist_after_step2)) {
                    pos1 = new_pos1;
                } else {
                    pos2 = new_pos2;
                }
                cert.addPoint({CPoint(pos1, 0), CPoint(pos2, 0)});
                // step *= 2;
                increase(step);
            } else if (step1_possible) {
                pos1 = new_pos1;
                // for some reason, not increasing the stepsize here is slightly
                // faster step *= 2; increase(step);
                cert.addPoint({CPoint(pos1, 0), CPoint(pos2, 0)});
            } else if (step2_possible) {
                pos2 = new_pos2;
                // for some reason, not increasing the stepsize here is slightly
                // faster step *= 2; increase(step);
                cert.addPoint({CPoint(pos1, 0), CPoint(pos2, 0)});
            } else {
                decrease(step);
            }
        }
    }

    cert.setAnswer(true);
    cert.validate();

    return true;
}

template <typename K>
bool Filter<K>::adaptiveSimultaneousGreedy()
{
    cert.reset();
    auto& curve1 = *curve1_pt;
    auto& curve2 = *curve2_pt;

    PointID pos1 = 0;
    PointID pos2 = 0;
    cert.addPoint({CPoint(pos1, 0), CPoint(pos2, 0)});

    if (possibly(Curve::distance(curve1[0], curve2[0], curve1.traits()) > distance) || // Uncertain (A)
        possibly(Curve::distance(curve1.back(), curve2.back(), curve1.traits()) > distance)) {
        return false;
    }

    std::size_t step = (std::max)(curve1.size(), curve2.size());
    while (pos1 + pos2 < curve1.size() + curve2.size() - 2) {
        // if we have to do the step on curve 2
        if (curve1.size() - 1 == pos1) {
            auto new_pos2 =
              (std::min<typename PointID::IDType>)(pos2 + step, curve2.size() - 1);
            if (isFree(curve1[pos1], curve2, pos2, new_pos2)) {
                pos2 = new_pos2;
                cert.addPoint({CPoint(pos1, 0), CPoint(pos2, 0)});
                // step *= 2;
                increase(step);
            } else if (step == 1) {
                return false;
            } else {
                decrease(step);
            }
        }
        // if we have to do the step on curve 1
        else if (curve2.size() - 1 == pos2) {
            auto new_pos1 =
                (std::min<typename PointID::IDType>)(pos1 + step, curve1.size() - 1);
            if (isFree(curve2[pos2], curve1, pos1, new_pos1)) {
                pos1 = new_pos1;
                cert.addPoint({CPoint(pos1, 0), CPoint(pos2, 0)});
                // step *= 2;
                increase(step);
            } else if (step == 1) {
                return false;
            } else {
                decrease(step);
            }
        }
        // if we cannot reduce the step size
        else if (step == 1) {
            auto dist1 = Curve::distance(curve1[pos1 + 1], curve2[pos2], curve1.traits());
            auto dist2 = Curve::distance(curve1[pos1], curve2[pos2 + 1], curve1.traits());
            auto dist12 = Curve::distance(curve1[pos1 + 1], curve2[pos2 + 1], curve1.traits());

            if (certainly(dist1 <= distance) && possibly(dist1 < dist2) && possibly(dist1 < dist12)) { // Uncertain (A)
                ++pos1;
            } else if (certainly(dist2 <= distance) && possibly(dist2 < dist12)) {
                ++pos2;
            } else if (certainly(dist12 <= distance)) {
                ++pos1;
                ++pos2;
            } else {
                return false;
            }
            cert.addPoint({CPoint(pos1, 0), CPoint(pos2, 0)});

            step = 2;
        } else {
            PointID new_pos1 =
                std::min<typename PointID::IDType>(pos1 + step, curve1.size() - 1);
            size_t step2 =
                (step * (curve2.size() - pos2)) / (curve1.size() - pos1);
            if (step2 < 1) {
                step2 = 1;
            }
            PointID new_pos2 =
                std::min<typename PointID::IDType>(pos2 + step2, curve2.size() - 1);
            if (isFree(curve1, pos1, new_pos1, curve2, pos2, new_pos2)) {
                if (pos1 != new_pos1 && pos2 != new_pos2) {
                    cert.addPoint({CPoint(pos1 + 1, 0), CPoint(pos2 + 1, 0)});
                }
                if (new_pos1 > pos1 + 1) {
                    cert.addPoint({CPoint(new_pos1, 0), CPoint(pos2 + 1, 0)});
                }
                if (new_pos2 > pos2 + 1) {
                    cert.addPoint({CPoint(new_pos1, 0), CPoint(new_pos2, 0)});
                }
                pos1 = new_pos1;
                pos2 = new_pos2;
                increase(step);
            } else {
                decrease(step);
            }
        }
    }

    cert.setAnswer(true);
    cert.validate();

    return true;
}

template <typename K>
bool Filter<K>::negative(PointID position1, PointID position2)
{
    cert.reset();
    auto& curve1 = *curve1_pt;
    auto& curve2 = *curve2_pt;

    if (certainly(Curve::distance(curve1[0], curve2[0], curve1.traits()) > distance) || // Uncertain (A)
        certainly(Curve::distance(curve1.back(), curve2.back(), curve1.traits()) > distance)) {
        return true;
    }

    size_t pos1 = position1;
    size_t pos2 = position2;
    for (size_t step = 1; pos1 + step <= curve1.size(); increase(step)) {
        size_t cur_pos1 = pos1 + step - 1;
        if (isPointTooFarFromCurve(curve1[cur_pos1], curve2)) {
            cert.setAnswer(false);
            cert.addPoint({CPoint(cur_pos1, 0), CPoint(0, 0)});
            cert.addPoint(
                {CPoint(cur_pos1, 0), CPoint(curve2.size() - 1, 0)});
            cert.validate();
            return true;
        }
    }
    for (size_t step = 1; pos2 + step <= curve2.size(); increase(step)) {
        size_t cur_pos2 = pos2 + step - 1;
        if (isPointTooFarFromCurve(curve2[cur_pos2], curve1)) {
            cert.setAnswer(false);
            cert.addPoint(
                {CPoint(curve1.size() - 1, 0), CPoint(cur_pos2, 0)});
            cert.addPoint({CPoint(0, 0), CPoint(cur_pos2, 0)});
            cert.validate();
            return true;
        }
    }

    return false;
}

} // namespace internal
} // namespace Frechet_distance_
} // namespace CGAL
