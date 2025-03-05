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
#include <CGAL/Frechet_distance/internal/filter.h>
#include <CGAL/Frechet_distance/internal/frechet_light_types.h>
#include <CGAL/Frechet_distance/internal/geometry_basics.h>
#include <CGAL/Frechet_distance/internal/id.h>
#include <CGAL/Frechet_distance/internal/high_level_predicates.h>


#include <algorithm>
#include <array>
#include <vector>

namespace CGAL {
namespace Frechet_distance {
namespace internal {


/*!
 * \ingroup PkgFrechetDistanceFunctions
 * A class representing a
*/
  template <typename C>
class FrechetLight
{
    using Curve = C;
//    using K = typename Curve::K;
    using Point = typename C::Point;
    using PointID = typename Curve::PointID;
    using distance_t = typename Curve::distance_t;
    using CPoint = CGAL::Frechet_distance::internal::CPoint<C>;
    using CInterval = CGAL::Frechet_distance::internal::CInterval<C>;
    using CIntervals = CGAL::Frechet_distance::internal::CIntervals<C>;
    using CIntervalsID = CGAL::Frechet_distance::internal::CIntervalsID<C>;
    using CPosition = CGAL::Frechet_distance::internal::CPosition<C>;
    using QSimpleInterval = CGAL::Frechet_distance::internal::QSimpleInterval<C>;
    using QSimpleIntervals = CGAL::Frechet_distance::internal::QSimpleIntervals<C>;
    using QSimpleOutputs = CGAL::Frechet_distance::internal::QSimpleOutputs<C>;
    using Certificate = CGAL::Frechet_distance::internal::Certificate<C>;
    using Filter = CGAL::Frechet_distance::internal::Filter<C>;
    using Inputs = CGAL::Frechet_distance::internal::Inputs<C>;
    using Outputs = CGAL::Frechet_distance::internal::Outputs<C>;
    using BoxData = CGAL::Frechet_distance::internal::BoxData<C>;
    using Lambda = CGAL::Frechet_distance::internal::Lambda<C>;
    using CurvePair = std::array<Curve const*, 2>;

public:
    FrechetLight() = default;
    void buildFreespaceDiagram(distance_t const& distance, Curve const& curve1,
                               Curve const& curve2);
    bool lessThan(distance_t const& distance, Curve const& curve1,
                  Curve const& curve2);
    bool lessThanWithFilters(distance_t const& distance, Curve const& curve1,
                             Curve const& curve2);
    std::pair<double,double> calcDistance(Curve const& curve1, Curve const& curve2, double epsilon);
    void clear();

    CurvePair getCurvePair() const;
    Certificate& computeCertificate();
    const Certificate& getCertificate() const { return cert; }

    void setPruningLevel(int pruning_level);
    void setRules(std::array<bool, 5> const& enable);

    std::size_t getNumberOfBoxes() const;

    std::size_t non_filtered = 0;

private:
    CurvePair curve_pair;

    typename Curve::IFT distance;

    std::vector<CIntervals> reachable_intervals_vec;
    QSimpleIntervals qsimple_intervals;
    std::size_t num_boxes;

    // 0 = no pruning ... 6 = full pruning
    int pruning_level = 6;
    // ... and additionally bools to enable/disable rules
    bool enable_box_shrinking = true;
    bool enable_empty_outputs = true;
    bool enable_propagation1 = true;
    bool enable_propagation2 = true;
    bool enable_boundary_rule = true;

#ifdef VIS
    CIntervals unknown_intervals;
    CIntervals connections;
    CIntervals free_non_reachable;
    CIntervals reachable_intervals;

    struct Cell {
        PointID i, j;
        Cell(PointID i, PointID j) : i(i), j(j) {}
    };

    std::vector<Cell> cells;
    std::vector<Cell> const& getCells() const { return cells; }
#endif

    Certificate cert;

  template <class IndexType>
  CInterval getInterval(Curve const& curve1,
                        IndexType const& center_id,
                        Curve const& curve2,
                        PointID seg_start) const
  {
    if constexpr (!std::is_same_v<IndexType, PointID>)
    {
      assert(false);
      return CInterval();
    }
    else
    {
      auto interval = HLPred::intersection_interval(curve1, center_id, curve2, seg_start, distance);
        return CInterval{seg_start, interval.begin, seg_start, interval.end};
    }
  }

    void merge(CIntervals& v, CInterval const& i) const;

    Outputs createFinalOutputs();
    Inputs computeInitialInputs();
    bool isClose(Curve const& curve1, PointID i, Curve const& curve) const;
    CPoint getLastReachablePoint(Curve const& curve1, PointID i,
                                 Curve const& curve2) const;
    bool isTopRightReachable(Outputs const& outputs) const;
    void computeOutputs(Box<C> const& initial_box, Inputs const& initial_inputs,
                        Outputs& final_outputs);

    void getReachableIntervals(BoxData& data);

    // subfunctions of getReachableIntervals
    bool emptyInputsRule(BoxData& data);
    void boxShrinkingRule(BoxData& data);
    void handleCellCase(BoxData& data);
    void getQSimpleIntervals(BoxData& data);
    void calculateQSimple1(BoxData& data);
    void calculateQSimple2(BoxData& data);
    bool boundaryPruningRule(BoxData& data);
    void splitAndRecurse(BoxData& data);

    // intervals used in getReachableIntervals and subfunctions
    CInterval const empty;
    CInterval const* firstinterval1;
    CInterval const* firstinterval2;
    Lambda min1_frac, min2_frac;
    QSimpleInterval qsimple1, qsimple2;
    CInterval out1, out2;
    // TODO: can those be made members of out1, out2?
    bool out1_valid = false, out2_valid = false;

    // qsimple interval calculation functions
    QSimpleInterval getTestFullQSimpleInterval(CPoint const& fixed,
                                               const Curve& fixed_curve,
                                               PointID min, PointID max,
                                               const Curve& curve) const;
    // bool updateQSimpleInterval(QSimpleInterval& qsimple, PointID const&
    // fixed, const Curve& fixed_curve, PointID min1, PointID max1, const Curve&
    // curve) const; internal use:
    class TestFullMode;  // mode to test whether simple interval is fill
    class RegularMode;   // mode to obtain complete information
    template <class MODE = RegularMode, typename IndexType = PointID>
    bool updateQSimpleInterval(QSimpleInterval& qsimple, const IndexType& fixed,
                               const Curve& fixed_curve, PointID min1,
                               PointID max1, const Curve& curve) const;
    template <class MODE = RegularMode, typename IndexType = PointID>
    void continueQSimpleSearch(QSimpleInterval& qsimple, const IndexType& fixed,
                               const Curve& fixed_curve, PointID min1,
                               PointID max1, const Curve& curve) const;

    bool isOnLowerRight(const CPosition& pt) const;
    bool isOnUpperLeft(const CPosition& pt) const;

    void initCertificate(Inputs const& initial_inputs);
    void certSetValues(CInterval& interval, CInterval const& parent,
                       PointID point_id, CurveID curve_id);
    void certAddEmpty(CPoint begin, CPoint end, CPoint fixed_point,
                      CurveID fixed_curve);

    // Those are empty function if VIS is not defined
    void visAddReachable(CInterval const& cinterval);
    void visAddUnknown(CPoint begin, CPoint end, CPoint fixed_point,
                       CurveID fixed_curve);
    void visAddConnection(CPoint begin, CPoint end, CPoint fixed_point,
                          CurveID fixed_curve);
    void visAddFreeNonReachable(CPoint begin, CPoint end, CPoint fixed_point,
                                CurveID fixed_curve);
    void visAddCell(Box<C> const& box);

    // Could also be done via getter member functions, but vis is a special
    // case in needing access to the internal structures.
    friend class FreespaceLightVis;
};

template <typename C>
void FrechetLight<C>::certSetValues(CInterval& interval, CInterval const& parent,
                                 PointID point_id, CurveID curve_id)
{
    interval.reach_parent = &parent;
    interval.fixed = CPoint(point_id, 0);
    interval.fixed_curve = curve_id;
}

template <typename C>
void FrechetLight<C>::visAddReachable(CInterval const& cinterval)
{
#ifdef VIS
    if (cinterval.is_empty()) {
        return;
    }
    reachable_intervals.push_back(cinterval);
#else
  CGAL_USE(cinterval);
#endif
}

template <typename C>
void FrechetLight<C>::visAddUnknown(CPoint begin, CPoint end, CPoint fixed_point,
                                 CurveID fixed_curve)
{
#ifdef VIS
    if (begin >= end) {
        return;
    }

    assert(begin >= 0);
    assert(end <= curve_pair[1 - fixed_curve]->size() - 1);

    unknown_intervals.emplace_back(begin, end);
    unknown_intervals.back().fixed = fixed_point;
    unknown_intervals.back().fixed_curve = fixed_curve;

    assert(unknown_intervals.back().begin >= 0);
    assert(unknown_intervals.back().end <=
           curve_pair[1 - unknown_intervals.back().fixed_curve]->size() - 1);
#else
  CGAL_USE(begin);
  CGAL_USE(end);
  CGAL_USE(fixed_point);
  CGAL_USE(fixed_curve);
#endif
}

template <typename C>
void FrechetLight<C>::visAddConnection(CPoint begin, CPoint end,
                                       CPoint fixed_point, CurveID fixed_curve)
{
#ifdef VIS
    if (begin >= end) {
        return;
    }

    assert(begin >= 0);
    assert(end <= curve_pair[1 - fixed_curve]->size() - 1);

    connections.emplace_back(begin, end);
    connections.back().fixed = fixed_point;
    connections.back().fixed_curve = fixed_curve;

    assert(connections.back().begin >= 0);
    assert(connections.back().end <=
           curve_pair[1 - connections.back().fixed_curve]->size() - 1);
#else
  CGAL_USE(begin);
  CGAL_USE(end);
  CGAL_USE(fixed_point);
  CGAL_USE(fixed_curve);
#endif
}

template <typename C>
void FrechetLight<C>::visAddFreeNonReachable(CPoint begin, CPoint end,
                                          CPoint fixed_point,
                                          CurveID fixed_curve)
{
#ifdef VIS
    if (begin >= end) {
        return;
    }

    assert(begin >= 0);
    assert(end <= curve_pair[1 - fixed_curve]->size() - 1);

    free_non_reachable.emplace_back(begin, end);
    free_non_reachable.back().fixed = fixed_point;
    free_non_reachable.back().fixed_curve = fixed_curve;

    assert(free_non_reachable.back().begin >= 0);
    assert(free_non_reachable.back().end <=
           curve_pair[1 - free_non_reachable.back().fixed_curve]->size() - 1);
#else
  CGAL_USE(begin);
  CGAL_USE(end);
  CGAL_USE(fixed_point);
  CGAL_USE(fixed_curve);
#endif
}


template <typename C>
void FrechetLight<C>::merge(CIntervals& intervals,
                                CInterval const& new_interval) const
{
    if (new_interval.is_empty()) {
        return;
    }

    if (!intervals.empty() && new_interval.begin == intervals.back().end) {
        intervals.back().end = new_interval.end;
    } else {
        intervals.emplace_back(new_interval);
    }
}

// TODO: better distinction -- returned SimpleIntervals are only useful to
// determine full/not full

template <typename C>
QSimpleInterval<C> FrechetLight<C>::getTestFullQSimpleInterval(
    CPoint const& fixed, const Curve& fixed_curve, PointID min, PointID max,
    const Curve& curve) const
{
    QSimpleInterval qsimple;
    updateQSimpleInterval<TestFullMode, CPoint>(qsimple, fixed, fixed_curve,
                                                min, max, curve);
    return qsimple;
}

template <typename C>
template <class MODE, class IndexType>
inline bool FrechetLight<C>::updateQSimpleInterval(QSimpleInterval& qsimple,
                                                const IndexType& fixed,
                                                const Curve& fixed_curve,
                                                PointID min, PointID max,
                                                const Curve& curve) const
{
    if (qsimple.is_valid() || (qsimple.hasPartialInformation() &&
                               qsimple.getLastValidPoint() >= max)) {
        qsimple.validate();
        qsimple.clamp(CPoint{min, 0}, CPoint{max, 0});
        return false;  // parent information already valid
    }

    if (qsimple.hasPartialInformation()) {
        if (qsimple.getLastValidPoint() < min) {
            qsimple = QSimpleInterval();  // we know nothing about this part
        } else {
            qsimple.clamp(CPoint{min, 0}, CPoint{max, 0});
            if (!qsimple.getFreeInterval().is_empty()) {
                return false;  // even restricted to current interval, parent
                               // information still gives invalidity
            }
        }
    }

    if (!qsimple.hasPartialInformation()) {  // fresh computation, start with
                                             // heuristics
        // heuristic check:
        auto mid = (min + max) / 2;
        auto maxdist = (CGAL::max)(curve.curve_length(min, mid),
                                   curve.curve_length(mid, max));
        // TODO: return upper bound on error of fixed_point displacement to then
        // overestimate the mid_dist distance.
        auto fixed_point = fixed_curve.interpolate_at(fixed);
        auto mid_dist = Curve::distance(fixed_point, curve[mid], curve.traits());

        // heuristic tests avoiding sqrts
        if (certainly(distance > maxdist)){
          auto comp_dist1 = distance - maxdist;
          if (certainly(mid_dist <= comp_dist1)) { // Uncertain (A)
            qsimple.setFreeInterval(min, max);  // full
            qsimple.validate();

            return true;
          }
        }
        auto comp_dist2 = distance + maxdist;
        if (certainly(mid_dist > comp_dist2)) { // Uncertain (A)
          qsimple.setFreeInterval(max, min);  // empty
          qsimple.validate();

          return true;
        }
    }

    continueQSimpleSearch<MODE, IndexType>(qsimple, fixed, fixed_curve, min,
                                           max, curve);

    return true;
}

template <typename C>
template <class MODE, class IndexType>
inline void FrechetLight<C>::continueQSimpleSearch(QSimpleInterval& qsimple,
                                                const IndexType& fixed,
                                                const Curve& fixed_curve,
                                                PointID min, PointID max,
                                                const Curve& curve) const
{
    assert(!qsimple.hasPartialInformation() ||
           (certainly(qsimple.getLastValidPoint() >= min) &&
            certainly(qsimple.getLastValidPoint() <= max)));  // Uncertain (A)

    PointID start;
    if (qsimple.hasPartialInformation()) {
        start = qsimple.getLastValidPoint();
    } else {
        start = min;
    }

    // Logical assert: there should be no free point in [min, start):
    //  if start > min, then the free interval must be empty

    auto fixed_point = fixed_curve.interpolate_at(fixed);

    // this is an Uncertain
    auto current_free = Curve::distance(fixed_point, curve[start], curve.traits()) <= distance;  // Uncertain (A)
    // If we are not certain if the current node is free, we cannot do anything
    // with interval arithmetic, so we return with an invalid interval.
    if (is_indeterminate(current_free)) {
        qsimple.invalidate();
        return;
    }

    //
    // NOTE: From here on we know that current_free is certain, so we don't have to check. (A)
    //

    // Work directly on the simple_interval of the boundary
    // CPoint first_free = current_free ? CPoint{min,0} : CPoint{max+1,0};
    assert(!current_free || start == min);
    CPoint first_free = current_free
                            ? CPoint{min, 0}
                            : CPoint{max + 1, 0};  // qsimple.free.begin;

    assert(first_free > max ||
           Curve::distance(fixed_point, curve.interpolate_at(first_free), curve.traits()) <= distance); // Uncertain (A)

    std::size_t stepsize = static_cast<std::size_t>(max - start) / 2;
    if (stepsize < 1 || qsimple.hasPartialInformation()) {
        stepsize = 1;
    }
    for (PointID cur = start; cur < max;) {
        // heuristic steps:

        stepsize = std::min<std::size_t>(stepsize, max - cur);
        auto mid = cur + (stepsize + 1) / 2;
        auto maxdist = (CGAL::max)(curve.curve_length(cur, mid),
                                  curve.curve_length(mid, cur + stepsize));
        auto mid_dist = Curve::distance(fixed_point, curve[mid], curve.traits());

        if (current_free && certainly(distance > maxdist) &&               // Uncertain (A)
            certainly(mid_dist <= distance - maxdist)) {
            cur += stepsize;

            stepsize *= 2;
            continue;
        }
        if (!current_free && certainly(mid_dist > distance + maxdist)) { // Uncertain (A)
            cur += stepsize;

            stepsize *= 2;
            continue;
        }
        // if heuristic steps don't work, then reduce stepsize:
        if (stepsize > 1) {
            stepsize /= 2;

            continue;
        }

        // from here on stepsize == 1:
        // mid == end holds in this case
        auto end_dist = mid_dist;
        // if last and next point are both free:
        if (current_free && certainly(end_dist <= distance)) {  // Uncertain (A)
            ++cur;
            stepsize *= 2;

            continue;
        }

        // Now the interval is either not free or the last comparison was inderterminate.
        // In both cases we say that the interval is not full in TestFullMode as we cannot
        // call getInterval with IndexType that is not PointID.
        if (std::is_same<MODE, TestFullMode>::value) {
            return;  // Simple interval is not full -- can return in TestFullMode
        }
        // from here on, regular mode -> IndexType = PointID

        // otherwise we have to compute the intersection interval:
        // TODO: bad style: stripping down information added by getInterval
        CInterval temp_interval = FrechetLight::getInterval<IndexType>(
            fixed_curve, fixed, curve, cur);
        Interval interval = Interval(temp_interval.begin.getPoint() == cur
                                         ? temp_interval.begin.getFraction()
                                         : 1,
                                     temp_interval.end.getPoint() == cur
                                         ? temp_interval.end.getFraction()
                                         : 1);

        // do previous check for fullness again, but now it is an exact decision
        if (is_zero(interval.begin) && is_one(interval.end)) {  // Uncertain (A)
            assert(current_free);
            ++cur;
            stepsize *= 2;

            continue;
        }
        // from now on we know that the interval is not full

        if (interval.is_empty()) {
            ++cur;
            stepsize *= 2;

            continue;
        }
        // from here on the intersection interval is non-trivial:
        if (!current_free) {
            if (qsimple.is_valid()) {  // we encountered a second free interval
                                       // -> not qsimple
                qsimple.invalidate();
                // still store information in qsimple interval for future use
                qsimple.setLastValidPoint(cur);

                return;
            }

            // if two changes in current interval:
            if ((! is_zero(interval.begin)) && (! is_one(interval.end))) {
                qsimple.setFreeInterval(CPoint{cur, interval.begin},
                                        CPoint{cur, interval.end});
                qsimple.validate();
            } else {  // if one change in current interval:
                first_free = CPoint{cur, interval.begin};
                current_free = true;
            }

            ++cur;
        } else {
            assert(is_zero(interval.begin));  // current_free holds
            assert(! is_one(interval.end));  // end_dist <= distance does not hold, otherwise we
                        // would have stopped before

            qsimple.setFreeInterval(first_free, CPoint{cur, interval.end});
            qsimple.validate();
            current_free = false;
            ++cur;
        }
    }

    // If it is still invalid but we didn't return yet, then the free interval
    // ends at max Note: if no free points were found, first_free defaults to an
    // invalid value (>> max), making the interval empty
    if (!qsimple.is_valid()) {
        qsimple.setFreeInterval(first_free, CPoint{max, 0});
        qsimple.validate();
    }

    return;
}

template <typename C>
typename CIntervals<C>::iterator getIntervalContainingNumber(
    const typename CIntervals<C>::iterator& begin, const typename CIntervals<C>::iterator& end,
    CPoint<C> const& x)
{
    auto it = std::upper_bound(
        begin, end,
        CInterval<C>{x, CPoint<C>{(std::numeric_limits<typename C::PointID::IDType>::max)(), 0}});
    if (it != begin) {
        --it;
        if (it->begin <= x && it->end >= x) {
            return it;
        }
    }
    return end;
}

template <typename C>
typename CIntervals<C>::iterator getIntervalContainingNumber(
    const typename CIntervals<C>::iterator& begin, const typename CIntervals<C>::iterator& end,
    typename C::PointID x)
{
    auto it = std::upper_bound(
        begin, end,
        CInterval<C>{x, 0, (std::numeric_limits<typename C::PointID::IDType>::max)(), 0});
    if (it != begin) {
        --it;
        if (it->begin <= x && it->end >= x) {
            return it;
        }
    }
    return end;
}


template <typename C>
void FrechetLight<C>::getReachableIntervals(BoxData& data)
{
    ++num_boxes;

    auto const& box = data.box;
    auto const& inputs = data.inputs;
    CInterval const empty;

    assert(box.max1 > box.min1 && box.max2 > box.min2);
    // assert(outputs.id1.valid() || outputs.id2.valid());

    firstinterval1 = (inputs.begin1 != inputs.end1) ? &*inputs.begin1 : &empty;
    firstinterval2 = (inputs.begin2 != inputs.end2) ? &*inputs.begin2 : &empty;

    if (emptyInputsRule(data)) {
        return;
    }

    min1_frac = 0, min2_frac = 0;
    boxShrinkingRule(data);

    if (box.isCell()) {
        visAddCell(box);
        handleCellCase(data);
        return;
    } else {
        getQSimpleIntervals(data);
        calculateQSimple1(data);
        calculateQSimple2(data);

        if (out1_valid && out2_valid) {
            return;
        }
        if (boundaryPruningRule(data)) {
            return;
        }

        assert(box.max1 >= box.min1 + 2 || box.max2 >= box.min2 + 2);
        assert(box.max1 >= box.min1 && box.max2 >= box.min2);

        splitAndRecurse(data);
    }
}

template <typename C>
bool FrechetLight<C>::emptyInputsRule(BoxData& data)
{
    auto const& box = data.box;

    // empty input intervals -> empty output intervals
    // Note: if we are currently handling a cell then even if we are not
    // pruning, we have to return with empty outputs for the subsequent code to
    // work.
    if (pruning_level > 0 ||
        (box.min1 == box.max1 - 1 && box.min2 == box.max2 - 1)) {
        if (firstinterval2->is_empty() && firstinterval1->is_empty()) {
            if (data.outputs.id2.valid()) {
                visAddUnknown(CPoint(box.min2, 0), CPoint(box.max2, 0),
                              CPoint(box.max1, 0), 0);
            }
            if (data.outputs.id1.valid()) {
                visAddUnknown(CPoint(box.min1, 0), CPoint(box.max1, 0),
                              CPoint(box.max2, 0), 1);
            }
            return true;
        }
    }

    return false;
}

template <typename C>
void FrechetLight<C>::boxShrinkingRule(BoxData& data)
{
    auto& box = data.box;

    // "movement of input boundary" if one of the inputs is empty
    if (pruning_level > 1 && enable_box_shrinking) {
        if (firstinterval2->is_empty() && firstinterval1->begin > box.min1) {
            auto old_min1 = box.min1;

            min1_frac = firstinterval1->begin.getFraction();
            box.min1 = firstinterval1->begin.getPoint();
            assert(box.min1 <= box.max1);
            if (box.min1 == box.max1) {
                box.min1 = box.max1 - 1;
                min1_frac = 1;
            }

            if (box.min1 != old_min1) {
                visAddUnknown(CPoint(box.min2, 0), CPoint(box.max2, 0),
                              CPoint(box.min1, 0), 0);
                if (data.outputs.id1.valid()) {
                    visAddUnknown(CPoint(old_min1, 0), CPoint(box.min1, 0),
                                  CPoint(box.max2, 0), 1);
                }
            }
        } else if (firstinterval1->is_empty() &&
                   firstinterval2->begin > box.min2) {
            auto old_min2 = box.min2;

            min2_frac = firstinterval2->begin.getFraction();
            box.min2 = firstinterval2->begin.getPoint();
            assert(box.min2 <= box.max2);
            if (box.min2 == box.max2) {
                box.min2 = box.max2 - 1;
                min2_frac = 1;
            }

            if (box.min2 != old_min2) {
                visAddUnknown(CPoint(box.min1, 0), CPoint(box.max1, 0),
                              CPoint(box.min2, 0), 1);
                if (data.outputs.id2.valid()) {
                    visAddUnknown(CPoint(old_min2, 0), CPoint(box.min2, 0),
                                  CPoint(box.max1, 0), 0);
                }
            }
        }
    }
}


template <typename C>
void FrechetLight<C>::handleCellCase(BoxData& data)
{
    auto const& curve1 = *curve_pair[0];
    auto const& curve2 = *curve_pair[1];
    auto const& box = data.box;

    if (data.outputs.id1.valid()) {
        CInterval output1 = getInterval(curve2, box.max2, curve1, box.min1);
        if (firstinterval2->is_empty()) {
            visAddFreeNonReachable(output1.begin,
                                   (std::min)(output1.end, firstinterval1->begin),
                                   {box.max2, 0}, 1);
            output1.begin.setFraction(
                                      (CGAL::max)(output1.begin.getFraction(),
                                                  firstinterval1->begin.getFraction()));
            certSetValues(output1, *firstinterval1, box.max2, 1);
        } else {
            certSetValues(output1, *firstinterval2, box.max2, 1);
        }
        merge(reachable_intervals_vec[data.outputs.id1], output1);
        visAddReachable(output1);
    }

    if (data.outputs.id2.valid()) {
        CInterval output2 = getInterval(curve1, box.max1, curve2, box.min2);
        if (firstinterval1->is_empty()) {
            visAddFreeNonReachable(output2.begin,
                                   (std::min)(output2.end, firstinterval2->begin),
                                   {box.max1, 0}, 0);
            output2.begin.setFraction(
                                      (CGAL::max)(output2.begin.getFraction(),
                                                  firstinterval2->begin.getFraction()));
            certSetValues(output2, *firstinterval2, box.max1, 0);
        } else {
            certSetValues(output2, *firstinterval1, box.max1, 0);
        }
        merge(reachable_intervals_vec[data.outputs.id2], output2);
        visAddReachable(output2);
    }
}


template <typename C>
void FrechetLight<C>::getQSimpleIntervals(BoxData& data)
{
    auto const& curve1 = *curve_pair[0];
    auto const& curve2 = *curve_pair[1];
    auto const& box = data.box;

    qsimple1.invalidate();
    qsimple2.invalidate();

    // Get qsimple intervals. Different cases depending on what has been
    // calculated yet.
    if (data.outputs.id1.valid()) {
        if (data.qsimple_outputs.id1.valid()) {
            qsimple1 = qsimple_intervals[data.qsimple_outputs.id1];
        } else {
            qsimple1 = QSimpleInterval();
        }
        bool changed = updateQSimpleInterval(qsimple1, box.max2, curve2,
                                             box.min1, box.max1, curve1);
        if (changed) {
            qsimple_intervals.push_back(qsimple1);
            data.qsimple_outputs.id1 = qsimple_intervals.size() - 1;
        }
    }
    if (data.outputs.id2.valid()) {
        if (data.qsimple_outputs.id2.valid()) {
            qsimple2 = qsimple_intervals[data.qsimple_outputs.id2];
        } else {
            qsimple2 = QSimpleInterval();
        }
        bool changed = updateQSimpleInterval(qsimple2, box.max1, curve1,
                                             box.min2, box.max2, curve2);
        if (changed) {
            qsimple_intervals.push_back(qsimple2);
            data.qsimple_outputs.id2 = qsimple_intervals.size() - 1;
        }
    }
}

template <typename C>
void FrechetLight<C>::calculateQSimple1(BoxData& data)
{
    auto const& curve1 = *curve_pair[0];
    auto const& curve2 = *curve_pair[1];
    auto const& box = data.box;

    out1_valid = false;
    out1.make_empty();

    // pruning rules depending on qsimple1
    if (qsimple1.is_valid()) {
        // output boundary is empty
        if (qsimple1.is_empty()) {
            // out1 is already empty due to initialization, so leave it.
            if (pruning_level > 2 && enable_empty_outputs) {
                out1_valid = true;
            }
        } else {
            CPoint x =
                (qsimple1.getFreeInterval().begin > CPoint{box.min1, min1_frac})
                    ? qsimple1.getFreeInterval().begin
                    : CPoint{box.min1, min1_frac};
            // check if beginning is reachable
            if (x == box.min1 && pruning_level > 3 && enable_propagation1) {
                typename CIntervals::iterator it =
                getIntervalContainingNumber<C>(
                    data.inputs.begin2, data.inputs.end2, box.max2);
                if (it != data.inputs.end2) {  //(box.min1, box.max2) is
                                               // reachable from interval *it
                    CInterval& parent = *it;
                    out1 = qsimple1.getFreeInterval();
                    out1_valid = true;
                    certSetValues(out1, parent, box.max2, 1);
                }
            }
            // check if nothing can be reachable
            if (x != box.min1 && x > qsimple1.getFreeInterval().end &&
                pruning_level > 4 && enable_propagation2) {
                // out1 is already empty due to initialization, so leave it.
                out1_valid = true;
            }
            // check if we can propagate reachability through the box to the
            // beginning of the free interval
            else if (x > box.min1 && pruning_level > 4 && enable_propagation2) {
                auto it = getIntervalContainingNumber(data.inputs.begin1,
                                                      data.inputs.end1, x);
                if (it !=
                    data.inputs
                        .end1) {  //(x, box.min2) is reachable from interval *it
                    auto interval = getTestFullQSimpleInterval(
                        x, curve1, box.min2, box.max2, curve2);
                    if (interval.is_valid()) {
                        CInterval& parent = *it;
                        out1 = qsimple1.getFreeInterval();
                        out1.begin = x;
                        out1_valid = true;
                        certSetValues(out1, parent, box.max2, 1);
                        visAddConnection({box.min2, 0}, {box.max2, 0}, x, 0);
                    }
                }
            }
        }
    }
    if (out1_valid) {
        merge(reachable_intervals_vec[data.outputs.id1], out1);
        visAddReachable(out1);
        if (!out1.is_empty()) {
            visAddFreeNonReachable(qsimple1.getFreeInterval().begin, out1.begin,
                                   {box.max2, 0}, 1);
        }
        data.outputs.id1.invalidate();
    }
    out1_valid = !data.outputs.id1.valid();
}

template <typename C>
void FrechetLight<C>::calculateQSimple2(BoxData& data)
{
    auto const& curve1 = *curve_pair[0];
    auto const& curve2 = *curve_pair[1];
    auto const& box = data.box;

    out2_valid = false;
    out2.make_empty();

    // pruning rules depending on qsimple2
    if (qsimple2.is_valid()) {
        // output boundary is empty
        if (qsimple2.is_empty()) {
            // out2 is already empty due to initialization, so leave it.
            if (pruning_level > 2 && enable_empty_outputs) {
                out2_valid = true;
            }
        } else {
            CPoint x =
                (qsimple2.getFreeInterval().begin > CPoint{box.min2, min2_frac})
                    ? qsimple2.getFreeInterval().begin
                    : CPoint{box.min2, min2_frac};
            // check if beginning is reachable
            if (x == box.min2 && pruning_level > 3 && enable_propagation1) {
                typename CIntervals::iterator it
                 = getIntervalContainingNumber<C>(
                     data.inputs.begin1, data.inputs.end1, box.max1);
                if (it != data.inputs.end1) {
                    CInterval& parent = *it;
                    out2 = qsimple2.getFreeInterval();
                    out2_valid = true;
                    certSetValues(out2, parent, box.max1, 0);
                }
            }
            // check if nothing can be reachable
            if (x != box.min2 && x > qsimple2.getFreeInterval().end &&
                pruning_level > 4 && enable_propagation2) {
                // out2 is already empty due to initialization, so leave it.
                out2_valid = true;
            }
            // check if we can propagate reachability through the box to the
            // beginning of the free interval
            else if (x > box.min2 && pruning_level > 4 && enable_propagation2) {
                auto it = getIntervalContainingNumber(data.inputs.begin2,
                                                      data.inputs.end2, x);
                if (it != data.inputs.end2) {
                    auto interval = getTestFullQSimpleInterval(
                        x, curve2, box.min1, box.max1, curve1);
                    if (interval.is_valid()) {
                        CInterval& parent = *it;
                        out2 = qsimple2.getFreeInterval();
                        out2.begin = x;
                        out2_valid = true;
                        certSetValues(out2, parent, box.max1, 0);
                        visAddConnection({box.min1, 0}, {box.max1, 0}, x, 1);
                    }
                }
            }
        }
    }
    if (out2_valid) {
        merge(reachable_intervals_vec[data.outputs.id2], out2);
        visAddReachable(out2);
        if (!out2.is_empty()) {
            visAddFreeNonReachable(qsimple2.getFreeInterval().begin, out2.begin,
                                   {box.max1, 0}, 0);
        }
        data.outputs.id2.invalidate();
    }
    out2_valid = !data.outputs.id2.valid();
}


template <typename C>
bool FrechetLight<C>::boundaryPruningRule(BoxData& data)
{
    auto const& curve1 = *curve_pair[0];
    auto const& curve2 = *curve_pair[1];
    auto const& box = data.box;

    // special cases for boxes which are at the boundary of the freespace
    // diagram
    if (pruning_level > 5 && enable_boundary_rule) {
        if (box.max1 == curve1.size() - 1 && out1_valid) {
            visAddUnknown({box.min2, 0}, {box.max2, 0}, {box.max1, 0}, 0);
            return true;
        }
        if (box.max2 == curve2.size() - 1 && out2_valid) {
            visAddUnknown({box.min1, 0}, {box.max1, 0}, {box.max2, 0}, 1);
            return true;
        }
    }

    return false;
}

template <typename C>
void FrechetLight<C>::splitAndRecurse(BoxData& data)
{
    auto const& box = data.box;

    if (box.max2 - box.min2 > box.max1 - box.min1) {  // horizontal split
        reachable_intervals_vec.emplace_back();
        CIntervalsID inputs1_middleID = reachable_intervals_vec.size() - 1;

        PointID split_position = (box.max2 + box.min2) / 2;
        assert(split_position > box.min2 && split_position < box.max2);

        auto bound = CInterval{split_position, 0,
          (std::numeric_limits<typename PointID::IDType>::max)(), 0};
        auto it = std::upper_bound(data.inputs.begin2, data.inputs.end2, bound);

        BoxData data_bottom{
            {box.min1, box.max1, box.min2, split_position},
            {data.inputs.begin1, data.inputs.end1, data.inputs.begin2, it},
            {inputs1_middleID, data.outputs.id2},
            {QSimpleID<C>(), data.qsimple_outputs.id2}};
        getReachableIntervals(data_bottom);

        if (it != data.inputs.begin2 && (it - 1)->end >= split_position) {
            --it;
        }
        CIntervals& inputs1_middle = reachable_intervals_vec[inputs1_middleID];

        BoxData data_top{{box.min1, box.max1, split_position, box.max2},
                         {inputs1_middle.begin(), inputs1_middle.end(), it,
                          data.inputs.end2},
                         {data.outputs.id1, data.outputs.id2},
                         {data.qsimple_outputs.id1, data.qsimple_outputs.id2}};
        getReachableIntervals(data_top);
    } else {  // vertical split
        reachable_intervals_vec.emplace_back();
        CIntervalsID inputs2_middleID = reachable_intervals_vec.size() - 1;

        PointID split_position = (box.max1 + box.min1) / 2;
        assert(split_position > box.min1 && split_position < box.max1);

        auto bound = CInterval{split_position, 0,
                               (std::numeric_limits<typename PointID::IDType>::max)(), 0};
        auto it = std::upper_bound(data.inputs.begin1, data.inputs.end1, bound);

        BoxData data_left{
            {box.min1, split_position, box.min2, box.max2},
            {data.inputs.begin1, it, data.inputs.begin2, data.inputs.end2},
            {data.outputs.id1, inputs2_middleID},
            {data.qsimple_outputs.id1, QSimpleID<C>()}};
        getReachableIntervals(data_left);

        if (it != data.inputs.begin1 && (it - 1)->end >= split_position) {
            --it;
        }
        CIntervals& inputs2_middle = reachable_intervals_vec[inputs2_middleID];

        BoxData data_right{
            {split_position, box.max1, box.min2, box.max2},
            {it, data.inputs.end1, inputs2_middle.begin(),
             inputs2_middle.end()},
            {data.outputs.id1, data.outputs.id2},
            {data.qsimple_outputs.id1, data.qsimple_outputs.id2}};
        getReachableIntervals(data_right);
    }
}

template <typename C>
CPoint<C> FrechetLight<C>::getLastReachablePoint(Curve const& curve1, PointID i,
                                                             Curve const& curve2) const
{
    Point const& point = curve1[i];
    PointID max = curve2.size() - 1;
    std::size_t stepsize = 1;
    for (PointID cur = 0; cur < max;) {
        // heuristic steps:
        stepsize = std::min<std::size_t>(stepsize, max - cur);

        auto mid = cur + (stepsize + 1) / 2;
        auto first_part = curve2.curve_length(cur + 1, mid);
        auto second_part = curve2.curve_length(mid, cur + stepsize);
        auto maxdist = (CGAL::max)(first_part, second_part);
        auto mid_dist = Curve::distance(point, curve2[mid], curve2.traits());


        if(certainly(distance >  maxdist) && (certainly(mid_dist <= distance - maxdist))) { // Uncertain (A)
            cur += stepsize;
            stepsize *= 2;
        }


        // if heuristic steps don't work, then reduce stepsize:
        else if (stepsize > 1) {
            stepsize /= 2;
        } else {
            auto interval = getInterval(curve1, i, curve2, cur);

            // interval is full, so continue
            if (is_one(interval.end.getFraction())) {
                assert(is_zero(interval.begin.getFraction()));
                ++cur;
                stepsize *= 2;
            }

            // the last reachable point is on the current interval
            return interval.end;
        }
    }
    return CPoint{max, 0};
}

/*
template <typename C>
void FrechetLight<C>::buildFreespaceDiagram(distance_t const& distance,
                                            Curve const& curve1,
                                            Curve const& curve2)
{
    this->curve_pair[0] = &curve1;
    this->curve_pair[1] = &curve2;
    this->distance = distance;

    clear();

    Box initial_box(0, curve1.size() - 1, 0, curve2.size() - 1);
    Inputs initial_inputs = computeInitialInputs();
    Outputs final_outputs = createFinalOutputs();

    initCertificate(initial_inputs);

    // this is the main computation of the decision problem
    computeOutputs(initial_box, initial_inputs, final_outputs);
}
*/

template <typename C>
bool FrechetLight<C>::lessThan(distance_t const& d, Curve const& curve1,
                               Curve const& curve2)
{
    this->curve_pair[0] = &curve1;
    this->curve_pair[1] = &curve2;
    this->distance = C::to_ift(d);

    // curves empty or start or end are already far
    if (curve1.empty() || curve2.empty()) {
        return false;
    }
    if (certainly(Curve::distance(curve1.front(), curve2.front(), curve1.traits()) > this->distance)) { // Uncertain (A)
        return false;
    }
    if (certainly(Curve::distance(curve1.back(), curve2.back(), curve1.traits()) > this->distance)) { // Uncertain (A)
        return false;
    }

    // cases where at least one curve has length 1
    if (curve1.size() == 1 && curve2.size() == 1) {
        return true;
    }
    if (curve1.size() == 1) {
        return isClose(curve1, 0, curve2);
    }
    if (curve2.size() == 1) {
        return isClose(curve2, 0, curve1);
    }

    clear();

    Box<C> initial_box(0, curve1.size() - 1, 0, curve2.size() - 1);
    Inputs initial_inputs = computeInitialInputs();
    Outputs final_outputs = createFinalOutputs();

    initCertificate(initial_inputs);

    // this is the main computation of the decision problem
    computeOutputs(initial_box, initial_inputs, final_outputs);

    return isTopRightReachable(final_outputs);
}

template <typename C>
bool FrechetLight<C>::lessThanWithFilters(distance_t const& d, Curve const& curve1,
                                          Curve const& curve2)
{
    this->curve_pair[0] = &curve1;
    this->curve_pair[1] = &curve2;
    this->distance = Curve::to_ift(d);

    assert(curve1.size());
    assert(curve2.size());

    if (certainly(Curve::distance(curve1[0], curve2[0], curve1.traits()) > this->distance) ||          // Uncertain (A)
        certainly(Curve::distance(curve1.back(), curve2.back(), curve1.traits()) > this->distance)) {
        return false;
    }
    if (curve1.size() == 1 && curve2.size() == 1) {
        return true;
    }

    Filter filter(curve1, curve2, this->distance);

    if (filter.bichromaticFarthestDistance()) {
        return true;
    }

    PointID pos1;
    PointID pos2;
    if (filter.adaptiveGreedy(pos1, pos2)) {
        return true;
    }
    if (filter.negative(pos1, pos2)) {
        return false;
    }
    if (filter.adaptiveSimultaneousGreedy()) {
        return true;
    }

    ++non_filtered;

    return lessThan(d, curve1, curve2);
}

template <typename C>
void FrechetLight<C>::computeOutputs(Box<C> const& initial_box,
                                         Inputs const& initial_inputs,
                                         Outputs& final_outputs)
{
    num_boxes = 0;

    BoxData box_data{initial_box, initial_inputs, final_outputs,
                     QSimpleOutputs()};
    getReachableIntervals(box_data);
}

template <typename C>
void FrechetLight<C>::visAddCell(Box<C> const& box)
{
#ifdef VIS
    cells.emplace_back(box.min1, box.min2);
#else
  CGAL_USE(box);
#endif
}

template <typename C>
bool FrechetLight<C>::isClose(Curve const& curve1, PointID i,
                                  Curve const& curve2) const
{
    return getLastReachablePoint(curve1, i, curve2) ==
           CPoint{PointID(curve2.size() - 1), 0};
}

template <typename C>
bool FrechetLight<C>::isTopRightReachable(Outputs const& outputs) const
{
    auto const& curve1 = *curve_pair[0];
    auto const& curve2 = *curve_pair[1];
    auto const& outputs1 = reachable_intervals_vec[outputs.id1];
    auto const& outputs2 = reachable_intervals_vec[outputs.id2];

    return (!outputs1.empty() &&
            (outputs1.back().end.getPoint() == curve1.size() - 1)) ||
           (!outputs2.empty() &&
            (outputs2.back().end.getPoint() == curve2.size() - 1));
}

template <typename C>
void FrechetLight<C>::initCertificate(Inputs const& initial_inputs)
{
    auto const& curve1 = *curve_pair[0];
    auto const& curve2 = *curve_pair[1];

    CInterval origin = CInterval(0, 0);
    certSetValues(origin, origin, 0, -1);
    certSetValues(*initial_inputs.begin1, origin, 0, 1);
    certSetValues(*initial_inputs.begin2, origin, 0, 0);

    visAddUnknown(initial_inputs.begin1->end, CPoint(curve1.size() - 1, 0),
                  {0, 0}, 1);
    visAddUnknown(initial_inputs.begin2->end, CPoint(curve2.size() - 1, 0),
                  {0, 0}, 0);
    visAddReachable(*initial_inputs.begin1);
    visAddReachable(*initial_inputs.begin2);
}

template <typename C>
auto FrechetLight<C>::createFinalOutputs() -> Outputs
{
    Outputs outputs;

    reachable_intervals_vec.emplace_back();
    outputs.id1 = reachable_intervals_vec.size() - 1;
    reachable_intervals_vec.emplace_back();
    outputs.id2 = reachable_intervals_vec.size() - 1;

    return outputs;
}

// this function assumes that the start points of the two curves are close
template <typename C>
auto FrechetLight<C>::computeInitialInputs() -> Inputs
{
    Inputs inputs;

    auto const& curve1 = *curve_pair[0];
    auto const& curve2 = *curve_pair[1];

    auto const first = CPoint(0, 0);

    auto last1 = getLastReachablePoint(curve2, 0, curve1);
    reachable_intervals_vec.emplace_back();
    reachable_intervals_vec.back().emplace_back(first, last1);
    inputs.begin1 = reachable_intervals_vec.back().begin();
    inputs.end1 = reachable_intervals_vec.back().end();

    auto last2 = getLastReachablePoint(curve1, 0, curve2);
    reachable_intervals_vec.emplace_back();
    reachable_intervals_vec.back().emplace_back(first, last2);
    inputs.begin2 = reachable_intervals_vec.back().begin();
    inputs.end2 = reachable_intervals_vec.back().end();

    return inputs;
}

template <typename C>
std::pair<double,double> FrechetLight<C>::calcDistance(Curve const& curve1, Curve const& curve2, double epsilon)
{
  //TODO: no interval here for split?
    double min = 0;
    double max = curve1.getUpperBoundDistance(curve2);

    while (max - min >= epsilon) {
        auto split = (max + min) / 2.;
        if (lessThan(distance_t(split), curve1, curve2)) {
            max = split;
        } else {
            min = split;
        }
    }

    return std::make_pair(min,max);
}

// This doesn't have to be called but is handy to make time measurements more
// consistent such that the clears in the lessThan call doen't have to do
// anything.

template <typename C>
void FrechetLight<C>::clear()
{
    reachable_intervals_vec.clear();
    qsimple_intervals.clear();

#ifdef VIS
    unknown_intervals.clear();
    connections.clear();
    free_non_reachable.clear();
    cells.clear();
    reachable_intervals.clear();
#endif
}

template <typename C>
bool FrechetLight<C>::isOnLowerRight(const CPosition& pt) const
{
    return pt[0] == curve_pair[0]->size() - 1 || pt[1] == 0;
}

template <typename C>
bool FrechetLight<C>::isOnUpperLeft(const CPosition& pt) const
{
    return pt[0] == 0 || pt[1] == curve_pair[1]->size() - 1;
}

// TODO by André: I am not sure whether this is currently fine using interval arithmetic, but this only affects the certificates, which are currently not active. This needs to be checked in case certificates are used in the future.
template <typename C>
Certificate<C>& FrechetLight<C>::computeCertificate()
{
    cert = Certificate();

    cert.setCurves(curve_pair[0], curve_pair[1]);
    cert.setDistance(distance);

    auto const& curve1 = *curve_pair[0];
    auto const& curve2 = *curve_pair[1];

    // TODO test handling of special cases!
    // special cases:
    if (certainly(Curve::distance(curve1.front(), curve2.front(), curve1.traits()) > distance)) { // Uncertain (A)
        cert.setAnswer(false);
        cert.addPoint({CPoint(0, 0), CPoint(0, 0)});
        cert.validate();
        return cert;
    }
    if (certainly(Curve::distance(curve1.back(), curve2.back(), curve1.traits()) > distance)) { // Uncertain (A)
        cert.setAnswer(false);
        cert.addPoint(
            {CPoint(curve1.size() - 1, 0), CPoint(curve2.size() - 1, 0)});
        cert.validate();
        return cert;
    }

    // cases where at least one curve has length 1
    if (curve1.size() == 1 && curve2.size() == 1) {
        cert.setAnswer(true);
        cert.addPoint({CPoint(0, 0), CPoint(0, 0)});
        cert.validate();
        return cert;
    }
    if (curve1.size() == 1) {
        auto last = getLastReachablePoint(curve1, 0, curve2);
        if (last == CPoint(curve2.size() - 1, 0)) {
            cert.setAnswer(true);
            cert.addPoint({CPoint(0, 0), CPoint(0, 0)});
            cert.addPoint({CPoint(0, 0), CPoint(curve2.size() - 1, 0)});
            cert.validate();
            return cert;
        } else {
            CInterval outer;
            (void)getInterval(curve1, (PointID)0, curve2, last.getPoint());
            CPoint safe_empty = outer.begin > CPoint(last.getPoint(), 0)
                                    ? outer.begin
                                    : outer.end;
            assert(safe_empty > CPoint(last.getPoint(), 0) ||
                   safe_empty < CPoint(last.getPoint() + 1, 0));
            cert.setAnswer(false);
            cert.addPoint({CPoint(0, 0), safe_empty});
            cert.validate();
            return cert;
        }
    }
    if (curve2.size() == 1) {
        auto last = getLastReachablePoint(curve2, 0, curve1);
        if (last == CPoint(curve1.size() - 1, 0)) {
            cert.setAnswer(true);
            cert.addPoint({CPoint(0, 0), CPoint(0, 0)});
            cert.addPoint({CPoint(curve1.size() - 1, 0), CPoint(0, 0)});
            cert.validate();
            return cert;
        } else {
            CInterval outer;
            (void)getInterval(curve2, (PointID)0, curve1, last.getPoint());
            CPoint safe_empty = outer.begin > CPoint(last.getPoint(), 0)
                                    ? outer.begin
                                    : outer.end;
            assert(safe_empty > CPoint(last.getPoint(), 0) ||
                   safe_empty < CPoint(last.getPoint() + 1, 0));
            cert.setAnswer(false);
            cert.addPoint({safe_empty, CPoint(0, 0)});
            cert.validate();
            return cert;
        }
    }

    // TODO: Check for case of a single point!

    CIntervals const& outputs1 = reachable_intervals_vec[2];
    CIntervals const& outputs2 = reachable_intervals_vec[3];

    bool answer = false;
    CInterval const* last_interval;
    if (outputs1.size() &&
        (outputs1.back().end.getPoint() == curve1.size() - 1)) {
        answer = true;
        last_interval = &outputs1.back();
    } else if (outputs2.size() &&
               (outputs2.back().end.getPoint() == curve2.size() - 1)) {
        answer = true;
        last_interval = &outputs2.back();
    }

    cert.setAnswer(answer);

    if (answer) {
        std::vector<CPosition> rev_traversal;
        CPosition cur_pos = {CPoint(curve1.size() - 1, 0),
                             CPoint(curve2.size() - 1, 0)};
        rev_traversal.push_back(cur_pos);
        CInterval const* interval = last_interval;
        while (cur_pos[0] > 0 || cur_pos[1] > 0) {
            CPosition next_pos = {CPoint(0, 0), CPoint(0, 0)};

            next_pos[interval->fixed_curve] = interval->fixed;
            next_pos[1 - interval->fixed_curve] =
                interval->end > cur_pos[1 - interval->fixed_curve]
                    ? cur_pos[1 - interval->fixed_curve]
                    : interval->end;
            assert(next_pos[0] <= cur_pos[0] && next_pos[1] <= cur_pos[1]);
            if (next_pos[0] != cur_pos[0] || next_pos[1] != cur_pos[1]) {
                rev_traversal.push_back(next_pos);
            }

            if (next_pos[1 - interval->fixed_curve] != interval->begin) {
                next_pos[1 - interval->fixed_curve] = interval->begin;
                rev_traversal.push_back(next_pos);
            }

            assert(next_pos[0] <= cur_pos[0] && next_pos[1] <= cur_pos[1]);
            cur_pos = next_pos;
            interval = interval->reach_parent;
        }

        cert.validate();

        for (int t = rev_traversal.size() - 1; t >= 0; t--) {
            cert.addPoint(rev_traversal[t]);
        }
    }

    return cert;
}

template <typename C>
auto FrechetLight<C>::getCurvePair() const -> CurvePair { return curve_pair; }

template <typename C>
void FrechetLight<C>::setPruningLevel(int pruning_level)
{
    this->pruning_level = pruning_level;
}


template <typename C>
void FrechetLight<C>::setRules(std::array<bool, 5> const& enable)
{
    enable_box_shrinking = enable[0];
    enable_empty_outputs = enable[1];
    enable_propagation1 = enable[2];
    enable_propagation2 = enable[3];
    enable_boundary_rule = enable[4];
}

template <typename C>
std::size_t FrechetLight<C>::getNumberOfBoxes() const { return num_boxes; }

} // namespace internal
} // namespace Frechet_distance
} // namespace CGAL
