// Copyright (c) 2019 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//                 Karl Bringmann <kbringma@mpi-inf.mpg.de>
//
// =============================================================================

#pragma once
#include <CGAL/license/Polyline_distance.h>
#include <vector>

#include <CGAL/Polyline_distance/internal/curve.h>
#include <CGAL/Polyline_distance/internal/geometry_basics.h>
#include <CGAL/Polyline_distance/internal/id.h>

namespace CGAL {
namespace Polyline_distance {
namespace internal {

//
// Box
//

/*!
 * \ingroup PkgPolylineDistanceFunctions
 * A class representing a
*/
template <typename K>
struct Box {
    typedef typename K::PointID PointID;
    PointID min1;
    PointID max1;
    PointID min2;
    PointID max2;

    Box(PointID min1, PointID max1, PointID min2, PointID max2)
        : min1(min1), max1(max1), min2(min2), max2(max2)
    {
    }

    bool isCell() const { return max1 - min1 == 1 && max2 - min2 == 1; }
};

template <typename K>
using Boxes = std::vector<Box<K>>;

//
// Inputs
//

/*!
 * \ingroup PkgPolylineDistanceFunctions
 * A class representing a
*/
template <typename K>
struct Inputs {
    typedef typename K::PointID PointID;
    typename CIntervals<K>::iterator begin1;
    typename CIntervals<K>::iterator end1;
    typename CIntervals<K>::iterator begin2;
    typename CIntervals<K>::iterator end2;

    bool haveDownInput() const { return begin1 != end1; }
    bool haveLeftInput() const { return begin2 != end2; }

    bool downContains(PointID point_id) const
    {
        for (auto it = begin1; it != end1; ++it) {
            if (it->begin <= point_id && it->end >= point_id) {
                return true;
            }
        }
        return false;
    }
    bool leftContains(PointID point_id) const
    {
        for (auto it = begin2; it != end2; ++it) {
            if (it->begin <= point_id && it->end >= point_id) {
                return true;
            }
        }
        return false;
    }
};

//
// Outputs
//

/*!
 * \ingroup PkgPolylineDistanceFunctions
 * A class representing a
*/

template <typename K>
struct Outputs {
    CIntervalsID<K> id1;
    CIntervalsID<K> id2;
};

//
// QSimpleInterval
//

/*!
 * \ingroup PkgPolylineDistanceFunctions
 * A class representing a
*/
template <typename K>
struct QSimpleInterval {
    typedef typename K::PointID PointID;
    using CPoint = CGAL::Polyline_distance::internal::CPoint<K>;
    using CInterval = CGAL::Polyline_distance::internal::CInterval<K>;

    QSimpleInterval() : valid(false) {}
    QSimpleInterval(CPoint const& begin, CPoint const& end)
        : valid(true), free(begin, end)
    {
    }

    void setFreeInterval(CPoint const& begin, CPoint const& end)
    {
        free = {begin, end};
        // valid = true;
    }
    void setFreeInterval(PointID begin, PointID end)
    {
        free = {begin, 0, end, 0};
        // valid = true;
    }
    void invalidateFreeInterval()
    {
        free.begin = CPoint{};
        free.end = CPoint{};
    }

    CInterval const& getFreeInterval() { return free; }

    void setLastValidPoint(PointID const& point)
    {
        last_valid_point = point;
        // valid = false;
    }
    PointID const& getLastValidPoint() const { return last_valid_point; }
    bool hasPartialInformation() const { return last_valid_point.valid(); }

    void validate() { valid = true; }
    void invalidate() { valid = false; }
    void clamp(CPoint const& min, CPoint const& max) { free.clamp(min, max); }

    bool is_empty() const { return free.is_empty(); }
    bool is_valid() const { return valid; }

private:
    bool valid;
    CInterval free;
    PointID last_valid_point;
};

template <typename K>
using QSimpleIntervals = std::vector<QSimpleInterval<K>>;

template <typename K>
using QSimpleID = ID<QSimpleInterval<K>>;

//
// QSimpleOutputs
//

/*!
 * \ingroup PkgPolylineDistanceFunctions
 * A class representing a
*/

template <typename K>
struct QSimpleOutputs {
    QSimpleID<K> id1;
    QSimpleID<K> id2;
};

//
// BoxData
//

/*!
 * \ingroup PkgPolylineDistanceFunctions
 * A class representing a
*/
template <typename K>
struct BoxData {
    Box<K> box;
    Inputs<K> inputs;
    Outputs<K> outputs;
    QSimpleOutputs<K> qsimple_outputs;
};

} // namespace internal
} // namespace Polyline_distance
} // namespace CGAL
