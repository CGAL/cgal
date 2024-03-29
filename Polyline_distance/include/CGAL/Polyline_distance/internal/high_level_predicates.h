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
#include <CGAL/Exact_rational.h>
#include <CGAL/Sqrt_extension/Sqrt_extension_type.h>
#include <CGAL/number_utils.h>

#include <iterator>

namespace CGAL {
namespace Polyline_distance {
namespace internal {
namespace HLPred
{

/*!
 * \ingroup PkgPolylineDistanceFunctions
 * computes
*/
inline
Interval intersection_interval(Point const & circle_center, distance_t radius,
                               Point const& line_start, Point const& line_end)
{
    Interval I;

    try {
        std::pair<RealType, RealType> II;
        // if not empty
        if (approximate_reals(circle_center, radius, line_start, line_end,
                              II)) {
            I = Interval(II.first, II.second);
        }
    } catch (const CGAL::Uncertain_conversion_exception& e) {
        std::pair<RealType, RealType> II;
        // if not empty
        if (exact_reals(circle_center, radius, line_start, line_end, II)) {
            I = Interval(II.first, II.second);
        }
    }

    return I;
}

}  // namespace HLPred

} // namespace internal
} // namespace Polyline_distance
} // namespace CGAL
