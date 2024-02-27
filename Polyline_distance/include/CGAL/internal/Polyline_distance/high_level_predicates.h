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

#include "low_level_predicates.h"
#include "predicate_types.h"

#include <CGAL/Exact_rational.h>
#include <CGAL/Sqrt_extension/Sqrt_extension_type.h>
#include <CGAL/number_utils.h>
#include <iterator>

namespace HLPred
{

// FIXME: Yes... bad style, but this is just temporary anyway.
using namespace LLPred;

bool dist_at_most(Point const& p, Point const& q, distance_t d)
{
        return compare_squared_distance(p, q, d*d) != LARGER;
}

Interval intersection_interval2(const Point& circle_center, distance_t radius, Point line_start, Point line_end, Interval* outer)
{
        using OutputType = std::pair<Circular_arc_point_2, unsigned>;

        auto circle_2 = Circle(circle_center, radius);
        auto line_arc_2 = LineArc(line_start, line_end);

        std::vector<OutputType> intersections;
        intersection(circle_2, line_arc_2, std::back_inserter(intersections), outer);

        std::vector<distance_t> ratios;
        for (auto const& intersection: intersections) {
                assert(line_start.x() != line_end.x() || line_start.y() != line_end.y());

                distance_t ratio;
                // TODO: replace by projection onto coordinate axis
                if (line_start.x() != line_end.x()) {
                        ratio = (intersection.first.x() - line_start.x())/(line_end.x() - line_start.x());
                }
                else {
                        ratio = (intersection.first.y() - line_start.y())/(line_end.y() - line_start.y());
                }
                ratio = std::max(0., std::min(1., ratio));
                ratios.push_back(ratio);
        }

        assert(intersections.size() <= 2);

        switch (intersections.size())
        {
        case 0:
                return Interval(); // empty interval
        case 1:
                if (intersections[0].second == 2) { return Interval(ratios[0], ratios[0]); }
                if (dist_at_most(line_start, circle_center, radius)) {
                        return Interval(0., ratios[0]);
                }
                if (dist_at_most(line_end, circle_center, radius)) {
                        return Interval(ratios[0], 1.);
                }
                assert(false);
        case 2:
                return Interval(std::min(ratios[0], ratios[1]), std::max(ratios[0], ratios[1]));
        }
}

// TODO: something is wrong here
Interval intersection_interval(const Point& circle_center, distance_t radius, Point line_start, Point line_end, Interval* outer)
{
        //////////////////////////
        // Just keeping the call to intersection to leave the assignment to "outer" as before
        // TODO: remove when removing all the "outer" stuff
        {
                // auto I = intersection_interval2(circle_center, radius, line_start, line_end, outer);
                // std::cout << "old: " << CGAL::to_double(I.begin) << " " << CGAL::to_double(I.end) << std::endl;
        }
        //////////////////////////

    // let a := v.x()^2 + v.y()^2,
        // let b := line_start.x() * v.x() + line_start.y() * v.y(),
        // let c := line_start.x()^2 + line_start.y()^2 - radius^2
    // <=> lambda^2 * a + lambda * 2 b + c = 0
    // <=> lambda^2 + (2 b / a) * lambda + (c / a) = 0
    // <=> lambda1/2 = - (b / a) +/- sqrt((b / a)^2 - c / a)
        Interval I;
#if 1
        try {
          std::pair<RealType,RealType> II;
          if(approximate_reals(circle_center, radius, line_start, line_end, II)){
            I = Interval(II.first, II.second);
          }else{
            std::cout << "empty" << std::endl;
          }

        } catch(const CGAL::Uncertain_conversion_exception& e){
          std::cout << "problem with interval arithmetic" << std::endl;
          std::pair<RealType,RealType> II;
          if(exact_reals(circle_center, radius, line_start, line_end, II)){
            I = Interval(II.first, II.second);
          }else{
            std::cout << "empty" << std::endl;
          }
        }

#else
        // TODO: use actual dimension here instead of 2!
        Rational a,b,c = 0;
        for (auto i = 0; i < 2; ++i) {
                Rational diff = line_end[i] - line_start[i];
                a += CGAL::square(diff);
                b += (line_start[i]-circle_center[i])*diff;
                c += CGAL::square(line_start[i]-circle_center[i]);
        }
        c -= CGAL::square(Rational(radius));


        if (CGAL::square(b/a) - c/a >= 0.) {
                // std::cout << CGAL::to_double(a) << " " << CGAL::to_double(b) << " " << CGAL::to_double(c) << std::endl;
                // std::cout << CGAL::to_double(CGAL::square(b/a) - c/a) << std::endl;

                // auto start = CGAL::max(RealType(0.), RealType(-b/a - CGAL::sqrt(CGAL::square(b/a) - c/a)));
                // auto end = CGAL::min(RealType(1.), RealType(-b/a + CGAL::sqrt(CGAL::square(b/a) - c/a)));
                auto start = CGAL::max(RealType(0.), RealType(-b/a, -1, CGAL::square(b/a) - c/a));
                auto end = CGAL::min(RealType(1.), RealType(-b/a, 1, CGAL::square(b/a) - c/a));

                if (start <= RealType(1.) && end >= RealType(0.)) {
                        I = Interval(start, end);
                }
        }
        else {
                // std::cout << "empty" << std::endl;
        }
#endif
        // std::cout << "new: " << CGAL::to_double(I.begin) << " " << CGAL::to_double(I.end) << std::endl << std::endl;
        if (outer != nullptr) { *outer = I; }

        return I;
}

} // end HLPred namespace
