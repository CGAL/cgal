/*
Copyright (C) 2026  TU Eindhoven

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once

#include "graph_curve_traits_2.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>

namespace CGAL::internal {
template <class K> struct Graph_segment_curve_traits_2 {
	using Kernel = K;
	using Curve_representation_2 = std::monostate;
	using Curve_2 = CGAL::Segment_2<K>;
	using Point_2 = CGAL::Point_2<K>;

	static Curve_representation_2 representation(const Curve_2& curve) {
		return std::monostate();
	}
	static Curve_2 curve(const Point_2& start, const Point_2& end, const Curve_representation_2& rep) {
		return Curve_2(start, end);
	}
	static Curve_2 curve_reversed(const Point_2& start, const Point_2& end,
	                              const Curve_representation_2& rep) {
		return Curve_2(start, end);
	}
	static CGAL::Bbox_2 bbox(const Point_2& start, const Point_2& end,
	                         const Curve_representation_2& rep) {
		return curve(start, end, rep).bbox();
	}

	static CGAL::Direction_2<Kernel> source_direction(const Point_2& start, const Point_2& end,
		const Curve_representation_2& rep) {
		return CGAL::Direction_2<Kernel>(end - start);
	}

	static CGAL::Direction_2<Kernel> target_direction(const Point_2& start, const Point_2& end,
	                                                  const Curve_representation_2& rep) {
		return CGAL::Direction_2<Kernel>(start - end);
	}

	static void reverse_representation(const Point_2& start, const Point_2& end,
	                                   Curve_representation_2& rep) {
		// skip
	}
	static void move_start(const Point_2& start, const Point_2& end,
	                                   Curve_representation_2& rep) {
		// skip
	}
	static void move_end(const Point_2& start, const Point_2& end, Curve_representation_2& rep) {
		// skip
	}
	static void transform(const CGAL::Aff_transformation_2<Kernel>& trans, const Point_2& start,
	                      const Point_2& end, Curve_representation_2& rep) {
		// skip
	}
};
static_assert(GraphCurveTraits_2<Graph_segment_curve_traits_2<CGAL::Epick>>);
static_assert(GraphCurveTraits_2<Graph_segment_curve_traits_2<CGAL::Epeck>>);
} // namespace CGAL::internal
