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

#include <CGAL/Bbox_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Line_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Simplification_2/geometries/Polyline_2.h>

namespace CGAL::internal {
// This file defines some aliases we use in the Simplification_2 package and functions to convert between Epick and Epeck.

using Exact = CGAL::Exact_predicates_exact_constructions_kernel;
using Inexact = CGAL::Exact_predicates_inexact_constructions_kernel;

template <class K> using Number = typename K::FT;
template <class K> using Point = CGAL::Point_2<K>;
template <class K> using Vector = CGAL::Vector_2<K>;
template <class K> using Direction = CGAL::Direction_2<K>;
template <class K> using Circle = CGAL::Circle_2<K>;
template <class K> using Line = CGAL::Line_2<K>;
template <class K> using Segment = CGAL::Segment_2<K>;
template <class K> using Ray = CGAL::Ray_2<K>;
template <class K> using Rectangle = CGAL::Iso_rectangle_2<K>;
template <class K> using Triangle = CGAL::Triangle_2<K>;
template <class K, class Container = std::vector<Point<K>>> using Polygon = CGAL::Polygon_2<K, Container>;
template <class K> using Polyline = CGAL::Polyline_2<K>;
template <class K, class Container = std::vector<Point<K>>> using PolygonWithHoles = CGAL::Polygon_with_holes_2<K, Container>;
template <class K, class Container = std::vector<Point<K>>, class DCEL = CGAL::Gps_default_dcel<CGAL::Gps_segment_traits_2<K, Container>>> using PolygonSet = CGAL::Polygon_set_2<K, Container, DCEL>;
using Box = CGAL::Bbox_2;

/// An epsilon value.
constexpr const Number<Inexact> M_EPSILON = 0.0000001;

namespace detail {
constexpr const CGAL::Cartesian_converter<Exact, Inexact> exact_to_inexact;
constexpr const CGAL::Cartesian_converter<Inexact, Exact> inexact_to_exact;

template <typename K, template <typename> class Type>
concept IsCGALPrimitiveType =
    std::is_same<Type<K>, Point<K>>::value || std::is_same<Type<K>, Vector<K>>::value ||
    std::is_same<Type<K>, Line<K>>::value || std::is_same<Type<K>, Ray<K>>::value ||
    std::is_same<Type<K>, Segment<K>>::value || std::is_same<Type<K>, Rectangle<K>>::value ||
    std::is_same<Type<K>, Triangle<K>>::value || std::is_same<Type<K>, Circle<K>>::value;
} // namespace detail

template <typename KernelOut, typename KernelIn, template <typename> class Type>
requires detail::IsCGALPrimitiveType<KernelIn,Type> Type<KernelOut> convert_kernel(const Type<KernelIn>& v) {

	if constexpr (std::is_same<KernelIn, KernelOut>::value) {
		return v;
	} else if constexpr (std::is_same<KernelIn, Exact>::value &&
	                     std::is_same<KernelOut, Inexact>::value) {
		return detail::exact_to_inexact(v);
	} else if constexpr (std::is_same<KernelIn, Inexact>::value &&
	                     std::is_same<KernelOut, Exact>::value) {
		return detail::inexact_to_exact(v);
	} else {
		CGAL::Cartesian_converter<KernelIn, KernelOut> conv;
		return conv(v);
	}
}

template <typename KernelOut, typename KernelIn>
Number<KernelOut> convert_kernel(const Number<KernelIn>& v) {

	if constexpr (std::is_same<KernelIn, KernelOut>::value) {
		return v;
	} else if constexpr (std::is_same<KernelIn, Exact>::value &&
	                     std::is_same<KernelOut, Inexact>::value) {
		return detail::exact_to_inexact(v);
	} else if constexpr (std::is_same<KernelIn, Inexact>::value &&
	                     std::is_same<KernelOut, Exact>::value) {
		return detail::inexact_to_exact(v);
	} else {
		CGAL::Cartesian_converter<KernelIn, KernelOut> conv;
		return conv(v);
	}
}

template <typename KernelOut>
Number<KernelOut> convert_kernel(const Number<Inexact>& v) {

	if constexpr (std::is_same<Inexact, KernelOut>::value) {
		return v;
	} else if constexpr (std::is_same<KernelOut, Exact>::value) {
		return detail::inexact_to_exact(v);
	} else {
		CGAL::Cartesian_converter<Inexact, KernelOut> conv;
		return conv(v);
	}
}

template <typename KernelOut> Number<KernelOut> convert_kernel(const Number<Exact>& v) {

	if constexpr (std::is_same<Exact, KernelOut>::value) {
		return v;
	} else if constexpr (std::is_same<KernelOut, Inexact>::value) {
		return detail::exact_to_inexact(v);
	} else {
		CGAL::Cartesian_converter<Exact, KernelOut> conv;
		return conv(v);
	}
}

template <typename KernelOut, typename KernelIn>
Polygon<KernelOut> convert_kernel(const Polygon<KernelIn>& v) {
	std::vector<Point<KernelOut>> points;
	for (Point<KernelIn> p : v.vertices()) {
		points.push_back(convert_kernel<KernelOut>(p));
	}
	return {points.begin(), points.end()};
}

template <typename KernelOut, typename KernelIn>
PolygonWithHoles<KernelOut> convert_kernel(const PolygonWithHoles<KernelIn>& v) {
	auto outer = convert_kernel<KernelOut>(v.outer_boundary());
	std::vector<Polygon<KernelOut>> holes;
	for (Polygon<KernelIn> h : v.holes()) {
		holes.push_back(convert_kernel<KernelOut>(h));
	}
	return {outer, holes.begin(), holes.end()};
}

template <typename KernelOut, typename KernelIn>
PolygonSet<KernelOut> convert_kernel(const PolygonSet<KernelIn>& v) {
	std::vector<PolygonWithHoles<KernelIn>> polygons;
	v.polygons_with_holes(std::back_inserter(polygons));
	PolygonSet<KernelOut> result;
	for (const PolygonWithHoles<KernelIn>& polygon : polygons) {
		result.insert(convert_kernel<KernelOut>(polygon));
	}
	return result;
}

template <typename KernelIn, template <typename> class Type>
Type<Exact> pretendExact(const Type<KernelIn>& v) {
	return convert_kernel<Exact>(v);
}

template <typename KernelIn, template <typename> class Type>
Type<Inexact> approximate(const Type<KernelIn>& v) {
	return convert_kernel<Inexact>(v);
}

template<typename KernelIn> 
Polygon<Exact> pretendExact(const Polygon<KernelIn>& v) {
	return convert_kernel<Exact>(v);
}

template <typename KernelIn> 
Polygon<Inexact> approximate(const Polygon<KernelIn>& v) {
	return convert_kernel<Inexact>(v);
}

template <typename KernelIn>
PolygonWithHoles<Exact> pretendExact(const PolygonWithHoles<KernelIn>& v) {
	return convert_kernel<Exact>(v);
}

template <typename KernelIn>
PolygonWithHoles<Inexact> approximate(const PolygonWithHoles<KernelIn>& v) {
	return convert_kernel<Inexact>(v);
}

template <typename KernelIn> 
PolygonSet<Exact> pretendExact(const PolygonSet<KernelIn>& v) {
	return convert_kernel<Exact>(v);
}

template <typename KernelIn> 
PolygonSet<Inexact> approximate(const PolygonSet<KernelIn>& v) {
	return convert_kernel<Inexact>(v);
}

template <typename KernelIn> 
Number<Exact> pretendExact(const Number<KernelIn>& v) {
	return convert_kernel<Exact>(v);
}

template <typename KernelIn> 
Number<Inexact> approximate(const Number<KernelIn>& v) {
	return convert_kernel<Inexact>(v);
}

inline Number<Exact> pretendExact(const Number<Inexact>& v) {
	return convert_kernel<Exact>(v);
}

inline Number<Exact> pretendExact(const Number<Exact>& v) {
	return convert_kernel<Exact>(v);
}

inline Number<Inexact> approximate(const Number<Inexact>& v) {
	return convert_kernel<Inexact>(v);
}

inline Number<Inexact> approximate(const Number<Exact>& v) {
	return convert_kernel<Inexact>(v);
}

} // namespace CGAL::internal
