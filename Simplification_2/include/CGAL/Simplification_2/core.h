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

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_with_history_2.h>
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

#include <numbers>

namespace CGAL::internal {

/// CGAL kernel for exact constructions (uses an exact number type).
using Exact = CGAL::Exact_predicates_exact_constructions_kernel;
/// CGAL kernel for inexact constructions.
using Inexact = CGAL::Exact_predicates_inexact_constructions_kernel;

/// The number type used for coordinates.
template <class K> using Number = typename K::FT;

/// A point in the plane. See \ref CGAL::Point_2.
template <class K> using Point = CGAL::Point_2<K>;
/// A vector in the plane. See \ref CGAL::Vector_2.
template <class K> using Vector = CGAL::Vector_2<K>;
/// A direction in the plane. See \ref CGAL::Direction_2.
template <class K> using Direction = CGAL::Direction_2<K>;
/// A circle in the plane. See \ref CGAL::Circle_2.
template <class K> using Circle = CGAL::Circle_2<K>;
/// A line in the plane. See \ref CGAL::Line_2.
template <class K> using Line = CGAL::Line_2<K>;
/// A line segment in the plane. See \ref CGAL::Segment_2.
template <class K> using Segment = CGAL::Segment_2<K>;
/// A ray in the plane. See \ref CGAL::Ray_2.
template <class K> using Ray = CGAL::Ray_2<K>;
/// An axis-aligned rectangle in the plane. See \ref CGAL::Iso_rectangle_2.
template <class K> using Rectangle = CGAL::Iso_rectangle_2<K>;
/// A triangle in the plane. See \ref CGAL::Triangle_2.
template <class K> using Triangle = CGAL::Triangle_2<K>;

/// A polygon in the plane. See \ref CGAL::Polygon_2.
template <class K, class Container = std::vector<Point<K>>> using Polygon = CGAL::Polygon_2<K, Container>;
/// A polyline in the plane. See \ref CGAL::Polyline_2.
template <class K> using Polyline = CGAL::Polyline_2<K>;
/// A polygon with holes in the plane. See \ref CGAL::Polygon_2.
template <class K, class Container = std::vector<Point<K>>> using PolygonWithHoles = CGAL::Polygon_with_holes_2<K, Container>;
/// A point set with polygonal boundaries. See \ref CGAL::Polygon_set_2.
template <class K, class Container = std::vector<Point<K>>, class DCEL = CGAL::Gps_default_dcel<CGAL::Gps_segment_traits_2<K, Container>>> using PolygonSet = CGAL::Polygon_set_2<K, Container, DCEL>;

/// Axis-aligned bounding box with inexact coordinates. See \ref
/// CGAL::Bbox_2.
using Box = CGAL::Bbox_2;

/// An epsilon value.
/**
 * Obviously this should be used only for non-exact computation. For this
 * reason we provide only a definition for `Number<Inexact>` and not for
 * `Number<Exact>`.
 */
constexpr const Number<Inexact> M_EPSILON = 0.0000001;

/// Returns the area of a polygon with holes.
/// TODO: move to more logical place
template <class K> Number<K> area(PolygonWithHoles<K> polygon) {
	Number<K> a = polygon.outer_boundary().area();
	for (const auto& h : polygon.holes())
		a -= h.area();
	return a; // = outer area minus area of each hole
}

/// An RGB color. Used for storing the color of elements to be drawn.
struct Color {
	/// Red component (integer 0-255).
	int r;
	/// Green component (integer 0-255).
	int g;
	/// Blue component (integer 0-255).
	int b;

	/// Returns a new color that is darker (\f$0 \le f < 1\f$) or lighter
	/// (\f$1 < f \le 2\f$).
	Color shaded(double f) const;
	/// Constructs the color black.
	Color();
	/// Constructs a color.
	Color(int r, int g, int b);
	/// Constructs a color from a single integer (useful combined with hexadecimal literals, e.g. 0xFFFFFF).
	Color(int rgb);
};

/// \f$2 \pi\f$, defined here for convenience.
constexpr Number<Inexact> two_pi = std::numbers::pi * 2;

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
