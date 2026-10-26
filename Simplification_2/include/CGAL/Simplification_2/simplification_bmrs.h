#pragma once

#include <CGAL/Simplification_2/graph_simplification/edge_moves.h>
#include <CGAL/Simplification_2/graph_simplification/orientation_restriction.h>
#include "straight_simplification_algorithm.h"

#include <CGAL/Direction_2.h>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace CGAL {
template <class K>
struct Simplification_BMRS_options {
	/// Points of different geometries within L-infinity distance prec are merged into one vertex.
	typename K::FT prec = 0;
	/// If not empty, the orientations of all edges are restricted to these orientations before simplifying,
	/// by replacing edges with other orientations by staircases. As the moves of the algorithm keep the orientations of edges,
	/// the simplification then only has edges with these orientations.
	/// An orientation includes its opposite direction: a direction d and its opposite -d give the same orientation.
	/// At least two distinct orientations are needed. See regular_orientations() for evenly spaced orientations.
	std::vector<Direction_2<K>> orientations;
	/// For restricting orientations: edges longer than lambda times the longest edge are first subdivided into edges of at most that length.
	double lambda = 1;
	/// For restricting orientations: the tolerance used for the staircases.
	double eps = 0.1;
};

/// Returns count orientations at equal angles, starting at the angle initial_angle (in radians) with the positive x-axis.
/// For example, count 2 gives rectilinear and count 4 octilinear orientations.
template <class K>
std::vector<Direction_2<K>> regular_orientations(int count, double initial_angle = 0) {
	std::vector<Direction_2<K>> orientations;
	for (int i = 0; i < count; ++i) {
		double angle = initial_angle + i * std::numbers::pi / count;
		orientations.emplace_back(std::cos(angle), std::sin(angle));
	}
	return orientations;
}

namespace detail {
template<class Kernel, bool History>
class BMRS_algorithm {
public:
	using Graph = CGAL::internal::EdgeMovesGraph<Kernel, History>;

private:
	static constexpr int DEPTH = 10;
	/// Quadtree child cells overlap their siblings by this fraction of the cell size
	static constexpr double FUZZ = 0.05;

	using BMRS = CGAL::internal::BuchinEtAl<Graph>;
	using VertexTree = BMRS::VertexTree;
	using EdgeTree = BMRS::EdgeTree;
	std::unique_ptr<VertexTree> m_vt;
	std::unique_ptr<EdgeTree> m_et;
	BMRS m_bmrs;

public:
	BMRS_algorithm(Graph& graph, std::unique_ptr<VertexTree> vt)
	    : m_vt(std::move(vt)), m_et(std::make_unique<EdgeTree>(m_vt->root_box(), DEPTH, FUZZ)), m_bmrs(graph, *m_et, *m_vt) {
		// Merges collinear degree-2 vertices, rebuilds the vertex tree to contain just the isolated vertices (points),
		// which are obstacles during simplification, and fills the edge tree.
		m_bmrs.initialize(true, true);
		if constexpr (History) {
			// The merges are part of preparing the input, so they should not be undone.
			graph.history().clear();
		}
	}

	bool run(std::optional<std::function<bool(int, typename Kernel::FT)>> stop = std::nullopt) {
		return m_bmrs.run(stop);
	}
};

/// Converts the orientations to unit vectors with angles in [0, pi), sorted by angle and without duplicates.
template <class K>
std::vector<CGAL::internal::Vector<CGAL::internal::Inexact>> orientation_vectors(const std::vector<Direction_2<K>>& orientations) {
	std::vector<double> angles;
	for (const auto& d : orientations) {
		double angle = std::atan2(CGAL::to_double(d.dy()), CGAL::to_double(d.dx()));
		angle = std::fmod(angle + 2 * std::numbers::pi, std::numbers::pi);
		if (angle >= std::numbers::pi - 1e-12) {
			angle = 0;
		}
		angles.push_back(angle);
	}
	std::sort(angles.begin(), angles.end());
	angles.erase(std::unique(angles.begin(), angles.end(), [](double a, double b) { return b - a < 1e-12; }), angles.end());
	if (angles.size() < 2) {
		throw std::invalid_argument("Restricting orientations requires at least two distinct orientations.");
	}
	std::vector<CGAL::internal::Vector<CGAL::internal::Inexact>> vectors;
	for (double angle : angles) {
		vectors.emplace_back(std::cos(angle), std::sin(angle));
	}
	return vectors;
}
}

/// The edge-move simplification algorithm for straight-edge geometries by Buchin, Meulemans, Van Renssen, and Speckmann.
/// The algorithm iteratively moves edges such that area is preserved, choosing the moves with the smallest symmetric difference.
/// Moves keep the orientations of edges, so restricting the orientations of the input (see Simplification_BMRS_options::orientations)
/// gives a schematization with only those orientations.
/// @tparam History whether to store the simplification operations for efficient backtracking.
/// @tparam Kernel 2D linear kernel.
/// @tparam Data data stored with each geometry, or std::monostate for none.
/// @tparam PolygonContainer the container of the polygon vertices.
template<bool History = false, class Kernel = CGAL::Epick, class Data = std::monostate, class PolygonContainer = std::vector<CGAL::Point_2<Kernel>>>
class Simplification_BMRS : public internal::Straight_simplification_algorithm<Kernel, detail::BMRS_algorithm<Kernel, History>, Data, PolygonContainer> {
	using Base = internal::Straight_simplification_algorithm<Kernel, detail::BMRS_algorithm<Kernel, History>, Data, PolygonContainer>;
	using Graph = typename detail::BMRS_algorithm<Kernel, History>::Graph;

	static std::function<void(Graph&)> restriction(const Simplification_BMRS_options<Kernel>& options) {
		if (options.orientations.empty()) {
			return {};
		}
		auto vectors = detail::orientation_vectors(options.orientations);
		return [vectors, lambda = options.lambda, eps = options.eps](Graph& graph) {
			CGAL::internal::restrict_orientations(graph, vectors, lambda, eps);
		};
	}

public:
	/// Inserts the geometries of the collection, restricts their orientations if options.orientations is not empty, and prepares the simplification.
	/// Throws std::invalid_argument if fewer than two distinct orientations are given.
	Simplification_BMRS(const typename Base::Collection& collection, const Simplification_BMRS_options<Kernel>& options = {})
	    : Base(collection, options.prec, restriction(options)) {}
};

template<class K, class Data, class PolygonContainer>
Simplification_BMRS(const CGAL::StraightGeometryCollection<K, Data, PolygonContainer>&)
    -> Simplification_BMRS<false, K, Data, PolygonContainer>;

template<class K, class Data, class PolygonContainer>
Simplification_BMRS(const CGAL::StraightGeometryCollection<K, Data, PolygonContainer>&, const Simplification_BMRS_options<K>&)
    -> Simplification_BMRS<false, K, Data, PolygonContainer>;

// ---- Simplifying a collection ---- //
/// Simplifies all geometries of the collection together, until the number of edges is at most the target complexity.
/// Boundaries and vertices shared by several geometries are simplified once, so they stay shared, and shared edges count once.
/// The result has the geometries in the same order and with the same data as the input, so handles of the input remain valid.
template<class K, class Data, class PolygonContainer>
CGAL::StraightGeometryCollection<K, Data, PolygonContainer> simplify_BMRS(const CGAL::StraightGeometryCollection<K, Data, PolygonContainer>& collection, size_t complexity, const Simplification_BMRS_options<K>& options = {}) {
	return internal::simplify_collection<Simplification_BMRS<false, K, Data, PolygonContainer>>(collection, complexity, options);
}

// ---- Simplifying a single geometry ---- //
template<class K, class PolygonContainer>
CGAL::Polygon_2<K, PolygonContainer> simplify_BMRS(const CGAL::Polygon_2<K, PolygonContainer>& polygon, size_t complexity, const Simplification_BMRS_options<K>& options = {}) {
	return internal::simplify_geometry<Simplification_BMRS<false, K, std::monostate, PolygonContainer>>(polygon, complexity, options);
}

template<class K, class PolygonContainer>
CGAL::Polygon_with_holes_2<K, PolygonContainer> simplify_BMRS(const CGAL::Polygon_with_holes_2<K, PolygonContainer>& polygon, size_t complexity, const Simplification_BMRS_options<K>& options = {}) {
	return internal::simplify_geometry<Simplification_BMRS<false, K, std::monostate, PolygonContainer>>(polygon, complexity, options);
}

template<class K, class PolygonContainer>
CGAL::Multipolygon_with_holes_2<K, PolygonContainer> simplify_BMRS(const CGAL::Multipolygon_with_holes_2<K, PolygonContainer>& multiPolygon, size_t complexity, const Simplification_BMRS_options<K>& options = {}) {
	CGAL::Multipolygon_with_holes_2<K, PolygonContainer> input;
	for (const auto& polygon : multiPolygon.polygons_with_holes()) {
		input.add_polygon_with_holes(polygon);
	}
	return internal::simplify_geometry<Simplification_BMRS<false, K, std::monostate, PolygonContainer>>(input, complexity, options);
}

template<class K>
CGAL::Polyline_2<K> simplify_BMRS(const CGAL::Polyline_2<K>& polyline, size_t complexity, const Simplification_BMRS_options<K>& options = {}) {
	return internal::simplify_geometry<Simplification_BMRS<false, K>>(polyline, complexity, options);
}

template<class K>
CGAL::MultiPolyline<K> simplify_BMRS(const CGAL::MultiPolyline<K>& multiPolyline, size_t complexity, const Simplification_BMRS_options<K>& options = {}) {
	return internal::simplify_geometry<Simplification_BMRS<false, K>>(multiPolyline, complexity, options);
}
}
