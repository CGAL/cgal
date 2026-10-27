#pragma once

#include <CGAL/Simplification_2/graph_simplification/edge_collapse.h>
#include "straight_simplification_algorithm.h"

namespace CGAL {
namespace internal {
template<class Kernel, bool History>
class KSBB_algorithm {
public:
	using Graph = CGAL::internal::EdgeCollapseGraph<Kernel, History>;

private:
	static constexpr int quadtree_depth = 10;
	/// Quadtree child cells overlap their siblings by this fraction of the cell size
	static constexpr double quadtree_fuzz = 0.05;

	using KSBB = CGAL::internal::KronenfeldEtAl<Graph>;
	using Vertex_tree = KSBB::VertexTree;
	using Edge_tree = KSBB::EdgeTree;
	std::unique_ptr<Vertex_tree> m_vt;
	std::unique_ptr<Edge_tree> m_et;
	KSBB m_ksbb;

public:
	KSBB_algorithm(Graph& graph, std::unique_ptr<Vertex_tree> vt)
	    : m_vt(std::move(vt)), m_et(std::make_unique<Edge_tree>(m_vt->root_box(), quadtree_depth, quadtree_fuzz)), m_ksbb(graph, *m_et, *m_vt) {
		// The edge collapse only uses isolated vertices (points) as obstacles,
		// so the vertex tree is rebuilt to contain just those, and the edge tree is filled.
		m_ksbb.initialize(true, true);
	}

	bool run(std::optional<std::function<bool(int, typename Kernel::FT)>> stop = std::nullopt) {
		return m_ksbb.run(stop);
	}
};
}

/// The edge-collapse simplification algorithm for straight-edge geometries by Kronenfeld, Stanislawski, Buttenfield, and Brockmeyer.
/// The algorithm iteratively replaces three consecutive edges by two, effectively collapsing an edge to a point.
/// The collapse is done such that area is preserved exactly (up to floating point errors).
/// The algorithm minimizes the area of symmetric difference of each step with the previous step.
/// Only edges whose endpoints both have degree two are collapsed.
/// @tparam History whether to store the simplification operations for efficient backtracking.
/// @tparam Kernel 2D linear kernel.
/// @tparam Data data stored with each geometry, or std::monostate for none.
/// @tparam PolygonContainer the container of the polygon vertices.
template<bool History = false, class Kernel = CGAL::Epick, class Data = std::monostate, class PolygonContainer = std::vector<CGAL::Point_2<Kernel>>>
class Simplification_KSBB : public internal::Straight_simplification_algorithm<Kernel, internal::KSBB_algorithm<Kernel, History>, Data, PolygonContainer> {
	// all methods are inherited from the base class
	using Base = internal::Straight_simplification_algorithm<Kernel, internal::KSBB_algorithm<Kernel, History>, Data, PolygonContainer>;
public:
	using Base::Base;
};

template<class K, class Data, class PolygonContainer>
Simplification_KSBB(const CGAL::Straight_geometry_collection_2<K, Data, PolygonContainer>&)
    -> Simplification_KSBB<false, K, Data, PolygonContainer>;

template<class K, class Data, class PolygonContainer>
Simplification_KSBB(const CGAL::Straight_geometry_collection_2<K, Data, PolygonContainer>&, const Simplification_options<K>&)
    -> Simplification_KSBB<false, K, Data, PolygonContainer>;

// ---- Simplifying a collection ---- //
/// Simplifies all geometries of the collection together, until the number of edges is at most the target complexity.
/// Boundaries and vertices shared by several geometries are simplified once, so they stay shared, and shared edges count once.
/// The result has the geometries in the same order and with the same data as the input, so handles of the input remain valid.
template<class K, class Data, class PolygonContainer>
CGAL::Straight_geometry_collection_2<K, Data, PolygonContainer> simplify_KSBB(const CGAL::Straight_geometry_collection_2<K, Data, PolygonContainer>& collection, size_t complexity, const Simplification_options<K>& options = {}) {
	return internal::simplify_collection<Simplification_KSBB<false, K, Data, PolygonContainer>>(collection, complexity, options);
}

// ---- Simplifying a single geometry ---- //
template<class K, class PolygonContainer>
CGAL::Polygon_2<K, PolygonContainer> simplify_KSBB(const CGAL::Polygon_2<K, PolygonContainer>& polygon, size_t complexity, const Simplification_options<K>& options = {}) {
	return internal::simplify_geometry<Simplification_KSBB<false, K, std::monostate, PolygonContainer>>(polygon, complexity, options);
}

template<class K, class PolygonContainer>
CGAL::Polygon_with_holes_2<K, PolygonContainer> simplify_KSBB(const CGAL::Polygon_with_holes_2<K, PolygonContainer>& polygon, size_t complexity, const Simplification_options<K>& options = {}) {
	return internal::simplify_geometry<Simplification_KSBB<false, K, std::monostate, PolygonContainer>>(polygon, complexity, options);
}

template<class K, class PolygonContainer>
CGAL::Multipolygon_with_holes_2<K, PolygonContainer> simplify_KSBB(const CGAL::Multipolygon_with_holes_2<K, PolygonContainer>& multipolygon, size_t complexity, const Simplification_options<K>& options = {}) {
	CGAL::Multipolygon_with_holes_2<K, PolygonContainer> input;
	for (const auto& polygon : multipolygon.polygons_with_holes()) {
		input.add_polygon_with_holes(polygon);
	}
	return internal::simplify_geometry<Simplification_KSBB<false, K, std::monostate, PolygonContainer>>(input, complexity, options);
}

template<class K>
CGAL::Polyline_2<K> simplify_KSBB(const CGAL::Polyline_2<K>& polyline, size_t complexity, const Simplification_options<K>& options = {}) {
	return internal::simplify_geometry<Simplification_KSBB<false, K>>(polyline, complexity, options);
}

template<class K>
CGAL::Multipolyline_2<K> simplify_KSBB(const CGAL::Multipolyline_2<K>& multipolyline, size_t complexity, const Simplification_options<K>& options = {}) {
	return internal::simplify_geometry<Simplification_KSBB<false, K>>(multipolyline, complexity, options);
}
}
