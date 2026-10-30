#pragma once

#include <CGAL/Simplification_2/graph_simplification/vertex_quad_tree.h>

#include <CGAL/Simplification_2/geometries/Straight_geometry_collection_2.h>

#include <functional>
#include <iostream>

namespace CGAL {
/// Options for the simplification algorithms.
template <class K>
struct Simplification_options {
	/// Points of different geometries within L-infinity distance prec are merged into one vertex.
	typename K::FT prec = 0;
};

namespace internal {
/// Common base class for algorithms that work on straight-edge geometries and have a graph as underlying subdivision representation.
///
/// All geometries of a collection are inserted so that boundaries shared by several geometries are simplified
/// only once and remain shared, and points are respected as obstacles.
///
/// \section complexity Complexity
/// The complexity is the number of edges where an edge shared by several geometries counts once, and points do not count.
///
/// @tparam K 2D linear kernel.
/// @tparam GraphSimplificationAlgorithm the simplification algorithm that operates on the graph.
/// @tparam Data data stored with each geometry, or std::monostate for none.
/// @tparam PolygonContainer the container of the polygon vertices.
template<class K, class GraphSimplificationAlgorithm, class Data = std::monostate, class PolygonContainer = std::vector<CGAL::Point_2<K>>>
class Straight_simplification_algorithm {
  public:
	using Collection = CGAL::Straight_geometry_collection_2<K, Data, PolygonContainer>;

  private:
	static const int m_quadtree_depth = 10;
	const typename K::FT m_prec;
	using Point = K::Point_2;

	// GraphTraits of the Graph should have decomposed set to true.
	using Graph = typename GraphSimplificationAlgorithm::Graph;
	using Vertex_handle = typename Graph::Vertex_handle;
	using Edge_handle = typename Graph::Edge_handle;
	using Path_handle = typename Graph::Path_handle;

	struct Directed_path {
		Path_handle path_handle;
		bool reversed;

		Edge_handle first_edge() const {
			return reversed ? path_handle->last_edge() : path_handle->first_edge();
		}

		Edge_handle last_edge() const {
			return reversed ? path_handle->first_edge() : path_handle->last_edge();
		}

		/// The vertex at which the traversal enters the given edge of this path.
		static Vertex_handle entry(Edge_handle e, bool reversed) {
			return reversed ? e->target() : e->source();
		}

		/// The vertex at which the traversal leaves the path.
		Vertex_handle end_vertex() const {
			return reversed ? last_edge()->source() : last_edge()->target();
		}

		/// Calls f on each vertex of the path in traversal order, except the vertex at which the path is left.
		/// For a cyclic path, every vertex is visited once.
		template <class F>
		void for_each_vertex_but_last(F f) const {
			Edge_handle e = first_edge();
			while (true) {
				f(entry(e, reversed));
				if (e == last_edge()) {
					break;
				}
				e = reversed ? e->prev() : e->next();
			}
		}
	};

	// The paths that each geometry consists of, mirroring the types of CGAL::Straight_geometry_2.
	// Points are not simplified, so they are stored as is.
	struct Polyline_paths {
		std::vector<Directed_path> chain;
	};
	struct Multipolyline_paths {
		std::vector<Polyline_paths> polylines;
	};
	struct Polygon_paths {
		std::vector<Directed_path> ring;
	};
	struct Polygon_with_holes_paths {
		/// The outer boundary followed by the holes.
		std::vector<Polygon_paths> rings;
	};
	struct Multipolygon_paths {
		std::vector<Polygon_with_holes_paths> polygons;
	};
	using Geometry_paths = std::variant<Point, CGAL::Multipoint_2<K>, Polyline_paths, Multipolyline_paths, Polygon_paths,
	                                   Polygon_with_holes_paths, Multipolygon_paths>;

	using Vertex_tree = CGAL::internal::VertexQuadTree<Graph>;
	std::unique_ptr<Vertex_tree> m_vqt;
	std::unique_ptr<GraphSimplificationAlgorithm> m_alg;
	Graph m_graph;
	std::vector<Geometry_paths> m_paths;
	std::vector<Data> m_data;

	/// Returns the vertex that the given point is merged into, or nullptr if there is none:
	/// the vertex with the smallest index among those within L-infinity distance m_prec.
	///
	/// Insertion uses the same rule and vertices are only appended, so vertices that are added after a point
	/// has been inserted have larger indices than the vertex it was merged into. Hence, looking up a point after
	/// all geometries are inserted returns the same vertex as during its insertion, without storing the merges.
	Vertex_handle find_vertex(const Point& point) {
		if (m_prec == 0) {
			// No two vertices coincide, so there is at most one match.
			return m_vqt->findElement(point, m_prec);
		}
		Vertex_handle best = nullptr;
		CGAL::Iso_rectangle_2<K> query(point.x() - m_prec, point.y() - m_prec, point.x() + m_prec, point.y() + m_prec);
		m_vqt->findContained(query, [&best](Vertex_handle v) {
			if (best == nullptr || v->graph_index() < best->graph_index()) {
				best = v;
			}
		});
		return best;
	}

	Vertex_handle find_or_add_vertex(const Point& point) {
		Vertex_handle v = find_vertex(point);
		if (v == nullptr) {
			v = m_graph.add_vertex(point);
			m_vqt->insert(v);
		}
		return v;
	}

	/// Walks the chain through the given points, calling on_edge(u, v) for every edge from u to v.
	/// Points are mapped to vertices with vertex_of; consecutive points that map to the same vertex are skipped.
	/// If closed, the last point is connected to the first.
	/// Both the insertion and the path computation use this, so that they see the same sequence of edges.
	template <class InputIterator,
		std::invocable<const Point&> VertexOf,
		std::invocable<Vertex_handle, Vertex_handle> OnEdge>
	void walk_chain(InputIterator begin, InputIterator end, bool closed, VertexOf vertex_of, OnEdge on_edge) {
		Vertex_handle first = nullptr;
		Vertex_handle prev = nullptr;
		for (auto it = begin; it != end; ++it) {
			Vertex_handle next = vertex_of(*it);
			if (first == nullptr) {
				first = next;
			}
			if (prev != nullptr && next != prev) {
				on_edge(prev, next);
			}
			prev = next;
		}
		if (closed && first != prev) {
			on_edge(prev, first);
		}
	}

	/// Adds the vertices and edges of the chain through the given points to the graph.
	template <class InputIterator>
	void add_chain(InputIterator begin, InputIterator end, bool closed) {
		walk_chain(begin, end, closed,
		    [this](const Point& p) { return find_or_add_vertex(p); },
		    [this](Vertex_handle u, Vertex_handle v) {
			    if (!u->is_neighbor_of(v)) {
				    m_graph.add_edge(u, v);
			    }
		    });
	}

	/// Returns the paths traversed by the chain through the given points, which must have been added to the graph.
	template <class InputIterator>
	std::vector<Directed_path> find_paths(InputIterator begin, InputIterator end, bool closed) {
		std::vector<Directed_path> paths;
		Edge_handle first_edge = nullptr;
		Vertex_handle first_vertex = nullptr; // where the chain enters the first edge
		Vertex_handle last_vertex = nullptr;  // where the chain leaves the last edge

		walk_chain(begin, end, closed,
		    [this](const Point& p) {
			    Vertex_handle v = find_vertex(p);
			    assert(v != nullptr);
			    return v;
		    },
		    [&](Vertex_handle u, Vertex_handle v) {
			    Edge_handle e = u->find_edge_to(v);
			    assert(e != nullptr);
			    if (first_edge == nullptr) {
				    first_edge = e;
				    first_vertex = u;
			    }
			    last_vertex = v;
			    // A new path starts at every vertex whose degree is not two.
			    if (u->degree() != 2) {
				    Path_handle p = e->path();
				    bool reversed = e->source() != u;
				    assert(e == (reversed ? p->last_edge() : p->first_edge()));
				    paths.push_back({p, reversed});
			    }
		    });

		if (first_edge == nullptr) {
			return paths;
		}
		if (first_vertex == last_vertex) {
			// A ring. If it starts inside a path, that path is found when the walk reaches the start of the path,
			// so the paths are still in cyclic order.
			if (paths.empty()) {
				// All vertices have degree two, so the ring is a single cyclic path.
				assert(first_edge->path()->cyclic());
				paths.push_back({first_edge->path(), first_edge->source() != first_vertex});
			}
		} else if (first_vertex->degree() == 2 || last_vertex->degree() == 2) {
			throw std::runtime_error("Polylines that start or end at a degree-2 vertex are not supported.");
		}
		return paths;
	}

	template <class Ring>
	void add_ring(const Ring& ring) {
		add_chain(ring.vertices_begin(), ring.vertices_end(), true);
	}

	void add_rings(const CGAL::Polygon_with_holes_2<K, PolygonContainer>& polygon) {
		add_ring(polygon.outer_boundary());
		for (const auto& hole : polygon.holes()) {
			add_ring(hole);
		}
	}

	/// Warns if the point coincides with a vertex of a polygon or polyline: that vertex is not pinned,
	/// so the simplification may remove or move it, after which the point no longer lies on the geometry.
	void warn_if_on_vertex(const Point& point) {
		Vertex_handle v = find_vertex(point);
		assert(v != nullptr);
		if (v->degree() > 0) {
			std::cerr << "Warning: point (" << point << ") coincides with a vertex of a polygon or polyline. "
			          << "The vertex is not pinned: the simplification may remove or move it, "
			          << "after which the point no longer lies on that geometry." << std::endl;
		}
	}

	// ---- Inserting geometries into the graph ---- //

	void add(const CGAL::Polygon_2<K, PolygonContainer>& polygon) {
		add_ring(polygon);
	}

	void add(const CGAL::Polygon_with_holes_2<K, PolygonContainer>& polygon) {
		add_rings(polygon);
	}

	void add(const CGAL::Multipolygon_with_holes_2<K, PolygonContainer>& multipolygon) {
		for (const auto& polygon : multipolygon.polygons_with_holes()) {
			add_rings(polygon);
		}
	}

	void add(const CGAL::Polyline_2<K>& polyline) {
		add_chain(polyline.vertices_begin(), polyline.vertices_end(), false);
	}

	void add(const CGAL::Multipolyline_2<K>& multipolyline) {
		for (const auto& polyline : multipolyline.polylines) {
			add_chain(polyline.vertices_begin(), polyline.vertices_end(), false);
		}
	}

	/// Points are added as isolated vertices, so that the simplification does not move geometry across them.
	void add(const CGAL::internal::Point<K>& point) {
		find_or_add_vertex(point);
	}

	void add(const CGAL::Multipoint_2<K>& multipoint) {
		for (const auto& point : multipoint.points) {
			find_or_add_vertex(point);
		}
	}

	// ---- Computing the paths of geometries after the graph is initialized ---- //

	Polygon_paths compute_paths(const CGAL::Polygon_2<K, PolygonContainer>& polygon) {
		return {find_paths(polygon.vertices_begin(), polygon.vertices_end(), true)};
	}

	Polygon_with_holes_paths compute_paths(const CGAL::Polygon_with_holes_2<K, PolygonContainer>& polygon) {
		Polygon_with_holes_paths paths;
		paths.rings.push_back(compute_paths(polygon.outer_boundary()));
		for (const auto& hole : polygon.holes()) {
			paths.rings.push_back(compute_paths(hole));
		}
		return paths;
	}

	Multipolygon_paths compute_paths(const CGAL::Multipolygon_with_holes_2<K, PolygonContainer>& multipolygon) {
		Multipolygon_paths paths;
		for (const auto& polygon : multipolygon.polygons_with_holes()) {
			paths.polygons.push_back(compute_paths(polygon));
		}
		return paths;
	}

	Polyline_paths compute_paths(const CGAL::Polyline_2<K>& polyline) {
		return {find_paths(polyline.vertices_begin(), polyline.vertices_end(), false)};
	}

	Multipolyline_paths compute_paths(const CGAL::Multipolyline_2<K>& multipolyline) {
		Multipolyline_paths paths;
		for (const auto& polyline : multipolyline.polylines) {
			paths.polylines.push_back(compute_paths(polyline));
		}
		return paths;
	}

	Point compute_paths(const CGAL::internal::Point<K>& point) {
		warn_if_on_vertex(point);
		return point;
	}

	CGAL::Multipoint_2<K> compute_paths(const CGAL::Multipoint_2<K>& multipoint) {
		for (const auto& point : multipoint.points) {
			warn_if_on_vertex(point);
		}
		return multipoint;
	}

	// ---- Converting paths back to geometries ---- //

	CGAL::Polygon_2<K, PolygonContainer> to_geometry(const Polygon_paths& paths) const {
		// Consecutive paths share a vertex, which is output once.
		CGAL::Polygon_2<K, PolygonContainer> polygon;
		for (const auto& path : paths.ring) {
			path.for_each_vertex_but_last([&polygon](Vertex_handle v) { polygon.push_back(v->point()); });
		}
		return polygon;
	}

	CGAL::Polygon_with_holes_2<K, PolygonContainer> to_geometry(const Polygon_with_holes_paths& paths) const {
		CGAL::Polygon_with_holes_2<K, PolygonContainer> polygon(to_geometry(paths.rings.front()));
		for (size_t i = 1; i < paths.rings.size(); ++i) {
			polygon.add_hole(to_geometry(paths.rings[i]));
		}
		return polygon;
	}

	CGAL::Multipolygon_with_holes_2<K, PolygonContainer> to_geometry(const Multipolygon_paths& paths) const {
		CGAL::Multipolygon_with_holes_2<K, PolygonContainer> multipolygon;
		for (const auto& polygon : paths.polygons) {
			multipolygon.add_polygon_with_holes(to_geometry(polygon));
		}
		return multipolygon;
	}

	/// If the polyline is closed, its first point is repeated at the end.
	CGAL::Polyline_2<K> to_geometry(const Polyline_paths& paths) const {
		CGAL::Polyline_2<K> polyline;
		for (const auto& path : paths.chain) {
			path.for_each_vertex_but_last([&polyline](Vertex_handle v) { polyline.push_back(v->point()); });
		}
		if (!paths.chain.empty()) {
			polyline.push_back(paths.chain.back().end_vertex()->point());
		}
		return polyline;
	}

	CGAL::Multipolyline_2<K> to_geometry(const Multipolyline_paths& paths) const {
		CGAL::Multipolyline_2<K> multipolyline;
		for (const auto& polyline : paths.polylines) {
			multipolyline.polylines.push_back(to_geometry(polyline));
		}
		return multipolyline;
	}

	Point to_geometry(const Point& point) const {
		return point;
	}

	CGAL::Multipoint_2<K> to_geometry(const CGAL::Multipoint_2<K>& multipoint) const {
		return multipoint;
	}

  protected:
	/// Like the public constructor, but preprocesses the graph after it is initialized but before the simplification algorithm is initialized.
	/// The result of the preprocessing is the input of the simplification, so it cannot be undone.
	/// The preprocessing must keep the paths of the graph intact.
	Straight_simplification_algorithm(const Collection& collection, typename K::FT prec, const std::function<void(Graph&)>& preprocess)
	    : m_prec(prec), m_vqt(std::make_unique<Vertex_tree>(collection.bbox(), m_quadtree_depth)) {
		// insert geometries into graph
		for (const auto& geometry : collection.geometries()) {
			std::visit([this](const auto& geom) { add(geom); }, geometry);
		}
		if constexpr (!std::is_same_v<Data, std::monostate>) {
			for (const auto& data : collection.data()) {
				m_data.push_back(data);
			}
		}

		// intialize everything and compute the paths by walking the geometries again
		m_graph.initialize();
		for (const auto& geometry : collection.geometries()) {
			m_paths.push_back(std::visit([this](const auto& geom) -> Geometry_paths { return compute_paths(geom); }, geometry));
		}

		if (preprocess) {
			preprocess(m_graph);
			if constexpr (Graph::Graph_traits::historic) {
				m_graph.history().clear();
			}
		}

		// Note that the vertex quad tree does not track changes in the graph during preprocessing.
		m_alg = std::make_unique<GraphSimplificationAlgorithm>(m_graph, std::move(m_vqt));
	};

  public:
	/// Inserts the geometries of the collection and prepares the simplification.
	Straight_simplification_algorithm(const Collection& collection, const Simplification_options<K>& options = {})
	    : Straight_simplification_algorithm(collection, options.prec, {}) {}

	/// Returns the simplified geometries, in the same order as the input, so handles of the input collection remain valid.
	Collection get_collection() const {
		Collection collection;
		for (size_t i = 0; i < m_paths.size(); ++i) {
			std::visit(
			    [&](const auto& paths) {
				    if constexpr (!std::is_same_v<Data, std::monostate>) {
					    collection.insert(to_geometry(paths), m_data[i]);
				    } else {
					    collection.insert(to_geometry(paths));
				    }
			    },
			    m_paths[i]);
		}
		return collection;
	}

	/// Returns the current complexity: the number of edges (edges shared among shapes count only once).
	size_t complexity() const {
		return m_graph.number_of_edges();
	}

	/// Simplifies until the complexity (see \ref complexity) is at most the target complexity.
	/// The result has a higher complexity if no further simplification step is possible.
	///
	/// If the graph stores its history, a target complexity above the current one undoes simplification steps,
	/// until the complexity is at least the target or the input is restored.
	/// Without history, a target complexity above the current one has no effect.
	void simplify(size_t target_complexity) {
		if constexpr (Graph::Graph_traits::historic) {
			auto& hist = m_graph.history();

			if (target_complexity > m_graph.number_of_edges()) {
				// revert
				while (hist.can_undo() && m_graph.number_of_edges() < target_complexity) {
					hist.undo();
				}
			}
			else if (target_complexity < m_graph.number_of_edges()) {
				while (hist.can_redo() && m_graph.number_of_edges() > target_complexity) {
					hist.redo();
				}

				if (!hist.can_redo() && target_complexity < m_graph.number_of_edges()) {
					// already at present, run algorithm further
					m_alg->run([&](int complexity, typename K::FT cost) {
						return complexity <= target_complexity;
					});
				}
			}
		}
		else {
			if (target_complexity < m_graph.number_of_edges()) {
				// already at present, run algorithm further
				m_alg->run([&](int complexity, typename K::FT cost) {
					return complexity <= target_complexity;
				});
			}
		}
	}
};

/// Simplifies a single geometry with the given algorithm until the number of edges is at most the target number.
template <class Algorithm, class Geometry, class Options>
Geometry simplify_geometry(const Geometry& geometry, size_t complexity, const Options& options) {
	typename Algorithm::Collection collection;
	auto handle = collection.insert(geometry);
	Algorithm algorithm(collection, options);
	algorithm.simplify(complexity);
	auto simplified = algorithm.get_collection();
	return simplified.get_geometry(handle);
}

/// Simplifies all geometries of a collection together with the given algorithm until the number of edges is at most the target number.
template <class Algorithm, class Options>
typename Algorithm::Collection simplify_collection(const typename Algorithm::Collection& collection, size_t complexity, const Options& options) {
	Algorithm algorithm(collection, options);
	algorithm.simplify(complexity);
	return algorithm.get_collection();
}
} // namespace internal
} // namespace CGAL
