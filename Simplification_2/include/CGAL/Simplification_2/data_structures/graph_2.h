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

#include <functional>

#include "graph_curve_traits_2.h"
#include "graph_map_2.h"
#include "graph_operation_2.h"
#include "graph_traits_2.h"

namespace CGAL::internal {

template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits, GraphTraits_2 GraphTraits>
class Graph_2_vertex;
template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits, GraphTraits_2 GraphTraits>
class Graph_2_edge;
template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits, GraphTraits_2 GraphTraits>
class Graph_2_path;

template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits, GraphTraits_2 GraphTraits>
class Graph_2 {
	friend class Graph_2_vertex<VertexData, EdgeData, CurveTraits, GraphTraits>;
	friend class Graph_2_edge<VertexData, EdgeData, CurveTraits, GraphTraits>;
	template <typename G, typename T> friend class Graph_vertex_map;
	template <typename G, typename T> friend class Graph_edge_map;
	template <typename G, typename T> friend class Graph_path_map;

	template <class InGraph, class OutGraph>
	friend void
	graph_2_copy(InGraph& input, OutGraph& output,
	             const std::function<typename OutGraph::Curve_traits::Curve_representation_2(
	                 const typename InGraph::Curve_traits::Curve_representation_2&)>& conversion,
	             bool requireResorting);

	// basic operations
	template <class G> friend class detail::AddVertex;
	template <class G> friend class detail::RemoveVertex;
	template <class G> friend class detail::AddEdge;
	template <class G> friend class detail::RemoveEdge;
	// geometric operations
	template <class G> friend class detail::ChangeCurve;
	template <class G> friend class detail::MoveVertex;
	// oriented operations
	template <class G> friend class detail::MergeVertex;
	template <class G> friend class detail::SplitVertex;
	template <class G> friend class detail::CollapseEdge;
	template <class G> friend class detail::SubdivideEdge;

  public:
	using Vertex = Graph_2_vertex<VertexData, EdgeData, CurveTraits, GraphTraits>;
	using Vertex_handle = Vertex*;
	using Vertex_const_handle = const Vertex*;

	using Edge = Graph_2_edge<VertexData, EdgeData, CurveTraits, GraphTraits>;
	using Edge_handle = Edge*;
	using Edge_const_handle = const Edge*;

	using Vertex_data = VertexData;
	using Edge_data = EdgeData;

	using Curve_traits = CurveTraits;
	using Kernel = CurveTraits::Kernel;
	using Point_2 = CurveTraits::Point_2;
	using Curve_2 = CurveTraits::Curve_2;
	using Curve_representation = CurveTraits::Curve_representation_2;

	using Graph_traits = GraphTraits;
	using History = std::conditional<GraphTraits::historic, GraphHistory, NoGraphHistory>::type;

  private:
	using Vertex_container = std::vector<Vertex_handle>;
	using Edge_container = std::vector<Edge_handle>;

  public:
	using Vertex_iterator = Vertex_container::iterator;
	using Vertex_const_iterator = Vertex_container::const_iterator;
	using Edge_iterator = Edge_container::iterator;
	using Edge_const_iterator = Edge_container::const_iterator;

  private:
	Vertex_container m_vertices;
	Edge_container m_edges;
	History m_history;
	bool m_initialized = false;

	std::vector<Graph_vertex_index_listener*> m_vertex_listeners;
	std::vector<Graph_edge_index_listener*> m_edge_listeners;
	std::vector<Graph_path_index_listener*> m_path_listeners;

	void insert_vertex_into_container(Vertex_handle v) {
		const size_t index = v->graph_index();
		if (index == m_vertices.size()) {
			m_vertices.push_back(v);
			for (Graph_vertex_index_listener* m : m_vertex_listeners) {
				m->add_index();
			}
		} else {
			m_vertices[index]->m_index = m_vertices.size();
			m_vertices.push_back(m_vertices[index]);
			m_vertices[index] = v;
			for (Graph_vertex_index_listener* m : m_vertex_listeners) {
				m->add_index(index);
			}
		}
	}

	void remove_vertex_from_container(Vertex_handle v) {
		const size_t index = v->m_index;
		const size_t last = m_vertices.size() - 1;
		if (index != last) {
			m_vertices[index] = m_vertices[last];
			m_vertices[index]->m_index = index;
			for (Graph_vertex_index_listener* m : m_vertex_listeners) {
				m->remove_index(index);
			}
		} else {
			for (Graph_vertex_index_listener* m : m_vertex_listeners) {
				m->remove_last_index();
			}
		}
		m_vertices.pop_back();
	}

	void insert_edge_into_container(Edge_handle e) {
		const size_t index = e->graph_index();
		if (index == m_edges.size()) {
			m_edges.push_back(e);
			for (Graph_edge_index_listener* m : m_edge_listeners) {
				m->add_index();
			}
		} else {
			m_edges[index]->m_index = m_edges.size();
			m_edges.push_back(m_edges[index]);
			m_edges[index] = e;
			for (Graph_edge_index_listener* m : m_edge_listeners) {
				m->add_index(index);
			}
		}
	}

	void remove_edge_from_container(Edge_handle e) {
		const size_t index = e->m_index;
		const size_t last = m_edges.size() - 1;
		if (index != last) {
			m_edges[index] = m_edges[last];
			m_edges[index]->m_index = index;
			for (Graph_edge_index_listener* m : m_edge_listeners) {
				m->remove_index(index);
			}
		} else {
			for (Graph_edge_index_listener* m : m_edge_listeners) {
				m->remove_last_index();
			}
		}
		m_edges.pop_back();
	}

	inline void ensure_traits(Vertex_handle v) {
		if constexpr (GraphTraits::oriented) {
			orient(v);
		}
		if constexpr (GraphTraits::sorted) {
			sort_incident_edges(v);
		}
	}

	bool verify_traits() const {
		bool correct = true;
		if constexpr (GraphTraits::oriented) {
			if (!verify_oriented())
				correct = false;
		}
		if constexpr (GraphTraits::sorted) {
			if (!verify_sorted())
				correct = false;
		}
		if constexpr (GraphTraits::decomposed) {
			if (!verify_decomposed()) {
				correct = false;
			}
		}
		return correct;
	}

	void ensure_oriented() requires Graph_traits::oriented {
		for (Vertex_handle v : m_vertices) {
			orient(v);
		}
	}

	void orient(Vertex_handle v) requires Graph_traits::oriented {
		if (v->degree() != 2)
			return; // irrelevant for orientation

		Edge_handle bwd = v->m_incident[0];
		Edge_handle fwd = v->m_incident[1];

		if (bwd->m_target == v && fwd->m_source == v) {
			return; // already satisfies orientation
		}

		if (fwd->m_target == v && bwd->m_source == v) {
			std::swap(v->m_incident[0], v->m_incident[1]);
			return;
		}

		if (bwd->m_target != v) {
			bwd->reverse();
		}

		if (fwd->m_source != v) {
			fwd->reverse();
		}

		Vertex_handle fv = fwd->m_target;
		while (fv->degree() == 2 && fv != v) {
			if (fv->m_incident[0] != fwd) {
				fv->m_incident[1] = fv->m_incident[0];
				fv->m_incident[0] = fwd;
			}

			fwd = fv->m_incident[1];
			if (fwd->m_source != fv) {
				fwd->reverse();
			}
			fv = fwd->m_target;
		}

		if (fv != v) {
			Vertex_handle bv = bwd->m_source;
			while (bv->degree() == 2) {
				if (bv->m_incident[1] != bwd) {
					bv->m_incident[0] = bv->m_incident[1];
					bv->m_incident[1] = bwd;
				}

				bwd = bv->m_incident[0];
				if (bwd->m_target != bv) {
					bwd->reverse();
				}
				bv = bwd->m_source;
			}
		}
	}

	bool verify_oriented() const requires Graph_traits::oriented {
		bool correct = true;
		for (Vertex_handle v : m_vertices) {
			if (!verify_oriented(v)) {
				std::cout << "! Incorrectly oriented vertex: " << *v << std::endl;
				correct = false;
			}
		}

		return correct;
	}
	bool verify_oriented(Vertex_const_handle v) const requires Graph_traits::oriented {
		if (v->degree() != 2)
			return true; // irrelevant for orientation

		auto bwd = v->m_incident[0];
		auto fwd = v->m_incident[1];

		return bwd->m_target == v && fwd->m_source == v;
	}

	void ensure_sorted() requires Graph_traits::sorted {
		for (Vertex_handle v : m_vertices) {
			sort_incident_edges(v);
		}
	}

	void sort_incident_edges(Vertex_handle v) requires Graph_traits::sorted {
		if (v->degree() > 2) {
			std::ranges::sort(v->m_incident, [&v](Edge_handle e, Edge_handle f) {
				CGAL::Direction_2<Kernel> dir_e =
				    CGAL::Direction_2<Kernel>(e->direction_at(v));
				CGAL::Direction_2<Kernel> dir_f =
				    CGAL::Direction_2<Kernel>(f->direction_at(v));
				return dir_e < dir_f;
			});
		}
	}
	bool verify_sorted() const requires Graph_traits::sorted {
		bool correct = true;
		for (Vertex_const_handle v : m_vertices) {
			if (!verify_sorted(v)) {
				std::cout << "! Unsorted vertex: " << *v << std::endl;
				correct = false;
			}
		}
		return correct;
	}
	bool verify_sorted(Vertex_const_handle v) const requires Graph_traits::sorted {
		const size_t d = v->degree();
		if (d > 2) {
			using Dir = CGAL::Direction_2<Kernel>;
			// TODO: this should use the tangents of the curve; in this form, it only works for straightline graphs
			size_t min_i = 0;
			Dir min_dir = v->incident_edge(0)->direction_at(v);
			for (size_t i = 1; i < d; ++i) {
				Dir dir = v->incident_edge(i)->direction_at(v);
				if (dir < min_dir) {
					min_dir = dir;
					min_i = i;
				}
			}

			Dir prev = min_dir;
			for (size_t i = 1; i < d; ++i) {
				Dir dir = v->incident_edge((min_i + i) % d)->direction_at(v);
				if (prev > dir) {
					return false;
				}
				prev = min_dir;
			}
		}
		return true;
	}

	// ---- ONLY WHEN GraphTraits::decomposed -------------------------------- //
  public:
	using Path = std::conditional<GraphTraits::decomposed,
	                              Graph_2_path<VertexData, EdgeData, CurveTraits, GraphTraits>,
	                              std::monostate>::type;
	using Path_handle = Path*;
	using Path_const_handle = const Path*;

	using Path_data = GraphTraits::PathData;

  private:
	using Path_container = std::vector<Path_handle>;

	Path_container m_paths;

	void ensure_decomposed() requires GraphTraits::decomposed {
		for (Edge_handle e : m_edges) {
			if (e->m_path != nullptr)
				continue; // already handled

			Path_handle p = new Path(e, m_paths.size());
			m_paths.push_back(p);

			e->m_path = p;

			if (e->m_source->degree() == 2) {
				//NB: we check first whether we can move at all so we can detect a cycle from the start pointer

				do {
					p->m_start = p->m_start->prev();
					p->m_start->m_path = p;
				} while (p->m_start->m_source->degree() == 2 && p->m_start != e);

				if (p->m_start == e) {
					p->m_cyclic = true;
					p->m_end = e->prev();
					continue; // cycle done, continues mainloop
				}
			}

			while (p->m_end->m_target->degree() == 2) {
				p->m_end = p->m_end->next();
				p->m_end->m_path = p;
			}
		}

		for (Graph_path_index_listener* list : m_path_listeners) {
			list->resize(number_of_paths());
		}

		assert(verify_decomposed());
	}
	bool verify_decomposed() const requires GraphTraits::decomposed {
		bool correct = true;
		for (Edge_handle e : m_edges) {
			if (e->m_path == nullptr) {
				std::cout << "! Edge without path: " << *e << std::endl;
				correct = false;
			} else {
				if (e->source()->degree() == 2 && e->prev()->m_path != e->m_path) {
					std::cout << "! Previous edge does not match path: " << *e << std::endl;
					correct = false;
				}
				if (e->target()->degree() == 2 && e->next()->m_path != e->m_path) {
					std::cout << "! Next edge does not match path: " << *e << std::endl;
					correct = false;
				}
			}
		}
		for (Path_handle p : m_paths) {
			if (p->m_start->m_path != p) {
				std::cout << "! First edge of path not on path: " << *p << std::endl;
				correct = false;
			}
			if (p->m_end->m_path != p) {
				std::cout << "! Last edge of path not on path: " << *p << std::endl;
				correct = false;
			}
			if (p->m_cyclic) {
				if (p->m_start->source()->degree() != 2) {
					std::cout << "! Cyclic path does not start with degree-2 vertex: " << *p
					          << std::endl;
					correct = false;
				}
				if (p->m_end->target() != p->m_start->source()) {
					std::cout << "! Cyclic path does not start and end at same vertex: " << *p
					          << std::endl;
					correct = false;
				}
			} else {
				if (p->m_start->source()->degree() == 2) {
					std::cout << "! Acyclic path starts at degree-2 vertex: " << *p << std::endl;
					correct = false;
				}
				if (p->m_end->target()->degree() == 2) {
					std::cout << "! Acyclic path ends at degree-2 vertex: " << *p << std::endl;
					correct = false;
				}
			}
		}
		return correct;
	}

	bool verify_graph_structure() const {
		for (Vertex_handle v : m_vertices) {
			if (v->m_index >= m_vertices.size()) {
				return false;
			} else if (m_vertices[v->m_index] != v) {
				return false;
			}

			for (Edge_handle e : v->m_incident) {
				if (m_edges[e->m_index] != e) {
					return false;
				}
			}
		}

		for (Edge_handle e : m_edges) {
			if (e->m_index >= m_edges.size()) {
				return false;
			} else if (m_edges[e->m_index] != e) {
				return false;
			}
		}

		return true;
	}

  public:
	using Path_iterator = Path_container::iterator;
	using Path_const_iterator = Path_container::const_iterator;

	class Path_range {
		friend class Graph_2<VertexData, EdgeData, CurveTraits, GraphTraits>;

	  private:
		Path_container& m_container;

		Path_range(Path_container& container) : m_container(container) {}

	  public:
		Path_iterator begin() const {
			return m_container.begin();
		}

		Path_iterator end() const {
			return m_container.end();
		}
	};

	class Path_const_range {
		friend class Graph_2<VertexData, EdgeData, CurveTraits, GraphTraits>;

	  private:
		const Path_container& m_container;

		Path_const_range(const Path_container& container) : m_container(container) {}

	  public:
		Path_const_iterator begin() const {
			return m_container.cbegin();
		}

		Path_const_iterator end() const {
			return m_container.cend();
		}
	};

	size_t number_of_paths() const requires GraphTraits::decomposed {
		return m_paths.size();
	}
	Path_handle path(size_t index) requires GraphTraits::decomposed {
		return m_paths[index];
	}
	Path_const_handle path(size_t index) const requires GraphTraits::decomposed {
		return m_paths[index];
	}
	Path_range paths() requires GraphTraits::decomposed {
		return Path_range(m_paths);
	}
	Path_const_range paths() const requires GraphTraits::decomposed {
		return Path_const_range(m_paths);
	}
	Path_iterator paths_begin() requires GraphTraits::decomposed {
		return m_paths.begin();
	}
	Path_iterator paths_end() requires GraphTraits::decomposed {
		return m_paths.end();
	}
	Path_const_iterator paths_begin() const requires GraphTraits::decomposed {
		return m_paths.cbegin();
	}
	Path_const_iterator paths_end() const requires GraphTraits::decomposed {
		return m_paths.cend();
	}

	// ----Vertex/edge listeners --------------------------------------------------- //
	void add_listener(Graph_vertex_index_listener* list) {
		m_vertex_listeners.push_back(list);
		list->resize(number_of_vertices());
	}

	void remove_listener(Graph_vertex_index_listener* list) {
		auto pos = std::find(m_vertex_listeners.begin(), m_vertex_listeners.end(), list);
		if (pos != m_vertex_listeners.end()) {
			m_vertex_listeners.erase(pos);
		}
	}

	void add_listener(Graph_edge_index_listener* list) {
		m_edge_listeners.push_back(list);
		list->resize(number_of_edges());
	}

	void remove_listener(Graph_edge_index_listener* list) {
		auto pos = std::find(m_edge_listeners.begin(), m_edge_listeners.end(), list);
		if (pos != m_edge_listeners.end()) {
			m_edge_listeners.erase(pos);
		}
	}

	void add_listener(Graph_path_index_listener* list) requires GraphTraits::decomposed {
		m_path_listeners.push_back(list);
		list->resize(number_of_paths());
	}

	void remove_listener(Graph_path_index_listener* list) requires GraphTraits::decomposed {
		auto pos = std::find(m_path_listeners.begin(), m_path_listeners.end(), list);
		if (pos != m_path_listeners.end()) {
			m_path_listeners.erase(pos);
		}
	}

  public:
	class Vertex_range {
		friend class Graph_2<VertexData, EdgeData, CurveTraits, GraphTraits>;

	  private:
		Vertex_container& m_container;

		Vertex_range(Vertex_container& container) : m_container(container) {}

	  public:
		Vertex_iterator begin() const {
			return m_container.begin();
		}

		Vertex_iterator end() const {
			return m_container.end();
		}
	};

	class Vertex_const_range {
		friend class Graph_2<VertexData, EdgeData, CurveTraits, GraphTraits>;

	  private:
		const Vertex_container& m_container;

		Vertex_const_range(const Vertex_container& container) : m_container(container) {}

	  public:
		Vertex_const_iterator begin() const {
			return m_container.cbegin();
		}

		Vertex_const_iterator end() const {
			return m_container.cend();
		}
	};

	class Edge_range {
		friend class Graph_2<VertexData, EdgeData, CurveTraits, GraphTraits>;

	  private:
		Edge_container& m_container;

		Edge_range(Edge_container& container) : m_container(container) {}

	  public:
		Edge_iterator begin() const {
			return m_container.begin();
		}

		Edge_iterator end() const {
			return m_container.end();
		}
	};

	class Edge_const_range {
		friend class Graph_2<VertexData, EdgeData, CurveTraits, GraphTraits>;

	  private:
		const Edge_container& m_container;

		Edge_const_range(const Edge_container& container) : m_container(container) {}

	  public:
		Edge_const_iterator begin() const {
			return m_container.cbegin();
		}

		Edge_const_iterator end() const {
			return m_container.cend();
		}
	};

	Graph_2& operator=(const Graph_2& other) {
		if (this == &other)
			return *this;

		clear();

		const size_t num_v = other.number_of_vertices();
		m_vertices.resize(num_v);
		for (Graph_vertex_index_listener* list : m_vertex_listeners) {
			list->resize(num_v);
		}

		const size_t num_e = other.number_of_edges();
		m_edges.resize(num_e);
		for (Graph_edge_index_listener* list : m_edge_listeners) {
			list->resize(num_e);
		}

		m_initialized = other.m_initialized;

		Graph_static_vertex_map<Graph_2, Vertex_handle> vmap(other, nullptr);

		for (auto vit : other.vertices()) {
			Vertex_handle new_v = new Vertex(vit->m_point, vit->m_index, vit->m_data);
			new_v->m_index = vit->m_index;
			m_vertices[vit->m_index] = new_v;

			vmap[vit] = new_v;
		}

		Graph_static_edge_map<Graph_2, Edge_handle> emap(other, nullptr);

		for (auto eit : other.edges()) {
			Vertex_handle new_source = vmap[eit->m_source];
			Vertex_handle new_target = vmap[eit->m_target];

			Edge_handle new_e =
			    new Edge(new_source, new_target, eit->m_index, eit->m_representation, eit->m_data);
			m_edges[eit->m_index] = new_e;

			emap[eit] = new_e;
		}

		for (auto vit : other.vertices()) {
			Vertex_handle new_v = vmap[vit];
			for (auto old_e : vit->m_incident) {
				new_v->m_incident.push_back(emap[old_e]);
			}
		}

		if constexpr (GraphTraits::decomposed) {
			const size_t num_p = other.number_of_paths();
			m_paths.resize(num_p);
			for (Graph_path_index_listener* list : m_path_listeners) {
				list->resize(num_p);
			}

			Graph_static_path_map<Graph_2, Path_handle> pmap(other, nullptr);

			for (auto pit : other.paths()) {
				Path_handle new_p = new Path(emap[pit->m_start], emap[pit->m_end], pit->m_cyclic,
				                             pit->m_index, pit->m_data);
				m_paths[new_p->m_index] = new_p;
			}

			for (auto eit : other.edges()) {
				emap[eit]->m_path = pmap[eit->m_path];
			}
		}

		return *this;
	}

	Graph_2() = default;

	Graph_2(const Graph_2& other) {
		*this = other;
	}

	~Graph_2() {
		for (Vertex_handle v : m_vertices) {
			delete v;
		}
		for (Edge_handle e : m_edges) {
			delete e;
		}
		if constexpr (GraphTraits::decomposed) {
			for (Path_handle p : m_paths) {
				delete p;
			}
		}
	}

	Graph_2 transform(CGAL::Aff_transformation_2<Kernel> trans) const {
		Graph_2 transformed;

		transformed.m_vertices.resize(number_of_vertices());
		transformed.m_edges.resize(number_of_edges());

		transformed.m_initialized = m_initialized;

		Graph_static_vertex_map<Graph_2, Vertex_handle> vmap(this, nullptr);

		for (auto vit : vertices()) {
			Vertex_handle new_v = new Vertex(vit->m_point.transform(trans), vit->m_data);
			new_v->m_index = vit->m_index;
			transformed.m_vertices[vit->m_index] = new_v;

			vmap[vit] = new_v;
		}

		Graph_static_edge_map<Graph_2, Edge_handle> emap(this, nullptr);

		for (auto eit : edges()) {
			Vertex_handle new_source = vmap[eit->m_source];
			Vertex_handle new_target = vmap[eit->m_target];

			Edge_handle new_e = new Edge(
			    new_source, new_target,
			    Curve_traits::transform(trans, new_source, new_target, eit->m_representation),
			    eit->m_data);
			new_e->m_index = eit->m_index;
			transformed.m_edges[eit->m_index] = new_e;

			emap[eit] = new_e;
		}

		for (auto vit : vertices()) {
			Vertex_handle new_v = vmap[vit];
			for (auto old_e : vit->m_incident) {
				new_v->m_incident.push_back(emap[old_e]);
			}
		}

		if constexpr (GraphTraits::decomposed) {
			transformed.m_paths.resize(number_of_paths());

			Graph_static_path_map<Graph_2, Path_handle> pmap(this, nullptr);

			for (auto pit : paths()) {
				Path_handle new_p = new Path(emap[pit->m_source], emap[pit->m_target], pit->m_cyclc,
				                             pit->m_index, pit->m_data);
				transformed.m_paths[new_p->m_index] = new_p;
			}

			for (auto eit : edges()) {
				emap[eit]->m_path = pmap[eit->m_path];
			}
		}

		return transformed;
	}

	void initialize() {
		if (!m_initialized) {
			if constexpr (GraphTraits::oriented || GraphTraits::sorted) {
				for (Vertex_handle v : m_vertices) {
					ensure_traits(v);
				}
			}
			if constexpr (GraphTraits::decomposed) {
				ensure_decomposed();
			}
			m_initialized = true;
		}
	}

	bool is_initialized() {
		assert(!m_initialized || verify_traits());
		return m_initialized;
	}

	size_t number_of_vertices() const {
		return m_vertices.size();
	}
	Vertex_handle vertex(size_t index) {
		return m_vertices[index];
	}
	Vertex_const_handle vertex(size_t index) const {
		return m_vertices[index];
	}
	Vertex_range vertices() {
		return Vertex_range(m_vertices);
	}
	Vertex_const_range vertices() const {
		return Vertex_const_range(m_vertices);
	}
	Vertex_iterator vertices_begin() {
		return m_vertices.begin();
	}
	Vertex_iterator vertices_end() {
		return m_vertices.end();
	}
	Vertex_const_iterator vertices_begin() const {
		return m_vertices.cbegin();
	}
	Vertex_const_iterator vertices_end() const {
		return m_vertices.cend();
	}

	size_t number_of_edges() const {
		return m_edges.size();
	}
	Edge_handle edge(size_t index) {
		return m_edges[index];
	}
	Edge_const_handle edge(size_t index) const {
		return m_edges[index];
	}
	Edge_range edges() {
		return Edge_range(m_edges);
	}
	Edge_const_range edges() const {
		return Edge_const_range(m_edges);
	}
	Edge_iterator edges_begin() {
		return m_edges.begin();
	}
	Edge_iterator edges_end() {
		return m_edges.end();
	}
	Edge_const_iterator edges_begin() const {
		return m_edges.cbegin();
	}
	Edge_const_iterator edges_end() const {
		return m_edges.cend();
	}

	void clear(bool clear_init = false) {
				
		if constexpr (GraphTraits::historic) {
			m_history.clear();
		}

		for (Graph_vertex_index_listener* list : m_vertex_listeners) {
			list->clear();
		}
		m_vertex_listeners.clear();
		for (Vertex_handle v : m_vertices) {
			delete v;
		}
		m_vertices.clear();

		for (Graph_edge_index_listener* list : m_edge_listeners) {
			list->clear();
		}
		m_edge_listeners.clear();
		for (Edge_handle e : m_edges) {
			delete e;
		}
		m_edges.clear();

		if constexpr (GraphTraits::decomposed) {
			for (Graph_path_index_listener* list : m_path_listeners) {
				list->clear();
			}
			m_path_listeners.clear();

			for (Path_handle p : m_paths) {
				delete p;
			}
			m_paths.clear();
		}
		m_initialized = !m_initialized || clear_init;
	}

	History& history() requires GraphTraits::historic {
		return m_history;
	}

	inline bool can_perform_operation() {
		if constexpr (GraphTraits::historic) {
			return !m_history.can_redo();
		} else {
			return true;
		}
	}

	inline void start_operation_group() {
		if constexpr (GraphTraits::historic) {
			if (m_initialized)
				m_history.start_group();
		}
	}

	inline void end_operation_group() {
		if constexpr (GraphTraits::historic) {
			if (m_initialized)
				m_history.end_group();
		}
	}

	Vertex_handle add_vertex(Point_2 p) {

		Vertex_handle v = new Vertex(p, m_vertices.size());

		if (GraphTraits::historic && m_initialized) {
			auto op = std::make_unique<detail::AddVertex<Graph_2>>(*this, v);
			op->redo();
			m_history.add_operation(std::move(op));
		} else {
			detail::AddVertex<Graph_2>::add_vertex(*this, v);
		}

		assert(verify_graph_structure());
		return v;
	}
	void remove_vertex(Vertex_handle v) {

		if constexpr (GraphTraits::decomposed) {
			assert(!m_initialized); // not yet supported
		}

		if (GraphTraits::historic && m_initialized) {
			auto op = std::make_unique<detail::RemoveVertex<Graph_2>>(*this, v);
			op->redo();
			m_history.add_operation(std::move(op));
		} else {
			detail::RemoveVertex<Graph_2>::remove_vertex(*this, v);
			delete v;
		}

		assert(verify_graph_structure());
	}
	Edge_handle add_edge(Vertex_handle source,
	                     Vertex_handle target) requires std::same_as<Curve_2, Segment<Kernel>> {
		return add_edge(source, target, Curve_2(source->m_point, target->m_point));
	}

	Edge_handle add_edge(Vertex_handle source, Vertex_handle target, const Curve_2& curve) {

		if constexpr (GraphTraits::decomposed) {
			assert(!m_initialized); // not yet supported
		}

		Edge_handle e = new Edge(source, target, m_edges.size(), curve);

		if (GraphTraits::historic && m_initialized) {
			auto op = std::make_unique<detail::AddEdge<Graph_2>>(*this, e);
			op->redo();
			m_history.add_operation(std::move(op));
		} else {
			detail::AddEdge<Graph_2>::add_edge(*this, e);
		}

		assert(verify_graph_structure());
		return e;
	}
	void remove_edge(Edge_handle e) {

		if constexpr (GraphTraits::decomposed) {
			assert(!m_initialized); // not yet supported
		}

		if (GraphTraits::historic && m_initialized) {
			auto op = std::make_unique<detail::RemoveEdge<Graph_2>>(*this, e);
			op->redo();
			m_history.add_operation(std::move(op));
		} else {
			detail::RemoveEdge<Graph_2>::remove_edge(*this, e);
			delete e;
		}

		assert(verify_graph_structure());
	}

	void change_curve(Edge_handle e, Curve_representation new_rep) {

		if (GraphTraits::historic && m_initialized) {
			auto op = std::make_unique<detail::ChangeCurve<Graph_2>>(e, new_rep);
			op->redo();
			m_history.add_operation(std::move(op));
		} else {
			detail::ChangeCurve<Graph_2>::change_curve(e, new_rep);
		}
	}

	void change_curve(Edge_handle e, Curve_2 new_curve) {
		change_curve(e, Curve_traits::representation(new_curve));
	}

	void move_vertex(
	    Vertex_handle v, Point_2 p,
	    std::optional<std::function<typename CurveTraits::Curve_representation_2(Edge_const_handle)>>
	        e_to_rep = std::nullopt) {

		if (GraphTraits::historic && m_initialized) {
			m_history.start_group();
		}

		if constexpr (!std::same_as<std::monostate, typename CurveTraits::Curve_representation_2>) {
			for (auto eh : v->incident_edges()) {
				if (e_to_rep.has_value()) {
					change_curve(eh, (*e_to_rep)(eh));
				} else {
					auto new_representation = eh->representation();
					if (eh->source() == v) {
						Curve_traits::move_start(p, eh->target()->point(), new_representation);
					} else {
						Curve_traits::move_end(eh->source()->point(), p, new_representation);
					}
					change_curve(eh, std::move(new_representation));
				}
			}
		}

		if (GraphTraits::historic && m_initialized) {
			auto op = std::make_unique<detail::MoveVertex<Graph_2>>(v, p);
			op->redo();
			m_history.add_operation(std::move(op));
			m_history.end_group();
		} else {
			detail::MoveVertex<Graph_2>::move_vertex(v, p);
		}
	}

	Edge_handle
	split_vertex(Vertex_handle v, Point_2 p0,
	             Point_2 p1) requires std::same_as<Curve_2, Segment<Kernel>>&& GraphTraits::oriented {
		return split_vertex(v, Curve_2(v->prev()->m_point, p0), Curve_2(p0, p1),
		                    Curve_2(p1, v->next()->m_point));
	}
	/// Split a vertex into two. Or equivalently, replace two curves by three new ones.
	/// Connect the two new vertices with three new curves.
	Edge_handle split_vertex(Vertex_handle v, Curve_2 c0, Curve_2 c1,
	                         Curve_2 c2) requires GraphTraits::oriented {

		assert(m_initialized);
		assert(v->degree() == 2);
		assert(c0.target() == c1.source());
		assert(c1.target() == c2.source());

		Vertex_handle new_v = new Vertex(c1.target(), m_vertices.size());
		Edge_handle new_e = new Edge(v, new_v, m_edges.size(), c1);
		if constexpr (GraphTraits::decomposed) {
			new_e->m_path = v->outgoing()->m_path;
		}
		new_v->add_incident(new_e);
		new_v->add_incident(v->outgoing());

		if constexpr (GraphTraits::historic) {
			auto op = std::make_unique<detail::SplitVertex<Graph_2>>(
			    *this, v, new_v, new_e, CurveTraits::representation(c0),
			    CurveTraits::representation(c2));
			op->redo();
			m_history.add_operation(std::move(op));
		} else {
			detail::SplitVertex<Graph_2>::split_vertex(*this, v, new_v, new_e,
			                                           CurveTraits::representation(c0),
			                                           CurveTraits::representation(c2));
		}

		assert(verify_graph_structure());
		return new_e;
	}

	Vertex_handle collapse_edge(Edge_handle e, Point_2 newPoint)
	    requires std::same_as<Curve_2, Segment<Kernel>>&& GraphTraits::oriented {
		return collapse_edge(e, Curve_2(e->prev()->source()->m_point, newPoint),
		                     Curve_2(newPoint, e->next()->target()->m_point));
	}

	/// Collapse an edge, removing it together with the edge that comes before and after it.
	/// Replace it with a vertex connected by two new curves.
	/// This is equivalent to replacing three curves by two new ones.
	Vertex_handle collapse_edge(Edge_handle e, Curve_2 toNewPoint,
	                            Curve_2 fromNewPoint) requires GraphTraits::oriented {
		assert(m_initialized);
		assert(toNewPoint.target() == fromNewPoint.source());
		assert(e->source()->degree() == 2);
		assert(e->target()->degree() == 2);

		Vertex_handle result = e->source();

		if constexpr (GraphTraits::historic) {
			auto op = std::make_unique<detail::CollapseEdge<Graph_2>>(
			    *this, e, toNewPoint.target(), CurveTraits::representation(toNewPoint),
			    CurveTraits::representation(fromNewPoint));
			op->redo();
			m_history.add_operation(std::move(op));
		} else {
			detail::CollapseEdge<Graph_2>::collapse_edge(*this, e, toNewPoint.target(),
			                                             CurveTraits::representation(toNewPoint),
			                                             CurveTraits::representation(fromNewPoint));
		}

		assert(verify_graph_structure());
		return result;
	}

	/// Merges the provided edge wit the next, removing the given edge and its target. The remaining edge is returned.
	Edge_handle merge_with_next(
	    Edge_handle e) requires std::same_as<Curve_2, Segment<Kernel>>&& GraphTraits::oriented {
		return merge_vertex(e->target(), false);
	}

	/// Merges the provided edge wit the previous, removing the given edge and its source. The remaining edge is returned.
	Edge_handle merge_with_prev(
	    Edge_handle e) requires std::same_as<Curve_2, Segment<Kernel>>&& GraphTraits::oriented {
		return merge_vertex(e->source(), true);
	}

	/// Merges the two edges of a vertex into a single edge, erasing the vertex.
	/// If keepIncoming = true, then the vertex's outgoing edge is removed; otherwise, its incoming edge is removed.
	/// \pre Vertex v has degree 2 and its two neighbors are not neighbors
	Edge_handle merge_vertex(Vertex_handle v, bool keepIncoming = true)
	    requires std::same_as<Curve_2, Segment<Kernel>>&& GraphTraits::oriented {
		return merge_vertex(v, Curve_2(v->prev()->point(), v->next()->point()), keepIncoming);
	}

	/// Merges the provided edge wit the next, removing the given edge and its target. The remaining edge is returned.
	Edge_handle merge_with_next(Edge_handle e, Curve_2 newCurve,
	                            bool keepIncoming = true) requires GraphTraits::oriented {
		return merge_vertex(e->target(), newCurve, false);
	}

	/// Merges the provided edge wit the previous, removing the given edge and its source. The remaining edge is returned.
	Edge_handle merge_with_prev(Edge_handle e, Curve_2 newCurve,
	                            bool keepIncoming = true) requires GraphTraits::oriented {
		return merge_vertex(e->source(), newCurve, true);
	}

	/// Merges the two edges of a vertex into a single edge, erasing the vertex.
	/// If keepIncoming = true, then the vertex's outgoing edge is removed; otherwise, its incoming edge is removed.
	/// \pre Vertex v has degree 2 and its two neighbors are not neighbors
	Edge_handle merge_vertex(Vertex_handle v, Curve_2 newCurve,
	                         bool keepIncoming = true) requires GraphTraits::oriented {
		assert(m_initialized);
		assert(v->degree() == 2);
		assert(!v->prev()->is_neighbor_of(v->next()));

		Edge_handle result = keepIncoming ? v->incoming() : v->outgoing();

		if constexpr (GraphTraits::historic) {
			auto op = std::make_unique<detail::MergeVertex<Graph_2>>(
			    *this, v, keepIncoming, CurveTraits::representation(newCurve));
			op->redo();
			m_history.add_operation(std::move(op));
		} else {
			detail::MergeVertex<Graph_2>::merge_vertex(*this, v, keepIncoming,
			                                           CurveTraits::representation(newCurve));
		}

		assert(verify_graph_structure());
		return result;
	}

	/// Replace an edge by two edges, and returns the handle of the newly created vertex.
	/// The outgoing edge of this vertex is the provided edge.
	/// The incoming edge of this vertex is a new edge.
	Vertex_handle subdivide_edge_at_source(Edge_handle e, Point_2 newPoint)
	    requires std::same_as<Curve_2, Segment<Kernel>>&& GraphTraits::oriented {
		return subdivide_edge(e, newPoint, false);
	}

	/// Replace an edge by two edges, and returns the handle of the newly created vertex.
	/// The incoming edge of this vertex is the provided edge.
	/// The outgoing edge of this vertex is a new edge.
	Vertex_handle subdivide_edge_at_target(Edge_handle e, Point_2 newPoint)
	    requires std::same_as<Curve_2, Segment<Kernel>>&& GraphTraits::oriented {
		return subdivide_edge(e, newPoint, true);
	}

	/// Replace an edge by two edges, and returns the handle of the newly created vertex.
	/// The incoming edge of this vertex is the provided edge, if newOut = true, and the new edge otherwise.
	/// The outgoing edge of this vertex is a new edge, if newOut = true, and the provided edge otherwise.
	Vertex_handle subdivide_edge(Edge_handle e, Point_2 newPoint, bool newOut = true)
	    requires std::same_as<Curve_2, Segment<Kernel>>&& GraphTraits::oriented {
		return subdivide_edge(e, Curve_2(e->source()->m_point, newPoint),
		                      Curve_2(newPoint, e->target()->m_point), newOut);
	}

	/// Replace an edge by two edges, and returns the handle of the newly created vertex.
	/// The outgoing edge of this vertex is the provided edge.
	/// The incoming edge of this vertex is a new edge.
	Vertex_handle subdivide_edge_at_source(Edge_handle e, Curve_2 toNewPoint, Curve_2 fromNewPoint)
	    requires std::same_as<Curve_2, Segment<Kernel>>&& GraphTraits::oriented {
		return subdivide_edge(e, toNewPoint, fromNewPoint, false);
	}

	/// Replace an edge by two edges, and returns the handle of the newly created vertex.
	/// The incoming edge of this vertex is the provided edge.
	/// The outgoing edge of this vertex is a new edge.
	Vertex_handle subdivide_edge_at_target(Edge_handle e, Curve_2 toNewPoint, Curve_2 fromNewPoint)
	    requires std::same_as<Curve_2, Segment<Kernel>>&& GraphTraits::oriented {
		return subdivide_edge(e, toNewPoint, fromNewPoint, true);
	}

	/// Replace an edge by two edges, and returns the handle of the newly created vertex.
	/// The incoming edge of this vertex is the provided edge, if newOut = true, and the new edge otherwise.
	/// The outgoing edge of this vertex is a new edge, if newOut = true, and the provided edge otherwise.
	Vertex_handle subdivide_edge(Edge_handle e, Curve_2 toNewPoint, Curve_2 fromNewPoint,
	                             bool newOut = true) requires GraphTraits::oriented {
		assert(m_initialized);
		assert(toNewPoint.target() == fromNewPoint.source());

		Vertex_handle new_v = new Vertex(toNewPoint.target(), m_vertices.size());

		Edge_handle new_e;
		if (newOut) {
			new_e = new Edge(new_v, e->target(), m_edges.size(), fromNewPoint);
			new_v->add_incident(e);
			new_v->add_incident(new_e);
		} else {
			new_e = new Edge(e->source(), new_v, m_edges.size(), toNewPoint);
			new_v->add_incident(new_e);
			new_v->add_incident(e);
		}
		if constexpr (GraphTraits::decomposed) {
			new_e->m_path = e->m_path;
		}

		if constexpr (GraphTraits::historic) {
			auto op = std::make_unique<detail::SubdivideEdge<Graph_2>>(
			    *this, e, new_e, newOut,
			    CurveTraits::representation(newOut ? toNewPoint : fromNewPoint));
			op->redo();
			m_history.add_operation(std::move(op));
		} else {
			detail::SubdivideEdge<Graph_2>::subdivide_edge(
			    *this, e, new_e, newOut,
			    CurveTraits::representation(newOut ? toNewPoint : fromNewPoint));
		}

		assert(verify_graph_structure());
		return new_v;
	}

	CGAL::Bbox_2 bbox() const {
		if (m_vertices.empty()) {
			return {0, 0, 0, 0};
		}
		std::vector<Point<Inexact>> points;
		for (Vertex_handle v : m_vertices) {
			points.push_back(approximate(v->point()));
		}
		auto box = CGAL::bbox_2(points.begin(), points.end());
		if constexpr (!std::same_as<Curve_2, Segment<Kernel>>) {
			for (const auto& e : m_edges) {
				box = box + Curve_traits::bbox(e->source()->point(), e->target()->point(),
				                               e->representation());
			}
		}
		return box;
	}

	CGAL::Iso_rectangle_2<Kernel> bounding_rectangle() const {

		typename Kernel::FT left = 0, right = 0, bottom = 0, top = 0;

		bool first = true;
		for (Vertex_handle v : m_vertices) {
			const Point<Kernel> pt = v->m_point;

			if (first) {
				left = right = pt.x();
				top = bottom = pt.y();
				first = false;
			} else {
				if (pt.x() < left) {
					left = pt.x();
				} else if (pt.x() > right) {
					right = pt.x();
				}

				if (pt.y() < bottom) {
					bottom = pt.y();
				} else if (pt.y() > top) {
					top = pt.y();
				}
			}
		}
		//  TODO: include edges

		return CGAL::Iso_rectangle_2<Kernel>(left, bottom, right, top);
	}
};

template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits, GraphTraits_2 GraphTraits>
class Graph_2_vertex {
	friend class Graph_2<VertexData, EdgeData, CurveTraits, GraphTraits>;
	friend class Graph_2_edge<VertexData, EdgeData, CurveTraits, GraphTraits>;

	template <class InGraph, class OutGraph>
	friend void
	graph_2_copy(InGraph& input, OutGraph& output,
	             const std::function<typename OutGraph::Curve_traits::Curve_representation_2(
	                 const typename InGraph::Curve_traits::Curve_representation_2&)>& conversion,
	             bool requireResorting);

	// basic operations
	template <class G> friend class detail::AddVertex;
	template <class G> friend class detail::RemoveVertex;
	template <class G> friend class detail::AddEdge;
	template <class G> friend class detail::RemoveEdge;
	// geometric operations
	template <class G> friend class detail::ChangeCurve;
	template <class G> friend class detail::MoveVertex;
	// oriented operations
	template <class G> friend class detail::MergeVertex;
	template <class G> friend class detail::SplitVertex;
	template <class G> friend class detail::CollapseEdge;
	template <class G> friend class detail::SubdivideEdge;

  public:
	using Graph = Graph_2<VertexData, EdgeData, CurveTraits, GraphTraits>;
	using Vertex = Graph::Vertex;
	using Vertex_handle = Graph::Vertex_handle;
	using Vertex_const_handle = Graph::Vertex_const_handle;
	using Edge = Graph::Edge;
	using Edge_handle = Graph::Edge_handle;
	using Edge_const_handle = Graph::Edge_const_handle;
	using Path = Graph::Path;
	using Path_handle = Graph::Path_handle;
	using Path_const_handle = Graph::Path_const_handle;

	using Vertex_data = VertexData;
	using Edge_data = EdgeData;

	using Curve_traits = CurveTraits;
	using Kernel = CurveTraits::Kernel;
	using Point_2 = CurveTraits::Point_2;
	using Curve_2 = CurveTraits::Curve_2;
	using Curve_representation = CurveTraits::Curve_representation_2;

  private:
	using Edge_container = std::vector<Edge_handle>;
	Point_2 m_point;
	Edge_container m_incident = Edge_container();
	Vertex_data m_data;
	size_t m_index;

	void remove_incident(const size_t index) {
		m_incident.erase(std::next(m_incident.begin(), index));
	}

	void remove_incident(Edge_handle e) {
		remove_incident(find_incident_index(e));
	}

	void insert_incident(Edge_handle e, const size_t index) {
		m_incident.insert(std::next(m_incident.begin(), index), e);
	}

	void add_incident(Edge_handle e) {
		m_incident.push_back(e);
	}

	Graph_2_vertex(const Point_2 point, const size_t index)
	    : m_point(std::move(point)), m_index(std::move(index)) {}
	Graph_2_vertex(const Point_2 point, const size_t index, const Vertex_data data)
	    : m_point(std::move(point)), m_index(std::move(index)), m_data(std::move(data)) {}

  public:
	size_t graph_index() const {
		return m_index;
	}
	const Point_2& point() const {
		return m_point;
	}
	size_t degree() const {
		return m_incident.size();
	}
	Edge_handle incident_edge(size_t i) {
		return m_incident[i];
	}
	Edge_const_handle incident_edge(const size_t i) const {
		return m_incident[i];
	}
	Vertex_handle neighbor(size_t i) {
		return m_incident[i]->other(this);
	}
	Vertex_const_handle neighbor(const size_t i) const {
		return m_incident[i]->other(this);
	}
	size_t find_incident_index(Edge_const_handle e) const {
		assert(e->m_source == this || e->m_target == this);
		const size_t deg = degree();
		for (size_t index = 0; index < deg; ++index) {
			if (m_incident[index] == e) {
				return index;
			}
		}
		assert(false);
		return -1; // NB: this will give some positive number, as size_t is unsigned
	}
	Edge_handle neighboring_incident_edge(Edge_handle e, const bool ccw) requires GraphTraits::sorted {
		return neighboring_incident_edge(find_incident_index(e), ccw);
	}
	Edge_handle neighboring_incident_edge(size_t i, const bool ccw) requires GraphTraits::sorted {
		if (ccw) {
			return m_incident[(i + 1) % m_incident.size()];
		} else {
			return m_incident[(i + m_incident.size() - 1) % m_incident.size()];
		}
	}
	Edge_handle find_edge_to(Vertex_const_handle v) const {
		for (Edge_handle e : m_incident) {
			if (e->other(this) == v) {
				return e;
			}
		}
		return nullptr;
	}
	bool is_neighbor_of(Vertex_const_handle v) const {
		for (Edge_const_handle e : m_incident) {
			if (e->other(this) == v) {
				return true;
			}
		}
		return false;
	}
	Edge_handle incoming() requires GraphTraits::oriented {
		assert(degree() == 2);
		return m_incident[0];
	}
	Edge_const_handle incoming() const requires GraphTraits::oriented {
		assert(degree() == 2);
		return m_incident[0];
	}
	Edge_handle outgoing() requires GraphTraits::oriented {
		assert(degree() == 2);
		return m_incident[1];
	}
	Edge_const_handle outgoing() const requires GraphTraits::oriented {
		assert(degree() == 2);
		return m_incident[1];
	}
	Vertex_handle prev() requires GraphTraits::oriented {
		return incoming()->source();
	}
	Vertex_const_handle prev() const requires GraphTraits::oriented {
		return incoming()->source();
	}
	Vertex_handle next() requires GraphTraits::oriented {
		return outgoing()->target();
	}
	Vertex_const_handle next() const requires GraphTraits::oriented {
		return outgoing()->target();
	}
	Edge_container::iterator incident_edges_begin() {
		return m_incident.begin();
	}
	Edge_container::iterator incident_edges_end() {
		return m_incident.end();
	}
	Edge_container::const_iterator incident_edges_begin() const {
		return m_incident.cbegin();
	}
	Edge_container::const_iterator incident_edges_end() const {
		return m_incident.cend();
	}
	Vertex_data& data() {
		return m_data;
	}
	const Vertex_data& data() const {
		return m_data;
	}

	class Incident_edges_range {
		friend class Graph_2_vertex<VertexData, EdgeData, CurveTraits, GraphTraits>;

	  private:
		Edge_container& m_container;

		Incident_edges_range(Edge_container& container) : m_container(container) {}

	  public:
		Edge_container::iterator begin() const {
			return m_container.begin();
		}

		Edge_container::iterator end() const {
			return m_container.end();
		}
	};

	class Incident_edges_const_range {
		friend class Graph_2_vertex<VertexData, EdgeData, CurveTraits, GraphTraits>;

	  private:
		const Edge_container& m_container;

		Incident_edges_const_range(const Edge_container& container) : m_container(container) {}

	  public:
		Edge_container::const_iterator begin() const {
			return m_container.cbegin();
		}

		Edge_container::const_iterator end() const {
			return m_container.cend();
		}
	};

	Incident_edges_range incident_edges() {
		return Incident_edges_range(m_incident);
	}
	Incident_edges_const_range incident_edges() const {
		return Incident_edges_const_range(m_incident);
	}

	friend std::ostream&
	operator<<(std::ostream& os,
	           const Graph_2_vertex<VertexData, EdgeData, CurveTraits, GraphTraits>& c) {
		return os << "V[" << c.m_index << " @ " << c.m_point << " ; d= " << c.m_incident.size()
		          << "]";
	}
};

template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits, GraphTraits_2 GraphTraits>
class Graph_2_edge {
	friend class Graph_2<VertexData, EdgeData, CurveTraits, GraphTraits>;
	friend class Graph_2_vertex<VertexData, EdgeData, CurveTraits, GraphTraits>;
	friend class Graph_2_path<VertexData, EdgeData, CurveTraits, GraphTraits>;

	template <class InGraph, class OutGraph>
	friend void
	graph_2_copy(InGraph& input, OutGraph& output,
	             const std::function<typename OutGraph::Curve_traits::Curve_representation_2(
	                 const typename InGraph::Curve_traits::Curve_representation_2&)>& conversion,
	             bool requireResorting);

	// basic operations
	template <class G> friend class detail::AddVertex;
	template <class G> friend class detail::RemoveVertex;
	template <class G> friend class detail::AddEdge;
	template <class G> friend class detail::RemoveEdge;
	// geometric operations
	template <class G> friend class detail::ChangeCurve;
	template <class G> friend class detail::MoveVertex;
	// oriented operations
	template <class G> friend class detail::MergeVertex;
	template <class G> friend class detail::SplitVertex;
	template <class G> friend class detail::CollapseEdge;
	template <class G> friend class detail::SubdivideEdge;

  public:
	using Graph = Graph_2<VertexData, EdgeData, CurveTraits, GraphTraits>;
	using Vertex = Graph::Vertex;
	using Vertex_handle = Graph::Vertex_handle;
	using Vertex_const_handle = Graph::Vertex_const_handle;
	using Edge = Graph::Edge;
	using Edge_handle = Graph::Edge_handle;
	using Edge_const_handle = Graph::Edge_const_handle;
	using Path = Graph::Path;
	using Path_handle = Graph::Path_handle;
	using Path_const_handle = Graph::Path_const_handle;

	using Vertex_data = VertexData;
	using Edge_data = EdgeData;

	using Curve_traits = CurveTraits;
	using Kernel = CurveTraits::Kernel;
	using Point_2 = CurveTraits::Point_2;
	using Curve_2 = CurveTraits::Curve_2;
	using Curve_representation = CurveTraits::Curve_representation_2;

  private:
	Vertex_handle m_source;
	Vertex_handle m_target;
	Curve_representation m_representation;
	size_t m_index;
	Edge_data m_data;
	Path_handle m_path = nullptr;

	Graph_2_edge(Vertex_handle source, Vertex_handle target,
	             const size_t index) requires std::same_as<Curve_representation, std::monostate>
	    : m_source(std::move(source)), m_target(std::move(target)), m_index(std::move(index)) {}

	Graph_2_edge(Vertex_handle source, Vertex_handle target, const size_t index, const Curve_2 curve)
	    : m_source(std::move(source)), m_target(std::move(target)), m_index(std::move(index)),
	      m_representation(CurveTraits::representation(curve)) {}

	Graph_2_edge(Vertex_handle source, Vertex_handle target, const size_t index,
	             const Curve_representation rep)
	    : m_source(std::move(source)), m_target(std::move(target)), m_index(std::move(index)),
	      m_representation(std::move(rep)) {}

	Graph_2_edge(Vertex_handle source, Vertex_handle target, const size_t index,
	             const Curve_representation rep, Edge_data data)
	    : m_source(std::move(source)), m_target(std::move(target)), m_index(std::move(index)),
	      m_representation(std::move(rep)), m_data(std::move(data)) {}

  public:
	size_t graph_index() const {
		return m_index;
	}
	Vertex_handle source() {
		return m_source;
	}
	Vertex_handle target() {
		return m_target;
	}
	Vertex_const_handle source() const {
		return m_source;
	}
	Vertex_const_handle target() const {
		return m_target;
	}
	Path_handle path() requires GraphTraits::decomposed {
		return m_path;
	}
	Path_const_handle path() const requires GraphTraits::decomposed {
		return m_path;
	}
	Edge_data& data() {
		return m_data;
	}
	const Edge_data& data() const {
		return m_data;
	}

	size_t find_source_incident_index() const {
		return m_source->find_incident_index(this);
	}
	size_t find_target_incident_index() const {
		return m_target->find_incident_index(this);
	}
	size_t find_other_incident_index(Vertex_const_handle v) const {
		return other(v)->find_incident_index(this);
	}

	Vertex_handle common_vertex(Edge_handle other) {
		if (other->m_source == m_source || other->m_target == m_source) {
			return m_source;
		} else if (other->m_source == m_target || other->m_target == m_target) {
			return m_target;
		} else {
			return nullptr;
		}
	}
	Vertex_handle other(Vertex_handle v) {
		assert(v == m_source || v == m_target);
		return v == m_source ? m_target : m_source;
	}
	Vertex_const_handle other(Vertex_const_handle v) const {
		assert(v == m_source || v == m_target);
		return v == m_source ? m_target : m_source;
	}
	const Curve_representation& representation() const {
		return m_representation;
	}
	Curve_2 curve() const {
		return CurveTraits::curve(m_source->m_point, m_target->m_point, m_representation);
	}
	void reverse() {
		std::swap(m_source, m_target);
		CurveTraits::reverse_representation(m_source->m_point, m_target->m_point, m_representation);
	}
	Edge_handle prev() requires GraphTraits::oriented {
		return m_source->incoming();
	}
	Edge_const_handle prev() const requires GraphTraits::oriented {
		return m_source->incoming();
	}
	Edge_handle next() requires GraphTraits::oriented {
		return m_target->outgoing();
	}
	Edge_const_handle next() const requires GraphTraits::oriented {
		return m_target->outgoing();
	}
	Edge_handle next(const bool ccw) requires GraphTraits::sorted {
		return m_target->neighboring_incident_edge(this, ccw);
	}
	Edge_handle prev(bool ccw) requires GraphTraits::sorted {
		return m_source->neighboring_incident_edge(this, ccw);
	}

	CGAL::Direction_2<Kernel> source_direction() const {
		return CurveTraits::source_direction(m_source->m_point, m_target->m_point, m_representation);
	}

	CGAL::Direction_2<Kernel> target_direction() const {
		return CurveTraits::target_direction(m_source->m_point, m_target->m_point, m_representation);
	}

	CGAL::Direction_2<Kernel> direction_at(Vertex_const_handle v) const {
		assert(v == m_source || v == m_target);
		if (v == m_source) {
			return source_direction();
		} else {
			return target_direction();
		}
	}

	friend std::ostream&
	operator<<(std::ostream& os,
	           const Graph_2_edge<VertexData, EdgeData, CurveTraits, GraphTraits>& c) {
		return os << "E[" << c.m_index << " : " << *(c.m_source) << " -> " << *(c.m_target) << "]";
	}
};

template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits, GraphTraits_2 GraphTraits>
class Graph_2_path {
	friend class Graph_2<VertexData, EdgeData, CurveTraits, GraphTraits>;
	friend class Graph_2_vertex<VertexData, EdgeData, CurveTraits, GraphTraits>;
	friend class Graph_2_edge<VertexData, EdgeData, CurveTraits, GraphTraits>;

	template <class InGraph, class OutGraph>
	friend void
	graph_2_copy(InGraph& input, OutGraph& output,
	             const std::function<typename OutGraph::Curve_traits::Curve_representation_2(
	                 const typename InGraph::Curve_traits::Curve_representation_2&)>& conversion,
	             bool requireResorting);

	// basic operations
	template <class G> friend class detail::AddVertex;
	template <class G> friend class detail::RemoveVertex;
	template <class G> friend class detail::AddEdge;
	template <class G> friend class detail::RemoveEdge;
	// geometric operations
	template <class G> friend class detail::ChangeCurve;
	template <class G> friend class detail::MoveVertex;
	// oriented operations
	template <class G> friend class detail::MergeVertex;
	template <class G> friend class detail::SplitVertex;
	template <class G> friend class detail::CollapseEdge;
	template <class G> friend class detail::SubdivideEdge;

  public:
	using Graph = Graph_2<VertexData, EdgeData, CurveTraits, GraphTraits>;
	using Vertex = Graph::Vertex;
	using Vertex_handle = Graph::Vertex_handle;
	using Vertex_const_handle = Graph::Vertex_const_handle;
	using Edge = Graph::Edge;
	using Edge_handle = Graph::Edge_handle;
	using Edge_const_handle = Graph::Edge_const_handle;
	using Path = Graph::Path;
	using Path_handle = Graph::Path_handle;
	using Path_const_handle = Graph::Path_const_handle;

	using Vertex_data = VertexData;
	using Edge_data = EdgeData;
	using Path_data = GraphTraits::PathData;

	using Curve_traits = CurveTraits;
	using Kernel = CurveTraits::Kernel;
	using Point_2 = CurveTraits::Point_2;
	using Curve_2 = CurveTraits::Curve_2;
	using Curve_representation = CurveTraits::Curve_representation_2;

  private:
	Edge_handle m_start, m_end;
	bool m_cyclic;
	size_t m_index;
	Path_data m_data;

	Graph_2_path(Edge_handle src, size_t index)
	    : m_start(src), m_end(src), m_cyclic(false), m_index(index) {}

	Graph_2_path(Edge_handle src, Edge_handle end, bool cyclic, size_t index)
	    : m_start(src), m_end(end), m_cyclic(cyclic), m_index(index) {}

	Graph_2_path(Edge_handle start, Edge_handle end, bool cyclic, size_t index, Path_data data)
	    : m_start(start), m_end(end), m_cyclic(cyclic), m_index(index), m_data(data) {}

  public:
	Edge_handle first_edge() {
		return m_start;
	}
	Edge_handle last_edge() {
		return m_end;
	}
	size_t graph_index() const {
		return m_index;
	}
	bool cyclic() const {
		return m_cyclic;
	}
	Path_data& data() {
		return m_data;
	}
	const Path_data& data() const {
		return m_data;
	}

	friend std::ostream&
	operator<<(std::ostream& os,
	           const Graph_2_path<VertexData, EdgeData, CurveTraits, GraphTraits>& c) {
		return os << "P[" << c.m_index << " : " << (c.m_cyclic ? "cyclic, " : "acyclic, ")
		          << *(c.m_start) << " -> " << *(c.m_end) << "]";
	}
};

template <class InGraph, class OutGraph>
requires std::is_same<typename InGraph::Curve_traits::Curve_representation_2,
                      typename OutGraph::Curve_traits::Curve_representation_2>::value void
graph_2_copy(InGraph& input, OutGraph& output) {
	using CurveRep = InGraph::Curve_traits::Curve_representation_2;
	graph_2_copy(input, output, [](const CurveRep& cr) { return cr; }, false);
}

/// General graph-copy function. A conversion function between the curve representations must be provided.
///
/// - The output graph will be cleared, and deinitialized
/// - All vertices and edges are copied
/// - The output graph will be initialized if the input graph is initialized.
/// - Indices of vertices, edges and paths (if both graphs have the decomposed GraphTrait) will match
/// - Any history of the input graph will NOT be copied
/// - Any associated data with vertices, edges or paths are NOT copied
///
/// Note that changes made to the curves may require re-sorting the edges. If it is ensured, by construction,
/// that order around vertices is not changed, setting the last parameter to false may give a performance gain.
template <class InGraph, class OutGraph>
void graph_2_copy(InGraph& input, OutGraph& output,
                  const std::function<typename OutGraph::Curve_traits::Curve_representation_2(
                      const typename InGraph::Curve_traits::Curve_representation_2&)>& conversion,
                  bool requireResorting = true) {

	using InVertex = InGraph::Vertex_handle;
	using InEdge = InGraph::Edge_handle;
	using OutVertex = OutGraph::Vertex_handle;
	using OutEdge = OutGraph::Edge_handle;

	output.clear(true);

	output.m_vertices.resize(input.m_vertices.size());
	for (InVertex v : input.m_vertices) {
		output.m_vertices[v->m_index] =
		    new OutGraph::Vertex(convert_kernel<typename OutGraph::Kernel>(v->m_point), v->m_index);
	}
	output.m_edges.resize(input.m_edges.size());
	for (InEdge e : input.m_edges) {
		output.m_edges[e->m_index] = new OutGraph::Edge(output.m_vertices[e->m_source->m_index],
		                                                output.m_vertices[e->m_target->m_index],
		                                                e->m_index, conversion(e->m_representation));
	}
	for (OutVertex v : output.m_vertices) {
		InVertex in_v = input.m_vertices[v->m_index];
		const size_t d = in_v->degree();
		v->m_incident.resize(d);
		for (size_t i = 0; i < d; ++i) {
			v->m_incident[i] = output.m_edges[in_v->m_incident[i]->m_index];
		}
	}

	if (input.is_initialized()) {
		output.m_initialized = true;

		if constexpr (OutGraph::Graph_traits::oriented) {
			if constexpr (InGraph::Graph_traits::oriented) {
				// nothing to do, oriented is copied by design
			} else {
				output.ensure_oriented();
			}
		}

		if constexpr (OutGraph::Graph_traits::sorted) {
			if constexpr (InGraph::Graph_traits::sorted) {
				// nothing to do, sorted is copied by design, unless the geometry has changed in the conversion
				if (requireResorting) {
					output.ensure_sorted();
				}
			} else {
				output.ensure_sorted();
			}
		}

		if constexpr (OutGraph::Graph_traits::decomposed) {
			if constexpr (InGraph::Graph_traits::decomposed) {

				using InPath = InGraph::Path_handle;
				using OutPath = OutGraph::Path_handle;

				output.m_paths.resize(input.m_paths.size());
				for (const InPath p : input.m_paths) {
					OutPath op = new OutGraph::Path(output.m_edges[p->m_start->m_index],
					                                output.m_edges[p->m_end->m_index], p->m_cyclic,
					                                p->m_index);
					output.m_paths[p->m_index] = op;
				}

				for (OutEdge e : output.m_edges) {
					e->m_path = output.m_paths[input.m_edges[e->m_index]->m_path->m_index];
				}
			} else {
				output.ensure_decomposed();
			}
		}

		assert(output.is_initialized());
	}
}

} // namespace CGAL::internal
