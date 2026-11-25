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

namespace CGAL::internal {
class Graph_map_base;

class Operation {

  public:
	virtual void undo() = 0;
	virtual void redo() = 0;
	virtual void forget_past() {}
	virtual void forget_future() {}
};

class OperationGroup {
  private:
	std::vector<std::unique_ptr<Operation>> ops;

  public:
	void add_operation(std::unique_ptr<Operation> op) {
		ops.push_back(std::move(op));
	}

	void undo() {
		for (auto it = ops.rbegin(); it != ops.rend(); ++it) {
			(*it)->undo();
		}
	}

	void redo() {
		for (auto& op : ops) {
			op->redo();
		}
	}

	void forget_past() {
		for (auto& op : ops) {
			op->forget_past();
		}
	}
	void forget_future() {
		for (auto& op : ops) {
			op->forget_future();
		}
	}

	size_t number_of_operations() {
		return ops.size();
	}
};

struct NoGraphHistory {
	void start_group() {
		assert(false); // should not be called
	}
	void end_group() {
		assert(false); // should not be called
	}
	void add_operation(std::unique_ptr<Operation> op) {
		assert(false); // should not be called
	}

	inline bool can_undo() const {
		assert(false); // should not be called
		return false;
	}
	inline void undo() {
		assert(false); // should not be called
	}
	void undo_all() {
		assert(false); // should not be called
	}
	inline bool can_redo() const {
		assert(false); // should not be called
		return false;
	}
	inline void redo() {
		assert(false); // should not be called
	}
	void redo_all() {
		assert(false); // should not be called
	}
	void clear() {
		assert(false); // should not be called
	}
	void forget_past() {
		assert(false); // should not be called
	}
	void forget_future() {
		assert(false); // should not be called
	}
};

class GraphHistory {
  private:
	using GroupContainer = std::list<OperationGroup>;
	using GroupHandle = GroupContainer::iterator;

	GroupContainer m_groups;
	GroupHandle m_curr = m_groups.end(); // points to the next redo-operation
	int m_build = 0;

  public:
	~GraphHistory() {
		clear();
	}

	void start_group() {
		if (m_build == 0) {
			assert(m_curr == m_groups.end());
			m_groups.emplace_back();
			--m_curr;
		}
		++m_build;
	}
	void end_group() {
		assert(m_build > 0);
		--m_build;

		if (m_build == 0) {
			if (m_curr->number_of_operations() == 0) {
				m_groups.pop_back();
				m_curr = m_groups.end();
			} else {
				++m_curr;
			}
		}
	}
	void add_operation(std::unique_ptr<Operation> op) {
		if (m_build == 0) {
			start_group();
			m_curr->add_operation(std::move(op));
			end_group();
		} else {
			m_curr->add_operation(std::move(op));
		}
	}

	inline bool can_undo() const {
		assert(m_build == 0);
		return m_curr != m_groups.begin();
	}
	inline void undo() {
		--m_curr;
		m_curr->undo();
	}
	void undo_all() {
		while (can_undo()) {
			undo();
		}
	}
	inline bool can_redo() const {
		return m_curr != m_groups.end();
	}
	inline void redo() {
		m_curr->redo();
		++m_curr;
	}
	void redo_all() {
		while (can_redo()) {
			redo();
		}
	}
	void clear() {
		forget_past();
		forget_future();
	}
	void forget_past() {
		while (m_curr != m_groups.begin()) {
			m_groups.front().forget_past();
			m_groups.pop_front();
		}
	}
	void forget_future() {
		if (m_curr != m_groups.end()) {
			while (true) {
				bool deleting_curr = m_curr == (--m_groups.end());
				m_groups.back().forget_future();
				m_groups.pop_back();
				if (deleting_curr)
					break;
			}

			m_curr = m_groups.end();
		}
	}
};

namespace detail {

template <class G> class AddVertex : public Operation {
  private:
	using Vertex_handle = G::Vertex_handle;

	G& m_graph;
	Vertex_handle m_vertex;

  public:
	AddVertex(G& graph, Vertex_handle v) : m_graph(graph), m_vertex(v) {}

	void forget_future() override {
		delete m_vertex;
	}

	static void add_vertex(G& g, Vertex_handle v) {
		g.insert_vertex_into_container(v);
		// NB: cannot violate traits
	}

	void undo() override {
		m_graph.remove_vertex_from_container(m_vertex);
		// NB: cannot violate traits
	}
	void redo() override {
		add_vertex(m_graph, m_vertex);
	}
};

template <class G> class RemoveVertex : public Operation {
  private:
	using Vertex_handle = G::Vertex_handle;
	using Edge_handle = G::Edge_handle;

	G& m_graph;
	Vertex_handle m_vertex;
	std::vector<size_t> m_indices;

  public:
	RemoveVertex(G& graph, Vertex_handle v) : m_graph(graph), m_vertex(v) {
		for (Edge_handle e : m_vertex->m_incident) {
			m_indices.push_back(e->find_other_incident_index(m_vertex));
		}
	}

	void forget_past() {
		delete m_vertex;
	}

	static void remove_vertex(G& g, Vertex_handle v) {
		for (Edge_handle e : v->m_incident) {
			e->other(v)->remove_incident(e);
			g.remove_edge_from_container(e);
		}
		g.remove_vertex_from_container(v);
		// NB: cannot violate sorted,oriented traits
		// TODO: maintain decomposed trait
	}

	void undo() override {
		for (size_t i = m_vertex->m_incident.size() - 1; i >= 0; --i) {
			Edge_handle e = m_vertex->m_incident[i];
			e->other(m_vertex)->insert_incident(e, i);
			m_graph.insert_edge_into_container(e);
		}
		m_graph.insert_vertex_into_container(m_vertex);
		// NB: cannot violate sorted,oriented traits
		// TODO: maintain decomposed trait
	}
	void redo() override {
		for (Edge_handle e : m_vertex->m_incident) {
			e->other(m_vertex)->remove_incident(e);
		}
		m_graph.remove_vertex_from_container(m_vertex);
		// NB: cannot violate sorted,oriented traits
		// TODO: maintain decomposed trait
	}
};

template <class G> class AddEdge : public Operation {
  private:
	using Edge_handle = G::Edge_handle;

	G& m_graph;
	Edge_handle m_edge;

  public:
	AddEdge(G& graph, Edge_handle e) : m_graph(graph), m_edge(e) {}

	void forget_future() override {
		delete m_edge;
	}

	static void add_edge(G& graph, Edge_handle e) {

		e->m_source->add_incident(e);
		e->m_target->add_incident(e);
		graph.insert_edge_into_container(e);

		if (graph.m_initialized) {
			// TODO: speed up
			graph.ensure_traits(e->m_source);
			graph.ensure_traits(e->m_target);

			// TODO: support decomposed
		}
	}

	void undo() override {
		m_graph.remove_edge_from_container(m_edge);

		m_edge->m_source->remove_incident(m_edge);
		m_edge->m_target->remove_incident(m_edge);

		// TODO: speed up
		m_graph.ensure_traits(m_edge->m_source);
		m_graph.ensure_traits(m_edge->m_target);

		// TODO: support decomposed
	}
	void redo() override {
		add_edge(m_graph, m_edge);
	}
};

template <class G> class RemoveEdge : public Operation {
  private:
	using Edge_handle = G::Edge_handle;

	G& m_graph;
	Edge_handle m_edge;

  public:
	RemoveEdge(G& graph, Edge_handle e) : m_graph(graph), m_edge(e) {}

	void forget_past() override {
		delete m_edge;
	}

	static void remove_edge(G& g, Edge_handle e) {
		g.remove_edge_from_container(e);

		e->m_source->remove_incident(e);
		e->m_target->remove_incident(e);

		// TODO: more efficient?
		g.ensure_traits(e->m_source);
		g.ensure_traits(e->m_target);

		// TODO: support decomposed
	}

	void undo() override {
		m_graph.insert_edge_into_container(m_edge);

		m_edge->m_source->m_incident.push_back(m_edge);
		m_edge->m_target->m_incident.push_back(m_edge);

		// TODO: more efficient?
		m_graph.ensure_traits(m_edge->source());
		m_graph.ensure_traits(m_edge->target());

		// TODO: support decomposed
	}
	void redo() override {
		remove_edge(m_graph, m_edge);
	}
};

template <class G> class ChangeCurve : public Operation {
  private:
	using Edge_handle = G::Edge_handle;
	using Curve_representation = typename G::Curve_traits::Curve_representation_2;

	Edge_handle m_edge;
	Curve_representation m_representation;

  public:
	ChangeCurve(Edge_handle e, Curve_representation new_representation)
	    : m_edge(e), m_representation(std::move(new_representation)) {}

	static void change_curve(Edge_handle e, Curve_representation new_rep) {
		e->m_representation = std::move(new_rep);
		// TODO: maintain sorted trait
	}

	void undo() override {
		std::swap(m_edge->m_representation, m_representation);
		// TODO: maintain sorted trait
	}

	void redo() override {
		std::swap(m_edge->m_representation, m_representation);
		// TODO: maintain sorted trait
	}
};

template <class G> class MoveVertex : public Operation {
  private:
	using Vertex_handle = G::Vertex_handle;
	using Point = G::Point_2;

	Vertex_handle m_vertex;
	Point m_point;

  public:
	MoveVertex(Vertex_handle v, Point p) : m_vertex(v), m_point(std::move(p)) {}

	static void move_vertex(Vertex_handle v, Point p) {
		v->m_point = std::move(p);
		// TODO: maintain sorted trait
	}

	void undo() override {
		std::swap(m_point, m_vertex->m_point);
		// TODO: maintain sorted trait
	}
	void redo() override {
		std::swap(m_point, m_vertex->m_point);
		// TODO: maintain sorted trait
	}
};

template <class G> class MergeVertex : public Operation {
  private:
	using Vertex_handle = G::Vertex_handle;
	using Edge_handle = G::Edge_handle;
	using Path_handle = G::Path_handle;
	using Point = G::Point_2;
	using Curve_representation = typename G::Curve_traits::Curve_representation_2;

	G& m_graph;
	Vertex_handle m_vertex;
	Curve_representation m_curve_rep;
	size_t m_index;
	bool m_keep_incoming;

	static void merge_vertex_no_curve(G& g, Vertex_handle v, const size_t index,
	                                  const bool keepIncoming) {

		Edge_handle incoming = v->incoming();
		Edge_handle outgoing = v->outgoing();

		// remove v and v.outgoing/incoming from graph
		g.remove_vertex_from_container(v);

		if (keepIncoming) {
			g.remove_edge_from_container(outgoing);

			// redirect v.incoming to v.next
			Vertex_handle target = outgoing->target();
			target->m_incident[index] = incoming;
			incoming->m_target = target;

			// update paths
			if constexpr (G::Graph_traits::decomposed) {
				Path_handle p = outgoing->m_path;
				if (p->m_start == outgoing) {
					assert(p->m_cyclic);
					p->m_start = outgoing->next();
				} else if (p->m_end == outgoing) {
					p->m_end = incoming;
				}
			}

		} else {
			g.remove_edge_from_container(incoming);

			// redirect v.outgoing to v.prev
			Vertex_handle source = incoming->source();
			source->m_incident[index] = outgoing;
			outgoing->m_source = source;

			// update paths
			if constexpr (G::Graph_traits::decomposed) {
				Path_handle p = incoming->m_path;
				if (p->m_end == incoming) {
					assert(p->m_cyclic);
					p->m_end = incoming->prev();
				} else if (p->m_start == incoming) {
					p->m_start = outgoing;
				}
			}
		}
	}

	static void maintain_traits(G& g, Vertex_handle v) {

		// TODO: maintain sorted trait
	}

  public:
	MergeVertex(G& g, Vertex_handle v, bool keepIncoming, Curve_representation curve_rep)
	    : m_graph(g), m_vertex(std::move(v)), m_curve_rep(std::move(curve_rep)),
	      m_keep_incoming(keepIncoming) {
		m_index = keepIncoming ? v->outgoing()->find_target_incident_index()
		                       : v->incoming()->find_source_incident_index();
	}

	void forget_past() override {
		if (m_keep_incoming)
			delete m_vertex->outgoing();
		else
			delete m_vertex->incoming();
		delete m_vertex;
	}

	static void merge_vertex(G& g, Vertex_handle v, const bool keepIncoming,
	                         Curve_representation curve_rep) {

		size_t index = keepIncoming ? v->outgoing()->find_target_incident_index()
		                            : v->incoming()->find_source_incident_index();
		merge_vertex_no_curve(g, v, index, keepIncoming);
		if (keepIncoming)
			v->incoming()->m_representation = curve_rep;
		else
			v->outgoing()->m_representation = curve_rep;
		maintain_traits(g, v);
	}

	void undo() override {

		Edge_handle incoming = m_vertex->incoming();
		Edge_handle outgoing = m_vertex->outgoing();

		if (m_keep_incoming) {

			// update curve
			std::swap(incoming->m_representation, m_curve_rep);

			// redirect v.incoming to v
			incoming->m_target = m_vertex;

			// add v and v.outgoing to graph
			m_graph.insert_vertex_into_container(m_vertex);
			m_graph.insert_edge_into_container(outgoing);

			outgoing->target()->m_incident[m_index] = outgoing;

			// TODO: maintain sorted trait

			if constexpr (G::Graph_traits::decomposed) {
				Path_handle p = outgoing->m_path;
				if (p->m_cyclic && p->m_start == outgoing->next()) {
					p->m_start = outgoing;
				} else if (p->m_end == incoming) {
					p->m_end = outgoing;
				}
			}
		} else {

			// update curve
			std::swap(outgoing->m_representation, m_curve_rep);

			// redirect v.outgoing to v
			outgoing->m_source = m_vertex;

			// add v and v.incoming to graph
			m_graph.insert_vertex_into_container(m_vertex);
			m_graph.insert_edge_into_container(incoming);

			incoming->source()->m_incident[m_index] = incoming;

			// TODO: maintain sorted trait

			if constexpr (G::Graph_traits::decomposed) {
				Path_handle p = outgoing->m_path;
				if (p->m_cyclic && p->m_end == incoming->prev()) {
					p->m_end = incoming;
				} else if (p->m_start == outgoing) {
					p->m_start = incoming;
				}
			}
		}
	}

	void redo() override {
		merge_vertex_no_curve(m_graph, m_vertex, m_index, m_keep_incoming);
		if (m_keep_incoming) {
			std::swap(m_vertex->incoming()->m_representation, m_curve_rep);
		} else {
			std::swap(m_vertex->outgoing()->m_representation, m_curve_rep);
		}
		maintain_traits(m_graph, m_vertex);
	}
};

template <class G> class SplitVertex : public Operation {
  private:
	using Vertex_handle = G::Vertex_handle;
	using Edge_handle = G::Edge_handle;
	using Path_handle = G::Path_handle;
	using Point = G::Point_2;
	using Curve_representation = typename G::Curve_traits::Curve_representation_2;

	G& m_graph;
	// TODO: optimize memory: store only new vertex or new edge
	Vertex_handle m_vertex, m_new_vertex;
	Edge_handle m_new_edge;
	Curve_representation m_curve_rep_inc, m_curve_rep_out;

	static void split_vertex_no_curve(G& g, Vertex_handle v, Vertex_handle new_v, Edge_handle new_e) {
		g.insert_vertex_into_container(new_v);
		g.insert_edge_into_container(new_e);

		v->m_incident[1] = new_e;
		new_v->outgoing()->m_source = new_v;

		if constexpr (G::Graph_traits::decomposed) {

			Path_handle p = new_e->m_path;
			// NB: these cases are not necessarily exclusive, but we need to execute only one
			if (p->m_start == new_v->outgoing()) {
				assert(p->m_cyclic);
				p->m_start = new_e;
			} else if (p->m_end == v->incoming()) {
				assert(p->m_cyclic);
				p->m_end = new_e;
			}
		}
	}

	static void maintain_traits(G& g, Vertex_handle v, Vertex_handle new_v) {

		// TODO: maintain sorted trait
	}

  public:
	SplitVertex(G& g, Vertex_handle v, Vertex_handle new_v, Edge_handle new_e,
	            Curve_representation crep_inc, Curve_representation crep_out)
	    : m_graph(g), m_vertex(v), m_new_vertex(new_v), m_new_edge(new_e),
	      m_curve_rep_inc(crep_inc), m_curve_rep_out(crep_out) {}

	void forget_future() override {
		delete m_new_vertex;
		delete m_new_edge;
	}

	static void split_vertex(G& g, Vertex_handle v, Vertex_handle new_v, Edge_handle new_e,
	                         Curve_representation crep_inc, Curve_representation crep_out) {
		split_vertex_no_curve(g, v, new_v, new_e);
		v->incoming()->m_representation = std::move(crep_inc);
		new_v->outgoing()->m_representation = std::move(crep_inc);
		maintain_traits(g, v, new_v);
	}

	void undo() override {

		m_graph.insert_vertex_into_container(m_new_vertex);
		m_graph.insert_edge_into_container(m_new_edge);

		m_vertex->m_incident[1] = m_new_vertex->outgoing();
		m_new_vertex->outgoing()->m_source = m_vertex;

		// TODO: maintain sorted trait

		if constexpr (G::Graph_traits::decomposed) {

			Path_handle p = m_new_edge->m_path;
			if (p->m_start == m_new_edge) {
				assert(p->m_cyclic);
				p->m_start = m_vertex->outgoing();
			} else if (p->m_end == m_new_edge) {
				assert(p->m_cyclic);
				p->m_end = m_vertex->incoming();
			}
		}
	}
	void redo() override {
		split_vertex_no_curve(m_graph, m_vertex, m_new_vertex, m_new_edge);
		std::swap(m_vertex->incoming()->m_representation, m_curve_rep_inc);
		std::swap(m_new_vertex->outgoing()->m_representation, m_curve_rep_out);
		maintain_traits(m_graph, m_vertex, m_new_vertex);
	}
};

template <class G> class CollapseEdge : public Operation {
  private:
	using Vertex_handle = G::Vertex_handle;
	using Edge_handle = G::Edge_handle;
	using Path_handle = G::Path_handle;
	using Point = G::Point_2;
	using Curve_representation = typename G::Curve_traits::Curve_representation_2;

	G& m_graph;
	Edge_handle m_edge;
	Point m_point;
	Curve_representation m_curve_rep_inc, m_curve_rep_out;

	static void collapse_edge_no_curves(G& g, Edge_handle e) {

		Vertex_handle v = e->source();
		Edge_handle out = e->next();

		// remove
		g.remove_vertex_from_container(e->target());
		g.remove_edge_from_container(e);

		// redirect
		out->m_source = v;
		v->m_incident[1] = out;

		// update paths
		if constexpr (G::Graph_traits::decomposed) {
			Path_handle p = e->m_path;
			if (p->m_start == e) {
				assert(p->m_cyclic);
				p->m_start = out;
			} else if (p->m_end == e) {
				assert(p->m_cyclic);
				p->m_end = v->incoming();
			}
		}
	}

	static void maintain_traits(G& g, Vertex_handle v) {
		// TODO: maintain sorted
	}

  public:
	CollapseEdge(G& g, Edge_handle e, Point p, Curve_representation crep_inc,
	             Curve_representation crep_out)
	    : m_graph(g), m_edge(e), m_point(p), m_curve_rep_inc(crep_inc), m_curve_rep_out(crep_out) {}

	void forget_past() override {
		delete m_edge->target();
		delete m_edge;
	}

	static void collapse_edge(G& g, Edge_handle e, Point p, Curve_representation crep_inc,
	                          Curve_representation crep_out) {
		collapse_edge_no_curves(g, e);
		Vertex_handle v = e->source();
		v->m_point = std::move(p);
		v->incoming()->m_representation = std::move(crep_inc);
		v->outgoing()->m_representation = std::move(crep_out);
		maintain_traits(g, v);
	}

	void undo() override {
		Vertex_handle v = m_edge->source();
		Edge_handle inc = v->incoming();
		Edge_handle out = m_edge->next();

		// remove
		m_graph.insert_vertex_into_container(m_edge->target());
		m_graph.insert_edge_into_container(m_edge);

		// redirect
		out->m_source = m_edge->target();
		v->m_incident[1] = m_edge;

		// update geometry
		std::swap(v->m_point, m_point);
		std::swap(inc->m_representation, m_curve_rep_inc);
		std::swap(out->m_representation, m_curve_rep_out);

		// TODO: maintain sorted

		// update paths
		if constexpr (G::Graph_traits::decomposed) {
			Path_handle p = m_edge->m_path;
			if (p->m_start == out) {
				assert(p->m_cyclic);
				p->m_start = m_edge;
			} else if (p->m_end == inc) {
				assert(p->m_cyclic);
				p->m_end = m_edge;
			}
		}
	}
	void redo() override {
		collapse_edge_no_curves(m_graph, m_edge);
		Vertex_handle v = m_edge->source();
		std::swap(v->m_point, m_point);
		std::swap(v->incoming()->m_representation, m_curve_rep_inc);
		std::swap(v->outgoing()->m_representation, m_curve_rep_out);
		maintain_traits(m_graph, v);
	}
};

template <class G> class SubdivideEdge : public Operation {
  private:
	using Vertex_handle = G::Vertex_handle;
	using Edge_handle = G::Edge_handle;
	using Path_handle = G::Path_handle;
	using Point = G::Point_2;
	using Curve_representation = typename G::Curve_traits::Curve_representation_2;

	G& m_graph;
	Edge_handle m_edge, m_new_edge;
	Curve_representation m_curve_rep;
	size_t m_index;
	bool m_new_out;

	static void subdivide_edge_no_curve(G& g, Edge_handle e, Edge_handle new_e, const size_t index,
	                                    bool newOut) {
		Vertex_handle v = newOut ? new_e->source() : new_e->target();

		g.insert_vertex_into_container(v);
		g.insert_edge_into_container(new_e);

		if (newOut) {
			e->m_target = v;
			new_e->target()->m_incident[index] = new_e;
		} else {
			e->m_source = v;
			new_e->source()->m_incident[index] = new_e;
		}

		if constexpr (G::Graph_traits::decomposed) {

			Path_handle p = new_e->m_path;
			if (newOut) {
				if (p->m_end == e) {
					p->m_end = new_e;
				}
			} else {
				if (p->m_start == e) {
					p->m_start = new_e;
				}
			}
		}
	}

	static void maintain_traits(G& g, Edge_handle e, Edge_handle new_edge) {
		// TODO: maintain sorted
	}

  public:
	SubdivideEdge(G& g, Edge_handle e, Edge_handle new_e, bool newOut, Curve_representation crep)
	    : m_graph(g), m_edge(e), m_new_edge(new_e), m_curve_rep(crep), m_new_out(newOut) {
		m_index =
		    newOut ? m_edge->find_target_incident_index() : m_edge->find_source_incident_index();
	}

	void forget_future() override {
		if (m_new_out)
			delete m_new_edge->source();
		else
			delete m_new_edge->target();
		delete m_new_edge;
	}

	static void subdivide_edge(G& g, Edge_handle e, Edge_handle new_e, bool newOut,
	                           Curve_representation crep) {
		size_t index = newOut ? e->find_target_incident_index() : e->find_source_incident_index();
		subdivide_edge_no_curve(g, e, new_e, index, newOut);
		e->m_representation = std::move(crep);
		maintain_traits(g, e, new_e);
	}

	void undo() override {
		Vertex_handle v = m_new_out ? m_new_edge->source() : m_new_edge->target();

		m_graph.remove_vertex_from_container(v);
		m_graph.remove_edge_from_container(m_new_edge);

		if (m_new_out) {
			m_edge->m_target = m_new_edge->target();
			m_edge->target()->m_incident[m_index] = m_edge;
		} else {
			m_edge->m_source = m_new_edge->source();
			m_edge->source()->m_incident[m_index] = m_edge;
		}

		std::swap(m_edge->m_representation, m_curve_rep);

		// TODO: maintain sorted

		if constexpr (G::Graph_traits::decomposed) {

			Path_handle p = m_new_edge->m_path;
			if (m_new_out) {
				if (p->m_end == m_new_edge) {
					p->m_end = m_edge;
				}
			} else {
				if (p->m_start == m_new_edge) {
					p->m_start = m_edge;
				}
			}
		}
	}
	void redo() override {
		subdivide_edge_no_curve(m_graph, m_edge, m_new_edge, m_index, m_new_out);
		std::swap(m_edge->m_representation, m_curve_rep);
		maintain_traits(m_graph, m_new_edge, m_new_edge);
	}
};
} // namespace detail

} // namespace CGAL::internal