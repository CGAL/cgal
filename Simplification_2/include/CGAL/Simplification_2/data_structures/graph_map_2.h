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
#include "graph_operation_2.h"
#include "graph_traits_2.h"

namespace CGAL::internal {

class Graph_vertex_index_listener {

	template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits, GraphTraits_2 GraphTraits>
	friend class Graph_2;

  protected:
	void clear() {
		resize(0);
	}

	virtual void resize(const size_t size) = 0;
	virtual void add_index() = 0;
	virtual void add_index(const size_t index) = 0;
	virtual void remove_index(const size_t index) = 0;
	virtual void remove_last_index() = 0;
};

class Graph_edge_index_listener {

	template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits, GraphTraits_2 GraphTraits>
	friend class Graph_2;

  protected:
	void clear() {
		resize(0);
	}

	virtual void resize(const size_t size) = 0;
	virtual void add_index() = 0;
	virtual void add_index(const size_t index) = 0;
	virtual void remove_index(const size_t index) = 0;
	virtual void remove_last_index() = 0;
};

class Graph_path_index_listener {

	template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits, GraphTraits_2 GraphTraits>
	friend class Graph_2;

  protected:
	void clear() {
		resize(0);
	}

	virtual void resize(const size_t size) = 0;
	virtual void add_index() = 0;
	virtual void add_index(const size_t index) = 0;
	virtual void remove_index(const size_t index) = 0;
	virtual void remove_last_index() = 0;
};

template <class G, class E, class Base> class Graph_set : public Base {

	friend G;

  private:
	G& m_graph;
	std::vector<size_t> m_index_map;
	std::vector<E> m_contained;

  protected:
	void resize(const size_t size) override {
		assert(size == 0 || size >= m_index_map.size());
		m_index_map.resize(size, 0);
		if (size == 0) {
			m_contained.clear();
		}
	}

	void add_index() override {
		m_index_map.push_back(0);
	}

	void add_index(const size_t index) override {
		m_index_map.push_back(m_index_map[index]);
		m_index_map[index] = 0;
	}

	void remove_index(const size_t index) override {
		const size_t i = m_index_map[index];
		if (i > 0) {
			if (i != m_contained.size()) {
				m_contained[i - 1] = m_contained.back();
				m_index_map[m_contained[i - 1]->graph_index()] = i;
			}
			m_contained.pop_back();
		}
		m_index_map[index] = m_index_map.back();
		m_index_map.pop_back();
	}

	void remove_last_index() override {
		const size_t i = m_index_map.back();
		if (i > 0) {
			if (i != m_contained.size()) {
				m_contained[i - 1] = m_contained.back();
				m_index_map[m_contained[i - 1]->graph_index()] = i;
			}
			m_contained.pop_back();
		}
		m_index_map.pop_back();
	}

  public:
	Graph_set(G& graph) : m_graph(graph) {
		m_graph.add_listener(this);
	}

	~Graph_set() {
		m_graph.remove_listener(this);
	}

	void add(E elt) {
		if (m_index_map[elt->graph_index()] == 0) {
			m_contained.push_back(elt);
			m_index_map[elt->graph_index()] = m_contained.size();
		}
	}
	template <typename Range> void addAll(const Range& elements) {
		for (E elt : elements) {
			add(elt);
		}
	}

	void remove(E elt) {
		const size_t i = m_index_map[elt->graph_index()];
		if (i > 0) {
			if (i != m_contained.size()) {
				m_contained[i - 1] = m_contained.back();
				m_index_map[m_contained[i - 1]->graph_index()] = i;
			}
			m_contained.pop_back();
			m_index_map[elt->graph_index()] = 0;
		}
	}

	void removeAll() {
		for (E elt : m_contained) {
			m_index_map[elt->graph_index()] = 0;
		}
		m_contained.clear();
	}

	bool contains(E elt) const {
		return m_index_map[elt->graph_index()] > 0;
	}

	size_t size() const {
		return m_contained.size();
	}

	E& operator[](const size_t index) {
		return m_contained[index];
	}

	std::vector<E>::iterator begin() {
		return m_contained.begin();
	}

	std::vector<E>::iterator end() {
		return m_contained.end();
	}
};

template <class G>
class Graph_vertex_set
    : public Graph_set<G, typename G::Vertex_handle, Graph_vertex_index_listener> {
  public:
   Graph_vertex_set(G& graph)
	   : Graph_set<G, typename G::Vertex_handle, Graph_vertex_index_listener>(graph) {}
};

template <class G>
class Graph_edge_set : public Graph_set<G, typename G::Edge_handle, Graph_edge_index_listener> {
  public:
	Graph_edge_set(G& graph)
	    : Graph_set<G, typename G::Edge_handle, Graph_edge_index_listener>(graph) {}
};

template <class G>
requires G::Graph_traits::decomposed class Graph_path_set
    : public Graph_set<G, typename G::Path_handle, Graph_path_index_listener> {

  public:
	Graph_path_set(G& graph)
	    : Graph_set<G, typename G::Path_handle, Graph_path_index_listener>(graph) {}
};

template <class G, class E, typename T, class Base> class Graph_map : public Base {

	friend G;

  private:
	G& m_graph;
	std::vector<T> m_vec;
	const T m_init;

  protected:
	void resize(const size_t size) override {
		m_vec.resize(size, m_init);
	}

	void add_index() override {
		m_vec.push_back(m_init);
	}

	void add_index(const size_t index) override {
		m_vec.push_back(m_vec[index]);
		m_vec[index] = m_init;
	}

	void remove_index(const size_t index) override {
		m_vec[index] = m_vec[m_vec.size() - 1];
		m_vec.pop_back();
	}

	void remove_last_index() override {
		m_vec.pop_back();
	}

  public:
	Graph_map(G& graph, const T init) : m_graph(graph), m_init(init) {
		m_graph.add_listener(this);
	}

	~Graph_map() {
		m_graph.remove_listener(this);
	}

	T& operator[](const E elt) {
		return m_vec[elt->graph_index()];
	}

	void assign(const T v) {
		m_vec.assign(m_vec.size(), v);
	}
};

template <class G, typename T>
class Graph_vertex_map : public Graph_map<G, typename G::Vertex_handle, T, Graph_vertex_index_listener> {	
	public:
	Graph_vertex_map(G& graph, const T init)
	    : Graph_map<G, typename G::Vertex_handle, T, Graph_vertex_index_listener>(graph, init) {}
};

template <class G, typename T> class Graph_static_vertex_map {

  public:
	using Vertex_handle = G::Vertex_handle;

  private:
	std::vector<T> m_vec;

  public:
	Graph_static_vertex_map(const G& graph) {
		m_vec.resize(graph.number_of_vertices());
	}

	Graph_static_vertex_map(const G& graph, const T init) {
		m_vec.resize(graph.number_of_vertices(), init);
	}

	std::vector<T>::reference operator[](const Vertex_handle vtx) {
		return m_vec[vtx->graph_index()];
	}

	void assign(const T v) {
		m_vec.assign(m_vec.size(), v);
	}

	void resize(size_t size) {
		m_vec.resize(size);
	}

	void resize(size_t size, const T init) {
		m_vec.resize(size, init);
	}

	size_t size() {
		return m_vec.size();
	}
};

template <class G> class Graph_vertex_index_map {
  public:
	using Vertex_handle = G::Vertex_handle;

	size_t operator[](const Vertex_handle vtx) {
		return vtx->graph_index();
	}
};

template <class G, typename T>
class Graph_edge_map : public Graph_map<G, typename G::Edge_handle, T, Graph_edge_index_listener> {

  public:
	Graph_edge_map(G& graph, const T init)
	    : Graph_map<G, typename G::Edge_handle, T, Graph_edge_index_listener>(graph, init) {
	}
};
	

template <class G, typename T> class Graph_static_edge_map {

  public:
	using Edge_handle = G::Edge_handle;

  private:
	std::vector<T> m_vec;

  public:
	Graph_static_edge_map(const G& graph) {
		m_vec.resize(graph.number_of_edges());
	}

	Graph_static_edge_map(const G& graph, const T init) {
		m_vec.resize(graph.number_of_edges(), init);
	}

	std::vector<T>::reference operator[](const Edge_handle edge) {
		return m_vec[edge->graph_index()];
	}

	void assign(const T v) {
		m_vec.assign(m_vec.size(), v);
	}

	void resize(size_t size) {
		m_vec.resize(size);
	}

	void resize(size_t size, const T init) {
		m_vec.resize(size, init);
	}

	size_t size() {
		return m_vec.size();
	}
};

template <class G> class Graph_edge_index_map {
  public:
	using Edge_handle = G::Edge_handle;

	size_t operator[](const Edge_handle edge) {
		return edge->graph_index();
	}
};

template <class G, typename T>
class Graph_path_map : public Graph_map<G, typename G::Path_handle, T, Graph_path_index_listener> {	
	
  public:
	Graph_path_map(G& graph, const T init)
	    : Graph_map<G, typename G::Path_handle, T, Graph_path_index_listener>(graph, init) {}
};

template <class G, typename T> class Graph_static_path_map {

  public:
	using Path_handle = G::Path_handle;

  private:
	std::vector<T> m_vec;

  public:
	Graph_static_path_map(const G& graph) {
		m_vec.resize(graph.number_of_edges());
	}

	Graph_static_path_map(const G& graph, const T init) {
		m_vec.resize(graph.number_of_edges(), init);
	}

	std::vector<T>::reference operator[](const Path_handle edge) {
		return m_vec[edge->graph_index()];
	}

	void assign(const T v) {
		m_vec.assign(m_vec.size(), v);
	}

	void resize(size_t size) {
		m_vec.resize(size);
	}

	void resize(size_t size, const T init) {
		m_vec.resize(size, init);
	}

	size_t size() {
		return m_vec.size();
	}
};

template <class G> class Graph_path_index_map {
  public:
	using Path_handle = G::Path_handle;

	size_t operator[](const Path_handle path) {
		return path->graph_index();
	}
};

} // namespace CGAL::internal