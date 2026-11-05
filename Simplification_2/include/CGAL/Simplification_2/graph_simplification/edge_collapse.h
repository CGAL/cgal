#pragma once

#include <CGAL/Simplification_2/core.h>
#include <CGAL/Simplification_2/data_structures/indexed_priority_queue.h>
#include <CGAL/Simplification_2/data_structures/straight_graph_2.h>
#include <CGAL/Simplification_2/data_structures/graph_traits_2.h>

#include "vertex_quad_tree.h"
#include "edge_quad_tree.h"
#include "common.h"

namespace CGAL::internal {

	namespace detail {

		template<bool H>
		using ECGraphTraits = DecomposedGraph<H, std::monostate>;

		template<typename K, bool H>
		struct ECData;

		template<class G>
		struct ECTraitsBase;

		template <class ECT>
		concept ECTraits = requires(typename ECT::Graph::Edge_handle e) {

			{ ECT::data(e) } -> std::same_as<ECData<typename ECT::Graph::Kernel, ECT::Graph::Graph_traits::historic>&>;
			{ ECT::determine_collapse(e) };

		};
	}

	template<typename K, bool H>
	using EdgeCollapseGraph = Straight_graph_2<std::monostate, detail::ECData<K, H>, K, detail::ECGraphTraits<H>>;

	template <detail::ECTraits ECT> class EdgeCollapse {
	public:
		using Graph = ECT::Graph;
		using Vertex_handle = Graph::Vertex_handle;
		using Edge_handle = Graph::Edge_handle;
		using Kernel = Graph::Kernel;
		using VertexTree = VertexQuadTree<Graph>;
		using EdgeTree = EdgeQuadTree<Graph>;

	private:
		using Queue = CGAL::internal::IndexedPriorityQueue<GraphQueueTraits<Edge_handle, Kernel>>;

		Graph& graph;
		EdgeTree& sqt;
		VertexTree& pqt;
		Queue queue;

		void update(Edge_handle e);
		bool blocks(Edge_handle edge, Edge_handle collapse);
		bool validateState();

		Edge_handle findNextStep();
		void performStep(Edge_handle e);
	public:
		EdgeCollapse(Graph& g, EdgeTree& sqt, VertexTree& pqt);

		void initialize(bool initSQT, bool initPQT);
		bool run(std::optional<std::function<bool(int,Number<Kernel>)>> stop = std::nullopt);
		bool step();

	};


	template <typename G, bool A> struct KronenfeldEtAlTraits : detail::ECTraitsBase<G> {
		using Kernel = G::Kernel;

		static void determine_collapse(typename G::Edge_handle e);
	};

	template <typename G> using KronenfeldEtAl = EdgeCollapse<KronenfeldEtAlTraits<G, false>>;

} // namespace CGAL::internal

#include "edge_collapse.hpp"