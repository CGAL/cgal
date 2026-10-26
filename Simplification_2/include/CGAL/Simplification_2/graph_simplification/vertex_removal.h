#pragma once

#include <CGAL/Simplification_2/core.h>
#include <CGAL/Simplification_2/data_structures/indexed_priority_queue.h>
#include <CGAL/Simplification_2/data_structures/straight_graph_2.h>
#include <CGAL/Simplification_2/data_structures/graph_traits_2.h>

#include "vertex_quad_tree.h"
#include "common.h"

namespace CGAL::internal {

	namespace detail {

		template<bool H>
		using VRGraphTraits = DecomposedGraph<H, std::monostate>;

		template<typename K, bool H>
		struct VRData;

		template<class G>
		struct VRTraitsBase;

		template<class VRT>
		concept VRTraits = requires(typename VRT::Graph::Vertex_handle v) {

			typename VRT::Graph;

			{ VRT::data(v) } -> std::same_as<VRData<typename VRT::Graph::Kernel, VRT::Graph::Graph_traits::historic>&>;
			{ VRT::compute_cost(v) } -> std::same_as<Number<typename VRT::Graph::Kernel>>;
		};
	}

	template<typename K, bool H>
	using VertexRemovalGraph = Straight_graph_2<detail::VRData<K, H>, std::monostate, K, detail::VRGraphTraits<H>>;

	/// <summary>
	/// The Vertex Removal algorithm. It is topologically safe, ensuring that vertices are only erased if they have degree 2 and the triangle spanned with its neighbors is empty. Can be configured with custom cost function, via the VertexRemovalTraits.
	/// </summary>
	/// <typeparam name="MG">Modifiable Graph type to be used; typically, will be one of VertexRemovalGraph or HistoricVertexRemovalGraph</typeparam>
	/// <typeparam name="VRT">VertexRemovalTraits, specifying the desired cost function</typeparam>
	template <detail::VRTraits VRT> class VertexRemoval {

	public:
		using Graph = VRT::Graph;
		using Vertex_handle = Graph::Vertex_handle;
		using Kernel = Graph::Kernel;
		using VertexTree = VertexQuadTree<Graph>;

	private:
		using Queue = CGAL::internal::IndexedPriorityQueue<GraphQueueTraits<Vertex_handle, Kernel>>;

		Graph& graph;
		VertexTree& pqt;
		Queue queue;

		void update(Vertex_handle v);

		Vertex_handle findNextStep();
		void performStep(Vertex_handle v);
	public:
		VertexRemoval(Graph& g, VertexTree& qt);

		void initialize(bool initQuadTree);
		bool run(std::optional<std::function<bool(int, Number<Kernel>)>> stop = std::nullopt);
		bool step();
	};


	/// <summary>
	/// Traits for running VisvalingamWhyatt vertex-removal algorithms. The cost is equal to the area of the spanned triangle.
	/// </summary>
	/// <typeparam name="G">The graph type for the algorithm</typeparam>
	template <class G>
	struct VisvalingamWhyattTraits : detail::VRTraitsBase<G> {
		using Kernel = G::Kernel;

		static Number<Kernel> compute_cost(typename G::Vertex_handle v) {
			return CGAL::abs(
				CGAL::area(v->point(), v->prev()->point(), v->next()->point()));
		}
	};

	/// <summary>
	/// Shorthand for the VisvalingamWhyatt vertex-removal algorithm.
	/// </summary>
	/// <typeparam name="G">The graph type for the algorithm</typeparam>
	template <class G>
	using VisvalingamWhyatt = VertexRemoval<VisvalingamWhyattTraits<G>>;

} // namespace CGAL::internal

#include "vertex_removal.hpp"