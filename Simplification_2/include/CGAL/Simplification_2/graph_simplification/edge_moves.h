#pragma once

#include <CGAL/Simplification_2/core.h>
#include <CGAL/Simplification_2/data_structures/indexed_priority_queue.h>
#include <CGAL/Simplification_2/data_structures/straight_graph_2.h>
#include <CGAL/Simplification_2/data_structures/graph_traits_2.h>

#include "vertex_quad_tree.h"
#include "edge_quad_tree.h"

namespace CGAL::internal {

	namespace detail {

		template<class G> class EdgeSet;

		template<class G> struct BaseMove;
		template<class G> struct SingleMove;
		template<class G> struct ComboMove;

		template <class K, bool H> struct EMData;
		template <class K, bool H> struct EMPathData;

		template<class K, bool H>
		using EMGraphTraits = CustomGraphTraits<H, true, true, true, EMPathData<K, H>>;

		template <class EMT>
		concept EMTraits = requires(typename EMT::Graph::Edge_handle e, typename EMT::Graph::Path_handle p, SingleMove<typename EMT::Graph>&sm, ComboMove<typename EMT::Graph>&cm) {

			typename EMT::Graph;

			{
				EMT::data(e)
			} -> std::same_as<EMData<typename EMT::Kernel, EMT::Graph::Graph_traits::historic>&>;

			{
				EMT::data(p)
			} -> std::same_as<EMPathData<typename EMT::Kernel, EMT::Graph::Graph_traits::historic>&>;

			{ EMT::determineSingleCost(sm) } -> std::same_as< Number<typename EMT::Kernel>>;
			{ EMT::determineComboCost(cm) } -> std::same_as< Number<typename EMT::Kernel>>;
		};

		template<class G>
		struct MoveQueueTraits;
	}

	template<typename K, bool H>
	using EdgeMovesGraph = Straight_graph_2<std::monostate, detail::EMData<K, H>, K, detail::EMGraphTraits<K, H>>;

	template <detail::EMTraits EMT> class EdgeMoves {
	public:
		using Graph = EMT::Graph;
		using Vertex_handle = Graph::Vertex_handle;
		using Edge_handle = Graph::Edge_handle;
		using Kernel = Graph::Kernel;
		using VertexTree = VertexQuadTree<Graph>;
		using EdgeTree = EdgeQuadTree<Graph>;

	private:
		using Data = detail::EMData<Kernel, Graph::Graph_traits::historic>;
		using Move = detail::BaseMove<Graph>;
		using Single = detail::SingleMove<Graph>;
		using Combo = detail::ComboMove<Graph>;
		using Queue = CGAL::internal::IndexedPriorityQueue<detail::MoveQueueTraits<Graph>>;
		
		using PairedSingles = std::pair<Single*, Single*>;
		using Operation = std::variant<Combo*, PairedSingles, Single*>;

		Graph& graph;
		EdgeTree& sqt;
		VertexTree& pqt;
		Queue queue;
		Graph_edge_set<Graph> edgeset;

		bool validate_state();

		template<bool LEFT>
		void update_single(Edge_handle e, Single& single);
		void update_singles(Edge_handle e);
		void update_combo(Edge_handle e);
		bool test_topology(Single& single);
		bool test_topology(Combo& combo);
		Single* find_compensate_move(Single& contract, Number<Kernel> area);

		std::optional<Operation> findNextStep();
		void start_move(Single& move);
		void contract(Single& move);
		void move(Single& move, const Number<Kernel> area);
		void move(Single& move, const Point<Kernel> source_dest, const Point<Kernel> target_dest);
		void performStep(Single& contract, Single& compensate);
		void performStep(Combo& combo);
		void performStep(Single& contact);

		void determineCheckout(Single& move, const bool contract);
		void checkOut();
		void postProcess();

	public:
		EdgeMoves(Graph& g, EdgeTree& sqt, VertexTree& pqt);

		void initialize(bool initSQT, bool initPQT);
		bool run(std::optional<std::function<bool(int, Number<Kernel>)>> stop = std::nullopt);
		bool step();
	};


	template <typename G> struct BuchinEtAlTraits {
		using Graph = G;
		using Kernel = Graph::Kernel;

		static detail::EMData<Kernel, Graph::Graph_traits::historic>& data(typename Graph::Edge_handle e) {
			return e->data();
		}
		static detail::EMPathData<Kernel, Graph::Graph_traits::historic>& data(typename Graph::Path_handle p) {
			return p->data();
		}

		static Number<Kernel> determineSingleCost(detail::SingleMove<Graph>& sm);
		static Number<Kernel> determineComboCost(detail::ComboMove<Graph>& cm);
	};

	template <typename G> using BuchinEtAl = EdgeMoves<BuchinEtAlTraits<G>>;

} // namespace CGAL::internal

#include "edge_moves.hpp"