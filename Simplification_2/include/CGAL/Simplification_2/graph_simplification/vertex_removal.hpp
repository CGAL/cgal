// -----------------------------------------------------------------------------
// IMPLEMENTATION OF TEMPLATE FUNCTIONS
// Do not include this file, but the .h file instead
// -----------------------------------------------------------------------------

#include "utils.h"

namespace CGAL::internal {

	namespace detail {

		template<typename K, bool H>
		struct VRData {
			Number<K> cost;
			std::vector<typename VertexRemovalGraph<K, H>::Vertex_handle> blocked_by;
			std::vector<typename VertexRemovalGraph<K, H>::Vertex_handle> blocking;
			int queue_index = -1;
		};

		template<class G>
		struct VRTraitsBase {

			using Graph = G;

			static VRData<typename G::Kernel, G::Graph_traits::historic>& data(typename G::Vertex_handle v) {
				return v->data();
			}
		};
	}

	template <detail::VRTraits VRT>
	VertexRemoval<VRT>::VertexRemoval(Graph& g, VertexTree& qt) : graph(g), pqt(qt) {
	}

	template <detail::VRTraits VRT>
	void VertexRemoval<VRT>::initialize(bool initQuadTree) {

		if (initQuadTree) {
			pqt.clear();
			for (Vertex_handle v : graph.vertices()) {
				pqt.insert(v);
			}
		}

		for (Vertex_handle v : graph.vertices()) {
			update(v);
		}
	}

	template <detail::VRTraits VRT>
	bool VertexRemoval<VRT>::run(std::optional<std::function<bool(int, Number<Kernel>)>> stop) {
		while (true) {
			Vertex_handle next = findNextStep();
			if (next == nullptr) {
				return false;
			}

			if (!stop.has_value() || (*stop)(graph.number_of_edges(), VRT::data(next).cost)) {
				return true;
			}

			performStep(next);
		}
	}

	template <detail::VRTraits VRT>
	VRT::Graph::Vertex_handle VertexRemoval<VRT>::findNextStep() {

		assert(graph.can_perform_operation());

		while (!queue.empty()) {
			Vertex_handle v = queue.peek();

			Vertex_handle u = v->prev();
			Vertex_handle w = v->next();

			// test whether the operation is blocked
			const Point<Kernel>& up = u->point();
			const Point<Kernel>& vp = v->point();
			const Point<Kernel>& wp = w->point();
			Triangle<Kernel> T(up, vp, wp);

			Rectangle<Kernel> rect = utils::boxOf(up, vp, wp);

			pqt.findContained(rect, [&T, &u, &v, &w](Vertex_handle b) {
				if (b != u && b != v && b != w && !T.has_on_unbounded_side(b->point())) {
					// blocked, record the pair
					VRT::data(b).blocking.push_back(v);
					VRT::data(v).blocked_by.push_back(b);
				}
				});

			if (VRT::data(v).blocked_by.empty()) {
				// not blocked, this is the next step
				return v;
			}
			else {
				// remove the element from the queue as it's not valid and continue searching
				queue.pop();
			}
		}

		// no valid steps exist
		return nullptr;
	}

	template <detail::VRTraits VRT>
	void VertexRemoval<VRT>::performStep(Vertex_handle v) {

		assert(graph.can_perform_operation());
		assert(queue.peek() == v);

		queue.pop();

		// remove from blocking lists and search structure
		pqt.remove(v);
		for (Vertex_handle b : VRT::data(v).blocking) {

			if (utils::listRemove(v, VRT::data(b).blocked_by)) {
				if (VRT::data(b).blocked_by.empty()) {
					queue.push(b);
				}
			}
		}

		// perform the removal
		Vertex_handle u = v->prev();
		Vertex_handle w = v->next();

		graph.merge_vertex(v);

		// update the neighbors
		update(u);
		update(w);

		// and their common neighbors (to avoid issues with triangles collapsing)
		for (int i = 0; i < u->degree(); i++) {
			Vertex_handle nbr = u->neighbor(i);
			if (nbr->is_neighbor_of(w)) {
				update(nbr);
			}
		}
	}

	template <detail::VRTraits VRT>
	bool VertexRemoval<VRT>::step() {

		assert(graph.can_perform_operation());

		Vertex_handle v = findNextStep();
		if (v == nullptr) {
			return false;
		}

		performStep(v);
		return true;
	}

	template <detail::VRTraits VRT>
	void VertexRemoval<VRT>::update(Vertex_handle v) {
		if (v->degree() != 2) {
			return;
		}

		// clear topology
		for (Vertex_handle b : VRT::data(v).blocked_by) {
			utils::listRemove(v, VRT::data(b).blocking);
		}
		VRT::data(v).blocked_by.clear();

		Vertex_handle u = v->prev();
		Vertex_handle w = v->next();

		if (u->is_neighbor_of(w)) {
			queue.remove(v);
			return;
		}

		VRT::data(v).cost = VRT::compute_cost(v);

		if (queue.contains(v)) {
			queue.update(v);
		}
		else {
			queue.push(v);
		}
	}
}