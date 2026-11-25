// -----------------------------------------------------------------------------
// IMPLEMENTATION OF TEMPLATE FUNCTIONS
// Do not include this file, but the .h file instead
// -----------------------------------------------------------------------------

#include "utils.h"
#include "precision_safe_tests.h"

namespace CGAL::internal {

	namespace detail {
		template <typename K, bool H> struct ECData {

			// collapse specification
			bool erase_both; // special case: both endpoints are to be removed
			Point<K> point; // general case: endpoints merge onto this point
			bool creates_difference; // special case: when the edge is collinear with its neighbors, there are no difference-triangles
			Triangle<K> T1, T2; // the two triangles of difference
			Number<K> cost; // the cost of the collapse

			// algorithm 
			bool blocked_by_degzero;
			std::vector<typename EdgeCollapseGraph<K, H>::Edge_handle> blocked_by;
			std::vector<typename EdgeCollapseGraph<K, H>::Edge_handle> blocking;
			int queue_index = -1;

		};

		template<class G>
		struct ECTraitsBase {
			using Graph = G;

			static ECData<typename G::Kernel, G::Graph_traits::historic>& data(typename G::Edge_handle e) {
				return e->data();
			}
		};
	}

	template <detail::ECTraits ECT>
	void EdgeCollapse<ECT>::update(Edge_handle e) {

		assert(graph.edge(e->graph_index()) == e);

		auto& edata = ECT::data(e);

		// clear topology
		for (Edge_handle b : edata.blocked_by) {
			utils::listRemove(e, ECT::data(b).blocking);
		}
		edata.blocked_by.clear();

		// last condition checks for a triangle
		if (e->source()->degree() != 2 || e->target()->degree() != 2 ||
			e->prev()->source() == e->next()->target()) {
			queue.remove(e);
			return;
		}

		ECT::determine_collapse(e);

		if (queue.contains(e)) {
			queue.update(e);
		}
		else {
			queue.push(e);
		}
	}

	template <detail::ECTraits ECT>
	bool EdgeCollapse<ECT>::blocks(Edge_handle edge, Edge_handle collapse) {
		Edge_handle prev = collapse->prev();
		Edge_handle next = collapse->next();

		if (edge == collapse || edge == prev || edge == next) {
			// involved in collapse
			return false;
		}

		Vertex_handle prev_v = collapse->prev()->source();
		Vertex_handle next_v = collapse->next()->target();;

		bool source_shared = edge->source() == prev_v || edge->source() == next_v;
		bool target_shared = edge->target() == prev_v || edge->target() == next_v;

		auto test_is = [&](std::optional<std::variant<Point<Kernel>, Segment<Kernel>>> is) {
			if (!is.has_value()) {
				// certainly no intersection
				return false;
			}


			if (Point<Kernel>* pt = std::get_if<Point<Kernel>>(&*is)) {
				// make sure it's not the common point
				if (source_shared && safe_test::same_point(*pt, edge->source()->point())) {
					return false;
				}
				if (target_shared && safe_test::same_point(*pt, edge->target()->point())) {
					return false;
				}

				return true;
			}
			else if constexpr (std::is_same<Inexact, Kernel>::value) {

				// running in inexact mode: may need to account for the intersection being a tiny segment near the start vertex

				Segment<Kernel>* ls = std::get_if<Segment<Kernel>>(&*is);

				if (source_shared
					&& safe_test::same_point(ls->source(), edge->source()->point())
					&& safe_test::same_point(ls->target(), edge->source()->point())) {
					return false;
				}

				if (target_shared
					&& safe_test::same_point(ls->source(), edge->target()->point())
					&& safe_test::same_point(ls->target(), edge->target()->point())) {
					return false;
				}

				// intersection doesnt reflect a shared endpoint
				return true;
			}
			else {
				// running in exact mode: 
				// segments always intersect
				return true;
			}
			};

		if (test_is(CGAL::intersection(ECT::data(collapse).T1, edge->curve()))) {
			return true;
		}

		if (test_is(CGAL::intersection(ECT::data(collapse).T2, edge->curve()))) {
			return true;
		}

		return false;
	}

	template <detail::ECTraits ECT>
	EdgeCollapse<ECT>::EdgeCollapse(Graph& g, EdgeTree& sqt, VertexTree& pqt)
		: graph(g), sqt(sqt), pqt(pqt) {

	}

	template <detail::ECTraits ECT>
	void EdgeCollapse<ECT>::initialize(bool initSQT, bool initPQT) {

		if (initSQT) {
			sqt.clear();
			for (Edge_handle e : graph.edges()) {
				sqt.insert(e);
			}
		}

		if (initPQT) {
			pqt.clear();
			for (Vertex_handle v : graph.vertices()) {
				if (v->degree() == 0) {
					pqt.insert(v);
				}
			}
		}

		queue.clear();

		for (Edge_handle e : graph.edges()) {
			auto& edata = ECT::data(e);
			edata.queue_index = -1;
			edata.blocked_by.clear();
			edata.blocking.clear();
			edata.blocked_by_degzero = false;

			update(e);
		}

		assert(validateState());
	}

	template <detail::ECTraits ECT>
	bool EdgeCollapse<ECT>::validateState() {
		bool ok = true;
		for (Edge_handle e : graph.edges()) {
			if (ECT::data(e).queue_index >= 0) {
				if (!queue.contains(e)) {
					std::cout << *e << " :: thinks it's in queue but isn't\n";
					ok = false;
				}
			}
			else if (queue.contains(e)) {
				std::cout << *e << " :: thinks it's not in queue, but is\n";
				ok = false;
			}

			for (Edge_handle b : ECT::data(e).blocked_by) {
				if (std::find(ECT::data(b).blocking.begin(), ECT::data(b).blocking.end(), e) == ECT::data(b).blocking.end()) {
					std::cout << *e << " :: thinks it's blocked by " << *b << ", but they don't agree\n";
					ok = false;
				}
			}

			for (Edge_handle b : ECT::data(e).blocking) {
				if (std::find(ECT::data(b).blocked_by.begin(), ECT::data(b).blocked_by.end(), e) == ECT::data(b).blocked_by.end()) {
					std::cout << *e << " :: thinks it's blocking " << *b << ", but they don't agree\n";
					ok = false;
				}
			}
		}
		for (Edge_handle e : queue.content()) {
			if (graph.edge(e->graph_index()) != e) {
				std::cout << *e << " :: in queue but not in graph!?\n";
				ok = false;
			}
		}
		return ok;
	}

	template <detail::ECTraits ECT>
	bool EdgeCollapse<ECT>::run(std::optional<std::function<bool(int, Number<Kernel>)>> stop) {
		while (true) {
			assert(validateState());

			Edge_handle next = findNextStep();
			if (next == nullptr) {
				return false;
			}

			if (!stop.has_value() || (*stop)(graph.number_of_edges(), ECT::data(next).cost)) {
				return true;
			}

			performStep(next);

			// useful for debugging: stop when an invalid state is reached, without erroring
			//if (!validateState()) {
			//	return false;
			//}
		}
	}

	template <detail::ECTraits ECT>
	ECT::Graph::Edge_handle EdgeCollapse<ECT>::findNextStep() {

		assert(graph.can_perform_operation());

		while (!queue.empty()) {
			Edge_handle e = queue.peek();

			auto& edata = ECT::data(e);

			if (edata.creates_difference) {
				// possibly blocked?

				Rectangle<Kernel> rect = utils::boxOf(edata.T1, edata.T2);

				edata.blocked_by_degzero = false;

				pqt.findContained(rect, [&edata](Vertex_handle b) {
					if (!edata.T1.has_on_unbounded_side(b->point()) ||
						!edata.T2.has_on_unbounded_side(b->point())) {
						// blocked, by an unmovable vertex
						edata.blocked_by_degzero = true;
					}
					});

				if (!edata.blocked_by_degzero) {

					sqt.findOverlapped(rect, [this, &e](Edge_handle b) {

						if (blocks(b, e)) {
							ECT::data(b).blocking.push_back(e);
							ECT::data(e).blocked_by.push_back(b);
						}
						});
				}

			} // else: no difference, cannot be blocked

			if (!edata.blocked_by_degzero && edata.blocked_by.empty()) {
				// not blocked, this is the next step
				return e;
			}
			else {
				// remove the element from the queue as it's not valid and continue searching
				queue.pop();
			}
		}

		// no steps exist
		return nullptr;
	}

	template <detail::ECTraits ECT>
	void EdgeCollapse<ECT>::performStep(Edge_handle e) {

		assert(graph.can_perform_operation());
		assert(queue.peek() == e);

		queue.pop();

		auto& edata = ECT::data(e);

		// remove from blocking lists and search structure
		Edge_handle prev = e->prev();
		Edge_handle next = e->next();
		sqt.remove(e);
		sqt.remove(prev);
		sqt.remove(next);

		queue.remove(prev);
		queue.remove(next);

		// clear blocking data from collapsing edge.
		// NB: this edge will be removed in the collapse, so we dont need to consider blocked_by
		for (Edge_handle b : edata.blocking) {
			auto& bdata = ECT::data(b);
			if (utils::listRemove(e, bdata.blocked_by)) {
				if (bdata.blocked_by.empty() && !bdata.blocked_by_degzero) {
					queue.push(b);
				}
			}
		}
		edata.blocking.clear();

		// clear blocking data from previous edge
		{
			auto& pdata = ECT::data(prev);
			for (Edge_handle b : pdata.blocking) {
				if (b == next) continue; // NB: next could be blocking prev, but will be reset or removed
				auto& bdata = ECT::data(b);
				if (utils::listRemove(prev, bdata.blocked_by)) {
					if (bdata.blocked_by.empty() && !bdata.blocked_by_degzero) {
						queue.push(b);
					}
				}
			}
			pdata.blocking.clear();
			for (Edge_handle b : pdata.blocked_by) {
				utils::listRemove(prev, ECT::data(b).blocking);
			}
			pdata.blocked_by.clear();
		}

		// clear blocking data from next edge
		{
			auto& ndata = ECT::data(next);
			for (Edge_handle b : ndata.blocking) {
				if (b == prev) continue; // NB: prev could be blocking next, but will be reset or removed
				auto& bdata = ECT::data(b);
				if (utils::listRemove(next, bdata.blocked_by)) {
					if (bdata.blocked_by.empty() && !bdata.blocked_by_degzero) {
						queue.push(b);
					}
				}
			}
			ndata.blocking.clear();
			for (Edge_handle b : ndata.blocked_by) {
				utils::listRemove(next, ECT::data(b).blocking);
			}
			ndata.blocked_by.clear();
		}

		Vertex_handle src = e->source();
		Vertex_handle tar = e->target();

		if (edata.erase_both) {

			graph.start_operation_group();

			graph.merge_vertex(src);
			Edge_handle ne = graph.merge_vertex(tar);

			graph.end_operation_group();

			// insert the one new edge
			sqt.insert(ne);

			// update it and its neighbors, if applicable
			update(ne);
			if (ne->source()->degree() == 2) {
				update(ne->prev());
			}
			if (ne->target()->degree() == 2) {
				update(ne->next());
			}
		}
		else {

			// perform the collapse
			Vertex_handle v = graph.collapse_edge(e, edata.point);

			// insert the two new edges
			sqt.insert(v->incoming());
			sqt.insert(v->outgoing());

			// update them and their neighbors, if applicable
			update(v->incoming());
			update(v->outgoing());

			if (v->prev()->degree() == 2) {
				update(v->prev()->incoming());
			}
			if (v->next()->degree() == 2) {
				update(v->next()->outgoing());
			}
		}
	}

	template <detail::ECTraits ECT>
	bool EdgeCollapse<ECT>::step() {

		assert(graph.can_perform_operation());

		Edge_handle e = findNextStep();
		if (e == nullptr) {
			return false;
		}

		performStep(e);
		return true;
	}

	template <typename G, bool A>
	void KronenfeldEtAlTraits<G, A>::determine_collapse(typename G::Edge_handle e) {
		auto& edata = e->data(); // NB: this needs to be the handle, otherwise it doesn't update... (which is fun, because data() already returns a handle...)

		Point<Kernel> a = e->prev()->source()->point();
		Point<Kernel> b = e->source()->point();
		Point<Kernel> c = e->target()->point();
		Point<Kernel> d = e->next()->target()->point();

		bool abc = safe_test::collinear(a, b, c);
		bool bcd = safe_test::collinear(b, c, d);

		if (abc && bcd) {
			edata.erase_both = true;
			edata.creates_difference = false;
			edata.cost = 0;
			return;
		}
		else if (abc) {
			edata.erase_both = false;
			edata.creates_difference = false;
			edata.cost = 0;
			edata.point = c;
			return;
		}
		else if (bcd) {
			edata.erase_both = false;
			edata.creates_difference = false;
			edata.cost = 0;
			edata.point = b;
			return;
		}

		// else, no consecutive collinear edges
		edata.creates_difference = true;

		Polygon<Kernel> P;
		P.push_back(a);
		P.push_back(b);
		P.push_back(c);
		P.push_back(d);

		Line<Kernel> ad(a, d);
		Line<Kernel> ab(a, b);
		Line<Kernel> bc(b, c);
		Line<Kernel> cd(c, d);

		// area = base * height / 2
		// height = 2*area / base
		// so, we're going to rotate the vector d-a, such that we get a normal of length |d-a| = base.
		// To get a vector of length height, we then multiply this vector with height_times_base / base^2.
		// This normalized the vector and makes it length height! (without squareroots...)
		Number<Kernel> height_times_base = 2 * P.area();

		Vector<Kernel> perpv = (d - a).perpendicular(CGAL::CLOCKWISE);

		CGAL::Aff_transformation_2<Kernel> s(CGAL::SCALING, height_times_base / perpv.squared_length());
		perpv = perpv.transform(s);

		CGAL::Aff_transformation_2<Kernel> t(CGAL::TRANSLATION, perpv);
		Line<Kernel> arealine = ad.transform(t);

		bool zero_area = safe_test::point_on_line(arealine.point(), ad);

		if (zero_area) {

			// these should be caught already by the collinearity checks earlier
			assert(!ad.has_on_boundary(b));
			assert(!ad.has_on_boundary(c));

			edata.erase_both = true;

			// implies that neither b nor c is on ad
			auto intersection = CGAL::intersection(bc, ad);
			Point<Kernel> pt = std::get<Point<Kernel>>(*intersection);

			edata.T1 = Triangle<Kernel>(a, b, pt);
			edata.T2 = Triangle<Kernel>(c, d, pt);
		}
		else {
			edata.erase_both = false;

			bool ab_determines_shape;
			// determine type
			if (ad.has_on_positive_side(b) == ad.has_on_positive_side(c)) {
				// same side of ab, so further point determines
				ab_determines_shape = CGAL::squared_distance(b, ad) > CGAL::squared_distance(c, ad);
			}
			else {
				// opposite sides of ad, so the one that is on same side as area line
				ab_determines_shape =
					ad.has_on_positive_side(b) == ad.has_on_positive_side(arealine.point());
			}

			// configure type
			if (ab_determines_shape) {

				auto intersection = CGAL::intersection(arealine, ab);
				edata.point = std::get<Point<Kernel>>(*intersection);

				Line<Kernel> ns = Line<Kernel>(edata.point, d);
				auto intersection2 = CGAL::intersection(bc, ns);
				Point<Kernel> is = std::get<Point<Kernel>>(*intersection2);

				edata.T1 = Triangle<Kernel>(b, is, edata.point);
				edata.T2 = Triangle<Kernel>(c, d, is);
			}
			else {

				auto intersection = CGAL::intersection(arealine, cd);
				edata.point = std::get<Point<Kernel>>(*intersection);

				Line<Kernel> ns = Line<Kernel>(edata.point, a);
				auto intersection2 = CGAL::intersection(bc, ns);
				Point<Kernel> is = std::get<Point<Kernel>>(*intersection2);

				edata.T1 = Triangle<Kernel>(a, b, is);
				edata.T2 = Triangle<Kernel>(c, is, edata.point);
			}

			if constexpr (A) {
				edata.point = convert_kernel<Kernel>(approximate(edata.point));
			}
		}

		// since it is an area preserving method, T1 and T2 have the same area
		edata.cost = 2 * CGAL::abs(edata.T1.area());
	}
} // namespace CGAL::internal